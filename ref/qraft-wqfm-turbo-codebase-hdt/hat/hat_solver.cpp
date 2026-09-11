#include "hat/hat_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>

// malloc_trim controls only glibc's allocator; jemalloc releases memory
// through its own decay/purge policy.
#if defined(__GLIBC__) && !defined(WQFM_LINKED_JEMALLOC)
#include <malloc.h>
#endif

namespace wqfm::internal {

void rethrow_hat_parallel_child_failure(
    const std::exception_ptr& first,
    const std::exception_ptr& second) {
  if (first) std::rethrow_exception(first);
  if (second) std::rethrow_exception(second);
}

namespace {

constexpr int kHatMinSide = 1;
constexpr int kHatParallelSubproblemMinIdentities = 64;

struct HatIdentityPayload {
  long long identities = 0;
  long long support_words = 0;
  long long taxa = 0;
  long long key_bytes = 0;

  HatIdentityPayload& operator+=(const HatIdentityPayload& other) {
    identities += other.identities;
    support_words += other.support_words;
    taxa += other.taxa;
    key_bytes += other.key_bytes;
    return *this;
  }
};

HatIdentityPayload hat_identity_payload(const Identity& identity) {
  return HatIdentityPayload{
      1,
      static_cast<long long>(identity.support.words().size()),
      static_cast<long long>(identity.support_taxa.size()),
      static_cast<long long>(identity.stable_key.size()),
  };
}

HatIdentityPayload hat_identity_payload(
    const std::vector<Identity>& identities) {
  HatIdentityPayload payload;
  for (const Identity& identity : identities) {
    payload += hat_identity_payload(identity);
  }
  return payload;
}

class HatPayloadLease {
 public:
  HatPayloadLease() = default;
  HatPayloadLease(const HatPayloadLease&) = delete;
  HatPayloadLease& operator=(const HatPayloadLease&) = delete;

  HatPayloadLease(HatPayloadLease&& other) noexcept {
    move_from(other);
  }

  HatPayloadLease& operator=(HatPayloadLease&& other) noexcept {
    if (this != &other) {
      release();
      move_from(other);
    }
    return *this;
  }

  ~HatPayloadLease() { release(); }

  void acquire(bool enabled,
               std::atomic<long long>& current,
               std::atomic<long long>& peak,
               long long amount) {
    release();
    if (!enabled || amount <= 0) return;
    current_ = &current;
    amount_ = amount;
    benchmark_live_payload_acquire(current, peak, amount);
  }

  void release() {
    if (current_ == nullptr) return;
    benchmark_live_payload_release(*current_, amount_);
    current_ = nullptr;
    amount_ = 0;
  }

 private:
  void move_from(HatPayloadLease& other) noexcept {
    current_ = other.current_;
    amount_ = other.amount_;
    other.current_ = nullptr;
    other.amount_ = 0;
  }

  std::atomic<long long>* current_ = nullptr;
  long long amount_ = 0;
};

class HatIdentityPayloadLease {
 public:
  HatIdentityPayloadLease() = default;
  HatIdentityPayloadLease(bool enabled, HatIdentityPayload payload)
      : enabled_(enabled), payload_(payload) {
    acquire_payload();
  }
  HatIdentityPayloadLease(const HatIdentityPayloadLease&) = delete;
  HatIdentityPayloadLease& operator=(const HatIdentityPayloadLease&) = delete;

  HatIdentityPayloadLease(HatIdentityPayloadLease&& other) noexcept {
    move_from(other);
  }

  HatIdentityPayloadLease& operator=(HatIdentityPayloadLease&& other) noexcept {
    if (this != &other) {
      release();
      move_from(other);
    }
    return *this;
  }

  ~HatIdentityPayloadLease() { release(); }

  const HatIdentityPayload& payload() const { return payload_; }

  void mark_continuation() {
    if (!enabled_ || continuation_) return;
    BenchmarkStats& stats = benchmark_stats();
    benchmark_live_gauge_acquire(
        stats.hat_outer_live_continuation_identity_entries,
        stats.hat_outer_peak_continuation_identity_entries,
        payload_.identities);
    benchmark_live_gauge_acquire(
        stats.hat_outer_live_continuation_support_words,
        stats.hat_outer_peak_continuation_support_words,
        payload_.support_words);
    benchmark_live_gauge_acquire(
        stats.hat_outer_live_continuation_taxon_entries,
        stats.hat_outer_peak_continuation_taxon_entries, payload_.taxa);
    benchmark_live_gauge_acquire(
        stats.hat_outer_live_continuation_key_bytes,
        stats.hat_outer_peak_continuation_key_bytes, payload_.key_bytes);
    continuation_ = true;
  }

  void unmark_continuation() {
    if (!continuation_) return;
    BenchmarkStats& stats = benchmark_stats();
    benchmark_live_gauge_release(
        stats.hat_outer_live_continuation_identity_entries,
        payload_.identities);
    benchmark_live_gauge_release(
        stats.hat_outer_live_continuation_support_words,
        payload_.support_words);
    benchmark_live_gauge_release(
        stats.hat_outer_live_continuation_taxon_entries, payload_.taxa);
    benchmark_live_gauge_release(
        stats.hat_outer_live_continuation_key_bytes, payload_.key_bytes);
    continuation_ = false;
  }

  void release() {
    if (enabled_) {
      unmark_continuation();
      BenchmarkStats& stats = benchmark_stats();
      benchmark_live_payload_release(stats.hat_outer_live_identity_key_bytes,
                                     payload_.key_bytes);
      benchmark_live_payload_release(
          stats.hat_outer_live_identity_taxon_entries, payload_.taxa);
      benchmark_live_payload_release(stats.hat_outer_live_support_words,
                                     payload_.support_words);
      benchmark_live_payload_release(stats.hat_outer_live_identity_entries,
                                     payload_.identities);
    }
    enabled_ = false;
    payload_ = HatIdentityPayload{};
  }

 private:
  void acquire_payload() {
    if (!enabled_) return;
    BenchmarkStats& stats = benchmark_stats();
    benchmark_live_payload_acquire(stats.hat_outer_live_identity_entries,
                                   stats.hat_outer_peak_identity_entries,
                                   payload_.identities);
    benchmark_live_payload_acquire(stats.hat_outer_live_support_words,
                                   stats.hat_outer_peak_support_words,
                                   payload_.support_words);
    benchmark_live_payload_acquire(
        stats.hat_outer_live_identity_taxon_entries,
        stats.hat_outer_peak_identity_taxon_entries, payload_.taxa);
    benchmark_live_payload_acquire(stats.hat_outer_live_identity_key_bytes,
                                   stats.hat_outer_peak_identity_key_bytes,
                                   payload_.key_bytes);
  }

  void move_from(HatIdentityPayloadLease& other) noexcept {
    enabled_ = other.enabled_;
    continuation_ = other.continuation_;
    payload_ = other.payload_;
    other.enabled_ = false;
    other.continuation_ = false;
    other.payload_ = HatIdentityPayload{};
  }

  bool enabled_ = false;
  bool continuation_ = false;
  HatIdentityPayload payload_;
};

struct HatTreePayload {
  long long nodes = 0;
  long long adjacencies = 0;
};

HatTreePayload hat_tree_payload(const SolverTree& tree) {
  HatTreePayload payload;
  payload.nodes = static_cast<long long>(tree.nodes.size());
  for (const SolverTree::Node& node : tree.nodes) {
    payload.adjacencies += static_cast<long long>(node.neighbors.size());
  }
  return payload;
}

HatTreePayload accounted_hat_tree_payload(const SolverTree& tree,
                                           bool accounting) {
  return accounting ? hat_tree_payload(tree) : HatTreePayload{};
}

HatTreePayload accounted_hat_glue_scratch_payload(
    const SolverTree& copied_tree,
    const AttachmentEdge& left_edge,
    const AttachmentEdge& right_edge,
    bool accounting) {
  HatTreePayload payload =
      accounted_hat_tree_payload(copied_tree, accounting);
  if (!accounting) return payload;
  const auto note_endpoint = [&](const AttachmentEdge& edge) {
    if (edge.a < 0) {
      // Empty-tree attachment materializes one isolated connector.
      ++payload.nodes;
    } else if (edge.b >= 0) {
      // Subdividing an ordinary edge adds one node and replaces one edge by
      // two, for a net increase of two directed adjacencies.
      ++payload.nodes;
      payload.adjacencies += 2;
    }
    // A singleton endpoint reuses its sole vertex and adds no local payload.
  };
  note_endpoint(left_edge);
  note_endpoint(right_edge);
  // Joining the two attachment locations adds one undirected edge.
  payload.adjacencies += 2;
  return payload;
}

class HatCompletedTreeLease {
 public:
  HatCompletedTreeLease() = default;
  HatCompletedTreeLease(bool enabled, HatTreePayload payload) {
    acquire(enabled, payload);
  }

  void acquire(bool enabled, HatTreePayload payload) {
    release();
    BenchmarkStats& stats = benchmark_stats();
    nodes_.acquire(enabled, stats.hat_outer_live_completed_tree_nodes,
                   stats.hat_outer_peak_completed_tree_nodes, payload.nodes);
    adjacencies_.acquire(
        enabled, stats.hat_outer_live_completed_tree_adjacencies,
        stats.hat_outer_peak_completed_tree_adjacencies,
        payload.adjacencies);
  }

  void release() {
    adjacencies_.release();
    nodes_.release();
  }

 private:
  HatPayloadLease nodes_;
  HatPayloadLease adjacencies_;
};

struct HatBankPayload {
  long long nodes = 0;
  long long other_entries = 0;
};

HatBankPayload hat_bank_payload(const std::vector<GeneIndex>& genes) {
  HatBankPayload payload;
  payload.other_entries = static_cast<long long>(genes.size());
  auto traversal_entries = [](const GeneTraversal& traversal) {
    return static_cast<long long>(
        traversal.parent.size() + traversal.preorder.size() +
        traversal.postorder.size() + traversal.tin.size() +
        traversal.tout.size());
  };
  for (const GeneIndex& gene : genes) {
    payload.nodes += static_cast<long long>(gene.nodes.size());
    payload.other_entries +=
        static_cast<long long>(gene.taxon_to_node.size()) +
        traversal_entries(gene.root_traversal) +
        traversal_entries(gene.reference_traversal);
    for (const GeneIndex::Node& node : gene.nodes) {
      payload.other_entries += static_cast<long long>(node.neighbors.size());
    }
  }
  return payload;
}

HatBankPayload accounted_hat_bank_payload(
    const std::vector<GeneIndex>& genes,
    bool accounting) {
  return accounting ? hat_bank_payload(genes) : HatBankPayload{};
}

class HatBankPayloadLease {
 public:
  HatBankPayloadLease() = default;
  HatBankPayloadLease(bool enabled, HatBankPayload payload) {
    acquire(enabled, payload);
  }

  void acquire(bool enabled, HatBankPayload payload) {
    release();
    BenchmarkStats& stats = benchmark_stats();
    nodes_.acquire(enabled, stats.hat_outer_live_bank_nodes,
                   stats.hat_outer_peak_bank_nodes, payload.nodes);
    other_.acquire(enabled, stats.hat_outer_live_bank_other_entries,
                   stats.hat_outer_peak_bank_other_entries,
                   payload.other_entries);
  }

  void release() {
    other_.release();
    nodes_.release();
  }

 private:
  HatPayloadLease nodes_;
  HatPayloadLease other_;
};

struct HatRecursionControl {
  bool accounting = false;
  bool force_shape = false;
  HatRecursionTestShape shape = HatRecursionTestShape::Balanced;
};

struct HatIdentityTaxonCatalog {
  struct Entry {
    int64_t identity_id = -1;
    std::size_t taxon_begin = 0;
    std::size_t taxon_count = 0;
  };

  std::vector<Entry> entries;
  std::vector<int> taxa;
  std::vector<std::pair<int64_t, std::size_t>> entry_by_identity;
};

HatIdentityTaxonCatalog make_hat_identity_taxon_catalog(
    const std::vector<Identity>& identities,
    bool accounting) {
  HatIdentityTaxonCatalog catalog;
  std::size_t taxon_count = 0;
  for (const Identity& identity : identities) {
    taxon_count += identity.support_taxa.size();
  }
  catalog.entries.reserve(identities.size());
  catalog.taxa.reserve(taxon_count);
  catalog.entry_by_identity.reserve(identities.size());
  for (const Identity& identity : identities) {
    const std::size_t begin = catalog.taxa.size();
    catalog.taxa.insert(catalog.taxa.end(), identity.support_taxa.begin(),
                        identity.support_taxa.end());
    const std::size_t entry = catalog.entries.size();
    catalog.entries.push_back(HatIdentityTaxonCatalog::Entry{
        identity.id, begin, identity.support_taxa.size()});
    catalog.entry_by_identity.emplace_back(identity.id, entry);
  }
  std::sort(catalog.entry_by_identity.begin(),
            catalog.entry_by_identity.end());
  catalog.entry_by_identity.erase(
      std::unique(catalog.entry_by_identity.begin(),
                  catalog.entry_by_identity.end(),
                  [](const auto& left, const auto& right) {
                    return left.first == right.first;
                  }),
      catalog.entry_by_identity.end());
  if (accounting) {
    BenchmarkStats& stats = benchmark_stats();
    benchmark_counter_add(stats.hat_outer_catalog_identity_entries,
                          static_cast<long long>(catalog.entries.size()));
    benchmark_counter_add(stats.hat_outer_catalog_taxon_entries,
                          static_cast<long long>(catalog.taxa.size()));
    benchmark_counter_add(stats.hat_outer_catalog_lookup_entries,
                          static_cast<long long>(
                              catalog.entry_by_identity.size()));
  }
  return catalog;
}

const HatIdentityTaxonCatalog::Entry* hat_catalog_entry(
    const HatIdentityTaxonCatalog& catalog,
    int64_t identity_id) {
  const auto it = std::lower_bound(
      catalog.entry_by_identity.begin(), catalog.entry_by_identity.end(),
      identity_id, [](const auto& entry, int64_t sought) {
        return entry.first < sought;
      });
  return it == catalog.entry_by_identity.end() || it->first != identity_id
             ? nullptr
             : &catalog.entries[it->second];
}

struct HatTrackedTaxa {
  std::vector<int> values;
  HatPayloadLease lease;
};

void radix_sort_hat_taxa(std::vector<int>& taxa, bool accounting) {
  if (taxa.size() <= 1) return;
  BenchmarkStats& stats = benchmark_stats();
  HatPayloadLease sort_payload;
  sort_payload.acquire(accounting, stats.hat_outer_live_taxon_scratch_entries,
                       stats.hat_outer_peak_taxon_scratch_entries,
                       static_cast<long long>(taxa.size()) * 2);
  std::vector<int> scratch(taxa.size());
  std::vector<int>* source = &taxa;
  std::vector<int>* destination = &scratch;
  for (int shift = 0; shift < 32; shift += 8) {
    std::array<std::size_t, 256> counts{};
    for (int taxon : *source) {
      ++counts[(static_cast<unsigned int>(taxon) >> shift) & 0xffU];
    }
    std::size_t offset = 0;
    for (std::size_t& count : counts) {
      const std::size_t bucket_size = count;
      count = offset;
      offset += bucket_size;
    }
    for (int taxon : *source) {
      (*destination)[counts[(static_cast<unsigned int>(taxon) >> shift) &
                            0xffU]++] = taxon;
    }
    std::swap(source, destination);
  }
  if (source != &taxa) taxa.swap(*source);
  if (accounting) {
    benchmark_counter_add(stats.hat_outer_taxon_radix_entries,
                          static_cast<long long>(taxa.size()) * 8);
    benchmark_counter_add(stats.hat_outer_taxon_radix_bucket_ops, 4 * 256);
  }
}

HatTrackedTaxa hat_taxa_from_solver_tree(
    const SolverTree& tree,
    const HatIdentityTaxonCatalog& catalog,
    bool accounting) {
  HatTrackedTaxa tracked;
  std::vector<int>& taxa = tracked.values;
  long long identity_lookups = 0;
  for (const SolverTree::Node& node : tree.nodes) {
    if (!node.active || node.identity_id < 0) continue;
    ++identity_lookups;
    const HatIdentityTaxonCatalog::Entry* entry =
        hat_catalog_entry(catalog, node.identity_id);
    if (entry == nullptr) continue;
    const auto begin = catalog.taxa.begin() +
                       static_cast<std::ptrdiff_t>(entry->taxon_begin);
    taxa.insert(taxa.end(), begin,
                begin + static_cast<std::ptrdiff_t>(entry->taxon_count));
  }
  radix_sort_hat_taxa(taxa, accounting);
  taxa.erase(std::unique(taxa.begin(), taxa.end()), taxa.end());
  if (accounting) {
    BenchmarkStats& stats = benchmark_stats();
    benchmark_counter_add(stats.hat_outer_tree_node_visits,
                          static_cast<long long>(tree.nodes.size()));
    benchmark_counter_add(stats.hat_outer_tree_identity_lookups,
                          identity_lookups);
    benchmark_counter_add(stats.hat_outer_tree_taxon_entries,
                          static_cast<long long>(taxa.size()));
    benchmark_counter_add(stats.hat_outer_taxon_unique_entries,
                          static_cast<long long>(taxa.size()));
    tracked.lease.acquire(
        true, stats.hat_outer_live_taxon_scratch_entries,
        stats.hat_outer_peak_taxon_scratch_entries,
        static_cast<long long>(taxa.size()));
  }
  return tracked;
}

void trim_hat_allocator_after_large_release(int active_identities) {
#if defined(__GLIBC__) && !defined(WQFM_LINKED_JEMALLOC)
  if (active_identities >= 128) {
    malloc_trim(0);
  }
#else
  (void)active_identities;
#endif
}

uint64_t hat_child_path(uint64_t path, int child_offset, int depth) {
  return stable_hash_u64(path, static_cast<uint64_t>(child_offset),
                         static_cast<uint64_t>(depth));
}

std::vector<int> singleton_taxa_from_identities(
    const std::vector<Identity>& identities) {
  std::vector<int> taxa(identities.size(), -1);
  for (int i = 0; i < static_cast<int>(identities.size()); ++i) {
    if (identities[static_cast<std::size_t>(i)].support_taxa.size() == 1) {
      taxa[static_cast<std::size_t>(i)] =
          identities[static_cast<std::size_t>(i)].support_taxa.front();
    }
  }
  return taxa;
}

template <typename T>
std::shared_ptr<const std::vector<T>> track_hat_identity_metadata(
    std::unique_ptr<std::vector<T>> values,
    bool accounting) {
  const long long entries = static_cast<long long>(values->size());
  BenchmarkStats& stats = benchmark_stats();
  if (accounting) {
    benchmark_live_payload_acquire(
        stats.hat_outer_live_bank_identity_entries,
        stats.hat_outer_peak_bank_identity_entries, entries);
  }
  return std::shared_ptr<const std::vector<T>>(
      values.release(), [accounting, entries](const std::vector<T>* ptr) {
        delete ptr;
        if (accounting) {
          benchmark_live_payload_release(
              benchmark_stats().hat_outer_live_bank_identity_entries,
              entries);
        }
      });
}

std::shared_ptr<const std::vector<int64_t>> hat_identity_ids(
    const std::vector<Identity>& identities,
    bool accounting) {
  auto ids = std::make_unique<std::vector<int64_t>>();
  ids->reserve(identities.size());
  for (const Identity& identity : identities) ids->push_back(identity.id);
  return track_hat_identity_metadata(std::move(ids), accounting);
}

std::shared_ptr<const std::vector<int>> hat_identity_taxa(
    const std::vector<Identity>& identities,
    bool accounting) {
  return track_hat_identity_metadata(
      std::make_unique<std::vector<int>>(
          singleton_taxa_from_identities(identities)),
      accounting);
}

int hat_source_identity_count(const SolverContext& ctx,
                              int active_identity_count) {
  if (ctx.active_gene_source_identity_count > 0) {
    return ctx.active_gene_source_identity_count;
  }
  if (ctx.active_gene_identity_ids != nullptr) {
    return static_cast<int>(ctx.active_gene_identity_ids->size());
  }
  if (ctx.owned_active_genes != nullptr || ctx.active_genes != nullptr) {
    return active_identity_count;
  }
  return static_cast<int>(context_taxa(ctx).size());
}

void note_hat_scoring_source(const SolverContext& ctx,
                             int active_identities) {
  const int source_identities =
      std::max(1, hat_source_identity_count(ctx, active_identities));
  benchmark_counter_max(benchmark_stats().hat_max_scoring_source_identities,
                        source_identities);
  if (active_identities > 0) {
    const long long inflation_ppm =
        (static_cast<long long>(source_identities) * 1000000ll +
         active_identities - 1) /
        active_identities;
    benchmark_counter_max(benchmark_stats().hat_max_scoring_source_inflation_ppm,
                          inflation_ppm);
    if (active_identities > 3 &&
        static_cast<long long>(source_identities) >
            static_cast<long long>(active_identities) * 2ll) {
      benchmark_counter_add(
          benchmark_stats().hat_scoring_source_inflation_violations);
    }
  }
  if (active_identities < static_cast<int>(context_taxa(ctx).size()) &&
      ctx.owned_active_genes == nullptr && ctx.active_genes == nullptr) {
    benchmark_counter_add(benchmark_stats().hat_root_source_scoring_calls);
  }
}

SolverContext make_hat_child_context_base(const SolverContext& parent) {
  SolverContext child;
  child.taxa_ref = &context_taxa(parent);
  child.genes_ref = &base_genes(parent);
  child.split_certificate_index = parent.split_certificate_index;
  child.consensus_ref = &context_consensus(parent);
  child.options = parent.options;
  child.workers = parent.workers;
  return child;
}

SolverContext make_hat_materialized_child_context(
    const SolverContext& parent,
    std::shared_ptr<const std::vector<GeneIndex>> active_genes,
    const std::vector<Identity>& child_identities,
    bool accounting) {
  SolverContext child = make_hat_child_context_base(parent);
  child.owned_active_genes = std::move(active_genes);
  child.active_genes_match_identities = child.owned_active_genes != nullptr;
  if (child.owned_active_genes != nullptr) {
    child.active_gene_identity_ids =
        hat_identity_ids(child_identities, accounting);
    child.active_gene_identity_taxa =
        hat_identity_taxa(child_identities, accounting);
    child.active_gene_source_identity_count =
        static_cast<int>(child_identities.size());
  }
  return child;
}

SolverContext make_hat_reused_child_context(const SolverContext& parent) {
  SolverContext child = make_hat_child_context_base(parent);
  child.owned_active_genes = parent.owned_active_genes;
  child.active_genes = parent.active_genes;
  if (child.owned_active_genes == nullptr && child.active_genes == nullptr) {
    child.active_genes = &base_genes(parent);
  }
  child.active_gene_identity_ids = parent.active_gene_identity_ids;
  child.active_gene_identity_taxa = parent.active_gene_identity_taxa;
  child.active_gene_source_identity_count =
      parent.active_gene_source_identity_count > 0
          ? parent.active_gene_source_identity_count
          : static_cast<int>(context_taxa(parent).size());
  child.active_genes_match_identities = false;
  benchmark_counter_add(benchmark_stats().hat_child_gene_bank_reuses);
  return child;
}

SolverContext make_hat_child_context(
    const SolverContext& parent,
    std::shared_ptr<const std::vector<GeneIndex>> active_genes,
    const std::vector<Identity>& child_identities,
    bool accounting) {
  if (active_genes != nullptr) {
    return make_hat_materialized_child_context(
        parent, std::move(active_genes), child_identities, accounting);
  }
  return make_hat_reused_child_context(parent);
}

long long hat_estimate_gene_bank_host_bytes(
    const std::vector<GeneIndex>& genes) {
  long long bytes = static_cast<long long>(genes.capacity()) *
                    static_cast<long long>(sizeof(GeneIndex));
  auto vector_bytes = [](const auto& values) -> long long {
    return static_cast<long long>(values.capacity()) *
           static_cast<long long>(
               sizeof(typename std::decay_t<decltype(values)>::value_type));
  };
  auto traversal_bytes = [&](const GeneTraversal& traversal) {
    return vector_bytes(traversal.parent) + vector_bytes(traversal.preorder) +
           vector_bytes(traversal.postorder) + vector_bytes(traversal.tin) +
           vector_bytes(traversal.tout);
  };
  for (const GeneIndex& gene : genes) {
    bytes += vector_bytes(gene.nodes) + vector_bytes(gene.taxon_to_node) +
             traversal_bytes(gene.root_traversal) +
             traversal_bytes(gene.reference_traversal);
    for (const GeneIndex::Node& node : gene.nodes) {
      bytes += vector_bytes(node.neighbors);
    }
  }
  return bytes;
}

long long hat_gene_bank_node_count(const std::vector<GeneIndex>& genes) {
  long long nodes = 0;
  for (const GeneIndex& gene : genes) {
    nodes += static_cast<long long>(gene.nodes.size());
  }
  return nodes;
}

std::shared_ptr<const std::vector<GeneIndex>> track_hat_child_gene_bank(
    std::unique_ptr<std::vector<GeneIndex>> bank,
    int identity_count,
    bool accounting) {
  BenchmarkStats& stats = benchmark_stats();
  const HatBankPayload payload =
      accounted_hat_bank_payload(*bank, accounting);
  if (accounting) {
    benchmark_live_payload_acquire(stats.hat_outer_live_bank_nodes,
                                   stats.hat_outer_peak_bank_nodes,
                                   payload.nodes);
    benchmark_live_payload_acquire(
        stats.hat_outer_live_bank_other_entries,
        stats.hat_outer_peak_bank_other_entries, payload.other_entries);
  }
  const long long active =
      stats.active_child_compressed_banks.fetch_add(1,
                                                    std::memory_order_relaxed) +
      1;
  benchmark_counter_max(stats.max_active_child_compressed_banks, active);
  const long long active_identity_slots =
      stats.active_child_compressed_bank_identity_slots.fetch_add(
          identity_count, std::memory_order_relaxed) +
      identity_count;
  benchmark_counter_max(stats.max_active_child_compressed_bank_identity_slots,
                        active_identity_slots);
  benchmark_counter_max(stats.max_child_compressed_bank_bytes,
                        hat_estimate_gene_bank_host_bytes(*bank));
  return std::shared_ptr<const std::vector<GeneIndex>>(
      bank.release(), [identity_count, accounting,
                       payload](const std::vector<GeneIndex>* ptr) {
        delete ptr;
        if (accounting) {
          BenchmarkStats& live_stats = benchmark_stats();
          benchmark_live_payload_release(
              live_stats.hat_outer_live_bank_other_entries,
              payload.other_entries);
          benchmark_live_payload_release(live_stats.hat_outer_live_bank_nodes,
                                         payload.nodes);
        }
        benchmark_stats().active_child_compressed_banks.fetch_sub(
            1, std::memory_order_relaxed);
        benchmark_stats().active_child_compressed_bank_identity_slots.fetch_sub(
            identity_count, std::memory_order_relaxed);
      });
}

bool should_attempt_hat_child_gene_bank(
    const std::vector<GeneIndex>& source,
    const std::vector<Identity>& child_identities) {
  return child_identities.size() > 3 && !source.empty() &&
         hat_gene_bank_node_count(source) > 0;
}

long long projected_hat_child_bank_identity_slots(
    const SolverContext& ctx,
    int parent_identity_count,
    const std::vector<Identity>& child_identities) {
  (void)parent_identity_count;
  const std::vector<GeneIndex>& source = scoring_genes(ctx);
  if (!should_attempt_hat_child_gene_bank(source, child_identities)) {
    return 0;
  }
  return static_cast<long long>(child_identities.size());
}

bool hat_source_bank_is_parallel_linear(const SolverContext& ctx,
                                        int active_identity_count) {
  constexpr long long kNodesPerGeneIdentity = 8;
  const std::vector<GeneIndex>& source = scoring_genes(ctx);
  const long long genes = static_cast<long long>(source.size());
  const long long source_identities = static_cast<long long>(
      std::max(1, hat_source_identity_count(ctx, active_identity_count)));
  if (genes == 0) return true;
  const long long nodes = hat_gene_bank_node_count(source);
  if (genes > std::numeric_limits<long long>::max() /
                  kNodesPerGeneIdentity) {
    return false;
  }
  const long long per_identity = kNodesPerGeneIdentity * genes;
  if (source_identities >
      std::numeric_limits<long long>::max() / per_identity) {
    return false;
  }
  return nodes <= per_identity * source_identities;
}

long long hat_structural_parallel_slot_budget(const SolverContext& ctx) {
  const double factor = std::max(1.0, ctx.options.hat_structural_budget_factor);
  const long double budget = std::ceil(
      static_cast<long double>(factor) *
      static_cast<long double>(
          std::max<std::size_t>(1, context_taxa(ctx).size())));
  if (!std::isfinite(budget) ||
      budget > static_cast<long double>(
                   std::numeric_limits<long long>::max())) {
    throw std::length_error(
        "HAT structural slot budget exceeds its integer domain");
  }
  return static_cast<long long>(budget);
}

class HatStructuralParallelSlotReservation {
 public:
  HatStructuralParallelSlotReservation() = default;
  HatStructuralParallelSlotReservation(
      const HatStructuralParallelSlotReservation&) = delete;
  HatStructuralParallelSlotReservation& operator=(
      const HatStructuralParallelSlotReservation&) = delete;

  ~HatStructuralParallelSlotReservation() {
    release();
  }

  bool acquire(const SolverContext& ctx, long long slots) {
    if (slots <= 0) return true;
    const long long budget = hat_structural_parallel_slot_budget(ctx);
    BenchmarkStats& stats = benchmark_stats();
    long long observed = stats.hat_active_structural_parallel_slots.load(
        std::memory_order_acquire);
    while (true) {
      const long long active_bank_slots =
          stats.active_child_compressed_bank_identity_slots.load(
              std::memory_order_acquire);
      if (observed < 0 || active_bank_slots < 0) {
        throw std::logic_error("negative HAT structural slot accounting");
      }
      if (observed > budget || slots > budget - observed) return false;
      const long long next = observed + slots;
      if (active_bank_slots > budget || next > budget - active_bank_slots) {
        return false;
      }
      if (stats.hat_active_structural_parallel_slots.compare_exchange_weak(
              observed, next, std::memory_order_acq_rel,
              std::memory_order_acquire)) {
        slots_ = slots;
        benchmark_counter_max(stats.hat_max_structural_parallel_slots, next);
        return true;
      }
    }
  }

  void release() {
    if (slots_ <= 0) return;
    benchmark_stats().hat_active_structural_parallel_slots.fetch_sub(
        slots_, std::memory_order_acq_rel);
    slots_ = 0;
  }

 private:
  long long slots_ = 0;
};

bool try_reserve_hat_parallel_children(
    const SolverContext& ctx,
    int parent_identity_count,
    const std::vector<Identity>& left,
    const std::vector<Identity>& right,
    HatStructuralParallelSlotReservation& reservation) {
  const long long projected_slots =
      projected_hat_child_bank_identity_slots(ctx, parent_identity_count,
                                               left);
  const long long projected_right = projected_hat_child_bank_identity_slots(
      ctx, parent_identity_count, right);
  if (projected_right >
      std::numeric_limits<long long>::max() - projected_slots) {
    benchmark_counter_add(
        benchmark_stats().hat_structural_parallel_rejections);
    return false;
  }
  if (!reservation.acquire(ctx, projected_slots + projected_right)) {
    benchmark_counter_add(
        benchmark_stats().hat_structural_parallel_rejections);
    return false;
  }
  benchmark_counter_add(
      benchmark_stats().hat_structural_parallel_child_subproblems);
  return true;
}

std::vector<int> source_identity_old_to_new(
    const SolverContext& ctx,
    const std::vector<int64_t>* derived_parent_identity_ids,
    const std::vector<Identity>& child_identities,
    std::vector<double>& old_scale) {
  std::vector<std::pair<int64_t, int>> child_index;
  child_index.reserve(child_identities.size());
  for (int i = 0; i < static_cast<int>(child_identities.size()); ++i) {
    child_index.emplace_back(
        child_identities[static_cast<std::size_t>(i)].id, i);
  }
  std::sort(child_index.begin(), child_index.end());
  child_index.erase(
      std::unique(child_index.begin(), child_index.end(),
                  [](const auto& left, const auto& right) {
                    return left.first == right.first;
                  }),
      child_index.end());
  const auto find_child = [&](int64_t identity_id) {
    const auto it = std::lower_bound(
        child_index.begin(), child_index.end(), identity_id,
        [](const auto& entry, int64_t sought) { return entry.first < sought; });
    return it == child_index.end() || it->first != identity_id ? -1
                                                               : it->second;
  };

  if (ctx.active_gene_identity_ids != nullptr) {
    std::vector<int> old_to_new(ctx.active_gene_identity_ids->size(), -1);
    old_scale.assign(ctx.active_gene_identity_ids->size(), 0.0);
    for (int old = 0;
         old < static_cast<int>(ctx.active_gene_identity_ids->size()); ++old) {
      const int64_t id =
          (*ctx.active_gene_identity_ids)[static_cast<std::size_t>(old)];
      const int child = find_child(id);
      if (child < 0) continue;
      old_to_new[static_cast<std::size_t>(old)] = child;
      old_scale[static_cast<std::size_t>(old)] = 1.0;
    }
    return old_to_new;
  }

  if (derived_parent_identity_ids == nullptr) {
    throw std::logic_error("identity-matched gene bank has no source identities");
  }
  std::vector<int> old_to_new(derived_parent_identity_ids->size(), -1);
  old_scale.assign(derived_parent_identity_ids->size(), 0.0);
  for (int i = 0;
       i < static_cast<int>(derived_parent_identity_ids->size()); ++i) {
    const int child = find_child(
        (*derived_parent_identity_ids)[static_cast<std::size_t>(i)]);
    if (child < 0) continue;
    old_to_new[static_cast<std::size_t>(i)] = child;
    old_scale[static_cast<std::size_t>(i)] = 1.0;
  }
  return old_to_new;
}

std::shared_ptr<const std::vector<int>> hat_scoring_identity_taxa(
    const SolverContext& ctx,
    const std::vector<Identity>& identities,
    bool accounting) {
  if (ctx.active_gene_identity_taxa != nullptr) {
    return ctx.active_gene_identity_taxa;
  }
  if (ctx.active_genes_match_identities || ctx.owned_active_genes != nullptr) {
    return hat_identity_taxa(identities, accounting);
  }
  return nullptr;
}

enum class HatChildBankDecision : std::uint8_t {
  NotNeeded,
  Materialized,
  HeavyReuse,
};

struct HatChildGeneBank {
  std::shared_ptr<const std::vector<GeneIndex>> bank;
  HatChildBankDecision decision = HatChildBankDecision::NotNeeded;
  bool attempted = false;
  long long source_nodes = 0;
  long long candidate_nodes = 0;
};

HatChildGeneBank build_hat_child_gene_bank(
    const SolverContext& ctx,
    int parent_identity_count,
    const std::vector<int64_t>* derived_parent_identity_ids,
    const std::vector<Identity>& child_identities,
    bool accounting) {
  (void)parent_identity_count;
  const std::vector<GeneIndex>& source = scoring_genes(ctx);
  HatChildGeneBank result;
  result.source_nodes = hat_gene_bank_node_count(source);
  result.attempted =
      should_attempt_hat_child_gene_bank(source, child_identities);
  if (!result.attempted) return result;

  auto bank = std::make_unique<std::vector<GeneIndex>>();
  bank->resize(source.size());
  BenchmarkStats& stats = benchmark_stats();
  benchmark_counter_add(stats.hat_child_bank_state_attempts);
  benchmark_counter_add(stats.hat_child_bank_state_input_nodes,
                        result.source_nodes);
  const auto started = std::chrono::steady_clock::now();
  // Remap tables and per-gene compression workspaces are construction-only:
  // all are destroyed before the completed bank is published and tracked.
  // parallel_for_inner bounds simultaneous per-gene workspaces by the active
  // worker budget; retained bank capacity is measured in
  // track_hat_child_gene_bank.
  if (ctx.active_gene_identity_ids != nullptr ||
      ctx.active_genes_match_identities) {
    std::vector<double> old_scale;
    std::vector<int> old_to_new = source_identity_old_to_new(
        ctx, derived_parent_identity_ids, child_identities, old_scale);
    parallel_for_inner(0, static_cast<int>(source.size()), [&](int g) {
      (*bank)[static_cast<std::size_t>(g)] = compress_gene_for_identity_remap(
          source[static_cast<std::size_t>(g)], child_identities, old_to_new,
          old_scale);
    }, 1);
  } else {
    const SubproblemWeights weights =
        build_subproblem_weights(child_identities, context_taxa(ctx).size());
    parallel_for_inner(0, static_cast<int>(source.size()), [&](int g) {
      (*bank)[static_cast<std::size_t>(g)] = compress_gene_for_subproblem(
          source[static_cast<std::size_t>(g)], child_identities, weights);
    }, 1);
  }
  result.candidate_nodes = hat_gene_bank_node_count(*bank);
  benchmark_counter_add(stats.hat_child_bank_state_output_nodes,
                        result.candidate_nodes);
  benchmark_timer_add(stats.compression_materialization_ns,
                      std::chrono::steady_clock::now() - started);
  const HatBankPayload candidate_payload =
      accounted_hat_bank_payload(*bank, accounting);
  HatBankPayloadLease candidate_live(accounting, candidate_payload);
  const long long candidate_bytes = hat_estimate_gene_bank_host_bytes(*bank);
  benchmark_counter_max(stats.max_child_compressed_bank_bytes,
                        candidate_bytes);

  // A newly retained bank must shrink physical topology, not merely owner
  // count. The additive per-gene vector headers are charged separately as
  // O(k) per recursion frame; all topology/traversal payload is linear in this
  // node state. The pair lemma checked below guarantees that at most one child
  // can miss this factor-3/4 retention gate.
  const long long state_limit =
      result.source_nodes - (result.source_nodes + 3ll) / 4ll;
  if (result.candidate_nodes > state_limit) {
    benchmark_counter_add(stats.hat_child_bank_state_heavy_reuses);
    result.decision = HatChildBankDecision::HeavyReuse;
    return result;
  }

  const long long accepted_ratio_ppm = static_cast<long long>(std::llround(
      static_cast<long double>(result.candidate_nodes) * 1000000.0L /
      static_cast<long double>(result.source_nodes)));
  benchmark_counter_max(stats.hat_child_bank_state_max_accepted_ratio_ppm,
                        accepted_ratio_ppm);
  benchmark_counter_add(stats.hat_child_bank_state_accepts);
  benchmark_counter_add(stats.compression_materializations);
  result.decision = HatChildBankDecision::Materialized;
  candidate_live.release();
  result.bank = track_hat_child_gene_bank(
      std::move(bank), static_cast<int>(child_identities.size()), accounting);
  return result;
}

void validate_hat_child_bank_pair(const HatChildGeneBank& left,
                                  const HatChildGeneBank& right) {
  if (!left.attempted || !right.attempted || left.source_nodes == 0 ||
      right.source_nodes == 0 ||
      left.decision != HatChildBankDecision::HeavyReuse ||
      right.decision != HatChildBankDecision::HeavyReuse) {
    return;
  }
  // In each reduced source tree, let A/B be the disjoint selected physical
  // leaves, U the inactive leaves, and I the internal vertices. Child
  // reduction can retain at most A+B selected leaves plus two copies of an
  // internal vertex, so |L|+|R| <= A+B+2I. Let L=A+B+U. A parsed source can
  // retain at most one unary root and all other unlabeled internal vertices
  // have degree at least three; a compressed source removes every unlabeled
  // degree-0/1/2 vertex. Thus I<=L in every source, which implies
  // 2(|L|+|R|) <= 2L+4I <= 3L+3I = 3*source. Summing over genes proves that
  // both integer children cannot exceed floor(3*source/4). Reaching this
  // branch means either the source invariant or the compressor violated the
  // space proof.
  BenchmarkStats& stats = benchmark_stats();
  benchmark_counter_add(stats.hat_child_bank_state_pair_violations);
  stats.complexity_contract_satisfied.store(0, std::memory_order_relaxed);
  throw std::logic_error(
      "HAT child-bank state gate produced two heavy children");
}

void split_identities(std::vector<Identity> identities,
                      const Candidate& candidate,
                      std::vector<Identity>& left,
                      std::vector<Identity>& right,
                      HatIdentityPayload& left_payload,
                      HatIdentityPayload& right_payload) {
  left.clear();
  right.clear();
  const int left_count = std::clamp(candidate.split.count(), 0,
                                    static_cast<int>(identities.size()));
  left.reserve(static_cast<std::size_t>(left_count));
  right.reserve(identities.size() - static_cast<std::size_t>(left_count));
  for (int i = 0; i < static_cast<int>(identities.size()); ++i) {
    Identity& identity = identities[static_cast<std::size_t>(i)];
    const bool goes_left = candidate.split.test(i);
    (goes_left ? left_payload : right_payload) +=
        hat_identity_payload(identity);
    (goes_left ? left : right).push_back(std::move(identity));
  }
}

void stable_singleton_split(std::vector<Identity> identities,
                            std::vector<Identity>& left,
                            std::vector<Identity>& right,
                            HatIdentityPayload& left_payload,
                            HatIdentityPayload& right_payload) {
  std::sort(identities.begin(), identities.end(),
            [](const Identity& a, const Identity& b) {
              if (a.stable_key != b.stable_key) return a.stable_key < b.stable_key;
              return a.id < b.id;
            });
  left.clear();
  right.clear();
  left.reserve(1);
  right.reserve(identities.size() - 1);
  Identity& singleton = identities.front();
  left_payload += hat_identity_payload(singleton);
  left.push_back(std::move(singleton));
  for (int i = 1; i < static_cast<int>(identities.size()); ++i) {
    Identity& identity = identities[static_cast<std::size_t>(i)];
    right_payload += hat_identity_payload(identity);
    right.push_back(std::move(identity));
  }
}

Candidate forced_hat_candidate(int identities, HatRecursionTestShape shape) {
  Candidate candidate;
  candidate.split = Bitset(identities);
  const int left = shape == HatRecursionTestShape::Balanced
                       ? identities / 2
                       : 1;
  for (int i = 0; i < left; ++i) candidate.split.set(i);
  return candidate;
}

bool legal_hat_split(const SolverContext& ctx,
                     int left_count,
                     int right_count) {
  (void)ctx;
  const int total = left_count + right_count;
  const int min_side =
      std::min(kHatMinSide, std::max(1, total / 2));
  return left_count >= min_side && right_count >= min_side;
}

AttachmentEdge first_hat_attachment_edge(const SolverTree& tree) {
  AttachmentEdge edge;
  for (int a = 0; a < static_cast<int>(tree.nodes.size()); ++a) {
    if (!tree.nodes[static_cast<std::size_t>(a)].active) continue;
    for (int b : tree.nodes[static_cast<std::size_t>(a)].neighbors) {
      if (a < b && b >= 0 && b < static_cast<int>(tree.nodes.size()) &&
          tree.nodes[static_cast<std::size_t>(b)].active) {
        edge.a = a;
        edge.b = b;
        return edge;
      }
    }
  }
  for (int a = 0; a < static_cast<int>(tree.nodes.size()); ++a) {
    if (tree.nodes[static_cast<std::size_t>(a)].active) {
      edge.a = a;
      edge.b = -1;
      return edge;
    }
  }
  return edge;
}

SolverTree solve_hat_recursive(SolverContext& ctx,
                               std::vector<Identity> identities,
                               uint64_t path,
                               int depth,
                               tf::Subflow* subflow,
                               const HatIdentityTaxonCatalog& identity_catalog,
                               const HatRecursionControl& control,
                               HatIdentityPayloadLease identity_payload) {
  DebugScope scope("hat.solve_recursive");
  const HatIdentityPayload active_payload = identity_payload.payload();
  const int active_identities = static_cast<int>(identities.size());
  const int active_taxa = static_cast<int>(active_payload.taxa);
  if (control.accounting) {
    BenchmarkStats& stats = benchmark_stats();
    benchmark_counter_add(stats.hat_outer_recursive_calls);
    benchmark_counter_add(stats.hat_outer_active_identity_mass,
                          active_payload.identities);
    benchmark_counter_add(stats.hat_outer_active_support_word_mass,
                          active_payload.support_words);
    benchmark_counter_add(stats.hat_outer_active_taxon_mass,
                          active_payload.taxa);
    benchmark_counter_max(stats.hat_outer_max_depth, depth);
  }
  if (identities.size() <= 3) {
    SolverTree result = hat_base_tree(identities);
    HatCompletedTreeLease result_live(
        control.accounting,
        accounted_hat_tree_payload(result, control.accounting));
    identity_payload.release();
    return result;
  }

  if (control.accounting) {
    benchmark_counter_add(benchmark_stats().hat_outer_nonbase_calls);
  }
  if (!ctx.options.benchmark_json.empty()) {
    benchmark_subproblem_start(path, depth, active_identities, active_taxa,
                               identities);
  }
  benchmark_counter_add(benchmark_stats().subproblems);
  benchmark_counter_add(benchmark_stats().hat_subproblems);
  benchmark_counter_max(benchmark_stats().max_depth, depth);
  benchmark_counter_max(benchmark_stats().max_active_identities,
                        active_identities);
  benchmark_counter_max(benchmark_stats().max_active_taxa, active_taxa);
  benchmark_counter_max(benchmark_stats().hat_max_active_taxa, active_taxa);
  note_hat_scoring_source(ctx, active_identities);
  progress_note_subproblem(path, depth, active_identities, active_taxa);

  progress_note_phase(path, "hat-split-score");
  std::vector<Identity> left;
  std::vector<Identity> right;
  std::shared_ptr<const std::vector<int>> parent_identity_taxa;
  std::vector<int64_t> derived_parent_identity_ids;
  HatPayloadLease derived_identity_payload;
  HatIdentityPayload left_payload;
  HatIdentityPayload right_payload;
  HatIdentityPayloadLease left_identity_payload;
  HatIdentityPayloadLease right_identity_payload;
  const bool needs_derived_parent_ids =
      ctx.active_gene_identity_ids == nullptr &&
      ctx.active_genes_match_identities;
  {
    HatSplitChoice split =
        hat_choose_active_split_hdt(ctx, identities, path, depth);
    if (control.force_shape) {
      split.candidate = forced_hat_candidate(active_identities, control.shape);
    }
    progress_note_phase(path, "hat-split-choose");
    parent_identity_taxa =
        hat_scoring_identity_taxa(ctx, identities, control.accounting);
    if (needs_derived_parent_ids) {
      derived_parent_identity_ids.reserve(identities.size());
      for (const Identity& identity : identities) {
        derived_parent_identity_ids.push_back(identity.id);
      }
      BenchmarkStats& stats = benchmark_stats();
      derived_identity_payload.acquire(
          control.accounting, stats.hat_outer_live_bank_identity_entries,
          stats.hat_outer_peak_bank_identity_entries,
          static_cast<long long>(derived_parent_identity_ids.size()));
    }

    int candidate_left_count = 0;
    for (int i = 0; i < active_identities; ++i) {
      candidate_left_count += split.candidate.split.test(i) ? 1 : 0;
    }
    if (control.accounting) {
      benchmark_counter_add(benchmark_stats().hat_outer_split_identity_tests,
                            active_identities);
    }
    const int candidate_right_count =
        active_identities - candidate_left_count;
    if (legal_hat_split(ctx, candidate_left_count, candidate_right_count)) {
      split_identities(std::move(identities), split.candidate, left, right,
                       left_payload, right_payload);
    } else if (ctx.options.hat_one_third_singleton) {
      throw std::logic_error(
          "HAT one-third/singleton mode selected an illegal split");
    } else {
      stable_singleton_split(std::move(identities), left, right, left_payload,
                             right_payload);
    }
    identity_payload.release();
    left_identity_payload =
        HatIdentityPayloadLease(control.accounting, left_payload);
    right_identity_payload =
        HatIdentityPayloadLease(control.accounting, right_payload);
    if (control.accounting) {
      BenchmarkStats& stats = benchmark_stats();
      benchmark_counter_add(stats.hat_outer_partition_identity_moves,
                            active_payload.identities);
      benchmark_counter_add(stats.hat_outer_partition_support_word_moves,
                            active_payload.support_words);
      benchmark_counter_add(stats.hat_outer_partition_taxon_moves,
                            active_payload.taxa);
      const int left_size = static_cast<int>(left.size());
      const int right_size = static_cast<int>(right.size());
      if (left_size == right_size) {
        benchmark_counter_add(stats.hat_outer_balanced_splits);
      }
      if (std::min(left_size, right_size) == 1 &&
          std::max(left_size, right_size) == active_identities - 1) {
        benchmark_counter_add(stats.hat_outer_maximally_skewed_splits);
      }
    }
    if (!ctx.options.benchmark_json.empty()) {
      benchmark_subproblem_set_final(path, left, right);
    }
  }

  if (ctx.options.debug_stop_after_root && depth == 0) {
    progress_note_split(path, depth, active_identities, active_taxa,
                        static_cast<int>(left.size()),
                        static_cast<int>(right.size()));
    SolverTree left_tree = hat_base_tree(left);
    SolverTree right_tree = hat_base_tree(right);
    HatCompletedTreeLease left_live(
        control.accounting,
        accounted_hat_tree_payload(left_tree, control.accounting));
    HatCompletedTreeLease right_live(
        control.accounting,
        accounted_hat_tree_payload(right_tree, control.accounting));
    left_identity_payload.release();
    right_identity_payload.release();
    const AttachmentEdge left_edge = first_hat_attachment_edge(left_tree);
    const AttachmentEdge right_edge = first_hat_attachment_edge(right_tree);
    HatTreePayload glue_payload = accounted_hat_glue_scratch_payload(
        right_tree, left_edge, right_edge, control.accounting);
    HatCompletedTreeLease glue_live(control.accounting, glue_payload);
    SolverTree result = hat_glue(std::move(left_tree), std::move(right_tree),
                                 left_edge, right_edge);
    left_live.release();
    right_live.release();
    glue_live.release();
    HatCompletedTreeLease result_live(
        control.accounting,
        accounted_hat_tree_payload(result, control.accounting));
    return result;
  }

  progress_note_split(path, depth, active_identities, active_taxa,
                      static_cast<int>(left.size()),
                      static_cast<int>(right.size()));

  const uint64_t left_path = hat_child_path(path, 1, depth);
  const uint64_t right_path = hat_child_path(path, 2, depth);
  SolverTree left_tree;
  SolverTree right_tree;
  HatCompletedTreeLease left_tree_payload;
  HatCompletedTreeLease right_tree_payload;

  HatStructuralParallelSlotReservation parallel_slot_reservation;
  const int scoring_parallelism = std::min<int>(
      scoring_genes(ctx).size(), get_num_workers());
  const bool base_parallel_candidate =
      subflow != nullptr &&
      active_identities >= kHatParallelSubproblemMinIdentities &&
      scoring_parallelism < get_num_workers();
  const bool physical_parallel_ok =
      !base_parallel_candidate ||
      hat_source_bank_is_parallel_linear(ctx, active_identities);
  if (base_parallel_candidate && !physical_parallel_ok) {
    benchmark_counter_add(
        benchmark_stats().hat_structural_parallel_physical_rejections);
  }
  const bool parallel_candidate =
      base_parallel_candidate && physical_parallel_ok;
  const bool parallel_children =
      parallel_candidate &&
      try_reserve_hat_parallel_children(ctx, active_identities, left, right,
                                        parallel_slot_reservation);
  if (parallel_children) {
    progress_note_phase(path, "hat-child-bank");
    const std::vector<int64_t>* derived_ids =
        needs_derived_parent_ids ? &derived_parent_identity_ids : nullptr;
    HatChildGeneBank left_bank =
        build_hat_child_gene_bank(ctx, active_identities, derived_ids, left,
                                  control.accounting);
    HatChildGeneBank right_bank =
        build_hat_child_gene_bank(ctx, active_identities, derived_ids, right,
                                  control.accounting);
    validate_hat_child_bank_pair(left_bank, right_bank);
    std::vector<int64_t>().swap(derived_parent_identity_ids);
    derived_identity_payload.release();
    parallel_slot_reservation.release();
    SolverContext left_ctx =
        make_hat_child_context(ctx, std::move(left_bank.bank), left,
                               control.accounting);
    SolverContext right_ctx =
        make_hat_child_context(ctx, std::move(right_bank.bank), right,
                               control.accounting);
    progress_note_phase(path, "hat-children");
    struct ParallelFrame {
      SolverTree first;
      SolverTree second;
      HatIdentityPayloadLease first_identities;
      HatIdentityPayloadLease second_identities;
      HatCompletedTreeLease first_tree;
      HatCompletedTreeLease second_tree;
      std::exception_ptr first_failure;
      std::exception_ptr second_failure;
    };
    auto frame = std::make_shared<ParallelFrame>();
    frame->first_identities = std::move(left_identity_payload);
    frame->second_identities = std::move(right_identity_payload);
    tf::Task lt = subflow->emplace(
        [frame, ids = std::move(left), child_ctx = std::move(left_ctx),
         left_path, depth, &identity_catalog,
         &control](tf::Subflow& sf) mutable {
          try {
            frame->first = solve_hat_recursive(
                child_ctx, std::move(ids), left_path, depth + 1, &sf,
                identity_catalog, control,
                std::move(frame->first_identities));
            frame->first_tree.acquire(control.accounting,
                                      accounted_hat_tree_payload(
                                          frame->first, control.accounting));
          } catch (...) {
            frame->first_failure = std::current_exception();
          }
        });
    tf::Task rt = subflow->emplace(
        [frame, ids = std::move(right), child_ctx = std::move(right_ctx),
         right_path, depth, &identity_catalog,
         &control](tf::Subflow& sf) mutable {
          try {
            frame->second = solve_hat_recursive(
                child_ctx, std::move(ids), right_path, depth + 1, &sf,
                identity_catalog, control,
                std::move(frame->second_identities));
            frame->second_tree.acquire(control.accounting,
                                       accounted_hat_tree_payload(
                                           frame->second,
                                           control.accounting));
          } catch (...) {
            frame->second_failure = std::current_exception();
          }
        });
    auto out = std::make_shared<SolverTree>();
    tf::Task join = subflow->emplace([&ctx, &identity_catalog, frame, out,
                                      parent_identity_taxa,
                                      active_identities, &control]() {
      rethrow_hat_parallel_child_failure(frame->first_failure,
                                         frame->second_failure);
      HatTrackedTaxa left_taxa = hat_taxa_from_solver_tree(
          frame->first, identity_catalog, control.accounting);
      HatTrackedTaxa right_taxa = hat_taxa_from_solver_tree(
          frame->second, identity_catalog, control.accounting);
      progress_note_stage("hat-one-x");
      progress_note_one_x_start(active_identities,
                                static_cast<int>(left_taxa.values.size()),
                                static_cast<int>(right_taxa.values.size()));
      const auto edges = hat_best_one_x_edges_hdt(
          ctx, frame->first, frame->second, left_taxa.values,
          right_taxa.values,
          parent_identity_taxa.get());
      progress_note_one_x_done(active_identities);
      HatTreePayload glue_payload = accounted_hat_glue_scratch_payload(
          frame->second, edges.first, edges.second, control.accounting);
      HatCompletedTreeLease glue_live(control.accounting, glue_payload);
      *out = hat_glue(std::move(frame->first), std::move(frame->second),
                      edges.first, edges.second);
      frame->first_tree.release();
      frame->second_tree.release();
      glue_live.release();
    });
    lt.precede(join);
    rt.precede(join);
    subflow->join();
    benchmark_counter_add(benchmark_stats().completed_subproblems);
    HatCompletedTreeLease result_live(
        control.accounting,
        accounted_hat_tree_payload(*out, control.accounting));
    return std::move(*out);
  }

  benchmark_counter_add(
      benchmark_stats().hat_structural_serial_child_subproblems);

  const std::vector<int64_t>* derived_ids =
      needs_derived_parent_ids ? &derived_parent_identity_ids : nullptr;
  progress_note_phase(path, "hat-left-child-bank");
  HatChildGeneBank prepared_left_bank = build_hat_child_gene_bank(
      ctx, active_identities, derived_ids, left, control.accounting);
  std::optional<HatChildGeneBank> prepared_right_bank;
  if (needs_derived_parent_ids ||
      prepared_left_bank.decision == HatChildBankDecision::HeavyReuse) {
    progress_note_phase(path, "hat-right-child-bank");
    prepared_right_bank.emplace(build_hat_child_gene_bank(
        ctx, active_identities, derived_ids, right, control.accounting));
    validate_hat_child_bank_pair(prepared_left_bank, *prepared_right_bank);
  }
  if (needs_derived_parent_ids) {
    std::vector<int64_t>().swap(derived_parent_identity_ids);
    derived_identity_payload.release();
  }

  right_identity_payload.mark_continuation();
  progress_note_phase(path, "hat-left-child-bank");
  {
    SolverContext left_ctx =
        make_hat_child_context(ctx, std::move(prepared_left_bank.bank), left,
                               control.accounting);
    progress_note_stage("hat-left-child");
    left_tree = solve_hat_recursive(left_ctx, std::move(left), left_path,
                                    depth + 1, nullptr, identity_catalog,
                                    control,
                                    std::move(left_identity_payload));
    left_tree_payload.acquire(
        control.accounting,
        accounted_hat_tree_payload(left_tree, control.accounting));
  }
  right_identity_payload.unmark_continuation();
  trim_hat_allocator_after_large_release(active_identities);

  progress_note_phase(path, "hat-right-child-bank");
  {
    if (!prepared_right_bank.has_value()) {
      prepared_right_bank.emplace(build_hat_child_gene_bank(
          ctx, active_identities, nullptr, right, control.accounting));
      validate_hat_child_bank_pair(prepared_left_bank,
                                   *prepared_right_bank);
    }
    SolverContext right_ctx =
        make_hat_child_context(ctx, std::move(prepared_right_bank->bank), right,
                               control.accounting);
    progress_note_stage("hat-right-child");
    right_tree = solve_hat_recursive(right_ctx, std::move(right), right_path,
                                     depth + 1, nullptr, identity_catalog,
                                     control,
                                     std::move(right_identity_payload));
    right_tree_payload.acquire(control.accounting,
                               accounted_hat_tree_payload(
                                   right_tree, control.accounting));
  }
  trim_hat_allocator_after_large_release(active_identities);
  HatTrackedTaxa left_taxa = hat_taxa_from_solver_tree(
      left_tree, identity_catalog, control.accounting);
  HatTrackedTaxa right_taxa = hat_taxa_from_solver_tree(
      right_tree, identity_catalog, control.accounting);
  progress_note_stage("hat-one-x");
  progress_note_one_x_start(active_identities,
                            static_cast<int>(left_taxa.values.size()),
                            static_cast<int>(right_taxa.values.size()));
  const auto edges = hat_best_one_x_edges_hdt(
      ctx, left_tree, right_tree, left_taxa.values, right_taxa.values,
      parent_identity_taxa.get());
  progress_note_one_x_done(active_identities);
  benchmark_counter_add(benchmark_stats().completed_subproblems);
  HatTreePayload glue_payload = accounted_hat_glue_scratch_payload(
      right_tree, edges.first, edges.second, control.accounting);
  HatCompletedTreeLease glue_live(control.accounting, glue_payload);
  SolverTree result = hat_glue(std::move(left_tree), std::move(right_tree),
                               edges.first, edges.second);
  left_tree_payload.release();
  right_tree_payload.release();
  glue_live.release();
  HatCompletedTreeLease result_live(
      control.accounting,
      accounted_hat_tree_payload(result, control.accounting));
  return result;
}

SolverTree solve_hat_v1_controlled(SolverContext& ctx,
                                   std::vector<Identity> identities,
                                   uint64_t path,
                                   int depth,
                                   tf::Subflow* subflow,
                                   const HatRecursionControl& control) {
  const HatIdentityPayload root_identity_payload =
      hat_identity_payload(identities);
  HatIdentityPayloadLease identities_live(control.accounting,
                                          root_identity_payload);
  const HatBankPayload root_bank_payload =
      accounted_hat_bank_payload(scoring_genes(ctx), control.accounting);
  HatBankPayloadLease root_bank_live(control.accounting, root_bank_payload);
  if (control.accounting) {
    benchmark_counter_add(benchmark_stats().hat_outer_root_bank_nodes,
                          root_bank_payload.nodes);
    benchmark_counter_add(
        benchmark_stats().hat_outer_root_bank_other_entries,
        root_bank_payload.other_entries);
  }
  const HatIdentityTaxonCatalog identity_catalog =
      make_hat_identity_taxon_catalog(identities, control.accounting);
  BenchmarkStats& stats = benchmark_stats();
  HatPayloadLease catalog_live;
  catalog_live.acquire(
      control.accounting, stats.hat_outer_live_catalog_entries,
      stats.hat_outer_peak_catalog_entries,
      static_cast<long long>(identity_catalog.entries.size() +
                             identity_catalog.taxa.size() +
                             identity_catalog.entry_by_identity.size()));
  SolverTree tree = solve_hat_recursive(
      ctx, std::move(identities), path, depth, subflow, identity_catalog,
      control, std::move(identities_live));
  HatCompletedTreeLease final_tree_live(control.accounting,
                                        accounted_hat_tree_payload(
                                            tree, control.accounting));
  if (depth == 0) {
    run_hat_nni_repair(ctx, tree, identity_catalog.taxa);
  }
  return tree;
}

}  // namespace

SolverTree solve_hat_v1(SolverContext& ctx,
                        std::vector<Identity> identities,
                        uint64_t path,
                        int depth,
                        tf::Subflow* subflow) {
  const HatRecursionControl control{
      !ctx.options.benchmark_json.empty(), false,
      HatRecursionTestShape::Balanced};
  return solve_hat_v1_controlled(ctx, std::move(identities), path, depth,
                                 subflow, control);
}

SolverTree solve_hat_v1_for_test(SolverContext& ctx,
                                 std::vector<Identity> identities,
                                 HatRecursionTestShape shape) {
  const HatRecursionControl control{true, true, shape};
  if (ctx.workers > 1 && get_num_workers() > 1) {
    SolverTree tree;
    tf::Taskflow solve_flow;
    solve_flow.emplace([&](tf::Subflow& sf) {
      tree = solve_hat_v1_controlled(ctx, std::move(identities),
                                     0x6a09e667f3bcc909ULL, 0, &sf,
                                     control);
    });
    run_inner(solve_flow);
    return tree;
  }
  return solve_hat_v1_controlled(ctx, std::move(identities),
                                 0x6a09e667f3bcc909ULL, 0, nullptr, control);
}

SolverTree solve_hat_consensus(SolverContext& ctx, bool run_nni) {
  const ConsensusTree& consensus = context_consensus(ctx);
  SolverTree tree;
  if (consensus.root < 0 || consensus.nodes.empty()) return tree;

  tree.nodes.resize(consensus.nodes.size());
  for (int i = 0; i < static_cast<int>(consensus.nodes.size()); ++i) {
    const ConsensusTree::Node& source =
        consensus.nodes[static_cast<std::size_t>(i)];
    SolverTree::Node& node = tree.nodes[static_cast<std::size_t>(i)];
    node.identity_id = source.taxon;
    node.active = true;
    node.neighbors.reserve(source.children.size() + (source.parent >= 0 ? 1 : 0));
    if (source.parent >= 0) node.neighbors.push_back(source.parent);
    for (int child : source.children) node.neighbors.push_back(child);
  }

  if (run_nni) {
    run_hat_nni_repair(
        ctx, tree,
        hat_taxa_from_identities(make_initial_identities(context_taxa(ctx))));
  }
  return tree;
}

}  // namespace wqfm::internal
