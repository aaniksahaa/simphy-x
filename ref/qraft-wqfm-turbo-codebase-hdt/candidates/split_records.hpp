#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace wqfm {
struct Options;
}

namespace wqfm::internal {

struct SplitKey256 {
  std::uint32_t size = 0;
  std::array<std::uint64_t, 4> h{};

  bool operator==(const SplitKey256&) const = default;
};

bool split_key_less(const SplitKey256& left, const SplitKey256& right);

struct SplitTaxonFingerprints256 {
  std::vector<std::array<std::uint64_t, 4>> word;
  std::array<std::uint64_t, 4> all{};
};

SplitTaxonFingerprints256 build_split_taxon_fingerprints256(
    int taxon_count,
    const Options& options);

// A physical eligible edge record. Sorting by (key, gene, edge) makes the
// first record of every (key, gene) run the deterministic representative.
struct SplitOccurrence256 {
  SplitKey256 key;
  std::uint32_t gene = 0;
  std::uint32_t edge = 0;
  std::uint32_t node = 0;
  std::uint32_t one_taxon = 0;
};

// One record per distinct key after per-gene deduplication. split_id is the
// zero-based rank of key in SplitKey256 lexicographic order.
struct SplitAggregate256 {
  SplitKey256 key;
  std::uint32_t split_id = 0;
  std::uint32_t count = 0;
  std::uint32_t first_tree = 0;
  std::uint32_t first_edge = 0;
  std::uint32_t rep_tree = 0;
  std::uint32_t rep_node = 0;
  std::uint32_t one_taxon = 0;
};

struct SplitRadixStats {
  long long occurrence_passes = 0;
  long long occurrence_items = 0;
  long long priority_passes = 0;
  long long priority_items = 0;
  long long emitted_occurrences = 0;
  long long deduplicated_gene_splits = 0;
  long long distinct_splits = 0;

  long long total_passes() const {
    return occurrence_passes + priority_passes;
  }
  long long total_items() const {
    return occurrence_items + priority_items;
  }
};

inline constexpr int kSplitOccurrenceRadixPasses = 22;
inline constexpr int kSplitPriorityRadixPasses = 8;

bool split_occurrence_less(const SplitOccurrence256& left,
                           const SplitOccurrence256& right);
bool split_priority_less_by_key(const SplitAggregate256& left,
                                const SplitAggregate256& right);

// Sorts occurrences by (key, gene, edge), removes repeated observations of a
// key in one gene, aggregates support, and assigns dense key-rank split ids.
// The returned vector remains in key/split-id order.
std::vector<SplitAggregate256> aggregate_split_occurrences256(
    std::vector<SplitOccurrence256>& occurrences,
    SplitRadixStats& stats);

// Produces the MRE priority order
// (count descending, first gene, first edge, key) using split_id for the final
// key tie-break. split_id is order-isomorphic to the complete key.
void radix_sort_split_priority256(
    std::vector<SplitAggregate256>& aggregates,
    SplitRadixStats& stats);

}  // namespace wqfm::internal
