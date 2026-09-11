#include "candidates/split_records.hpp"

#include "internal.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <tuple>

namespace wqfm::internal {
namespace {

constexpr int kRadixDigitBits = 16;
constexpr std::size_t kRadixBuckets = std::size_t{1} << kRadixDigitBits;

template <class Value, class Digit>
void stable_radix_pass(std::vector<Value>& values,
                       std::vector<Value>& scratch,
                       std::vector<std::size_t>& offsets,
                       Digit&& digit,
                       long long& pass_count,
                       long long& item_count) {
  std::fill(offsets.begin(), offsets.end(), std::size_t{0});
  for (const Value& value : values) {
    ++offsets[static_cast<std::size_t>(digit(value))];
  }
  std::size_t prefix = 0;
  for (std::size_t& offset : offsets) {
    const std::size_t count = offset;
    offset = prefix;
    prefix += count;
  }
  for (const Value& value : values) {
    scratch[offsets[static_cast<std::size_t>(digit(value))]++] = value;
  }
  values.swap(scratch);
  ++pass_count;
  item_count += static_cast<long long>(values.size());
}

template <class Value, class Getter>
void radix_u32(std::vector<Value>& values,
               std::vector<Value>& scratch,
               std::vector<std::size_t>& offsets,
               Getter&& getter,
               long long& pass_count,
               long long& item_count) {
  for (int shift = 0; shift < 32; shift += kRadixDigitBits) {
    stable_radix_pass(
        values, scratch, offsets,
        [&](const Value& value) {
          return static_cast<std::uint16_t>(getter(value) >> shift);
        },
        pass_count, item_count);
  }
}

template <class Value, class Getter>
void radix_u64(std::vector<Value>& values,
               std::vector<Value>& scratch,
               std::vector<std::size_t>& offsets,
               Getter&& getter,
               long long& pass_count,
               long long& item_count) {
  for (int shift = 0; shift < 64; shift += kRadixDigitBits) {
    stable_radix_pass(
        values, scratch, offsets,
        [&](const Value& value) {
          return static_cast<std::uint16_t>(getter(value) >> shift);
        },
        pass_count, item_count);
  }
}

void radix_sort_occurrences(std::vector<SplitOccurrence256>& occurrences,
                            SplitRadixStats& stats) {
  if (occurrences.empty()) return;
  std::vector<SplitOccurrence256> scratch(occurrences.size());
  std::vector<std::size_t> offsets(kRadixBuckets);

  // Stable LSD schedule for the lexicographic tuple
  // (size, h[0], h[1], h[2], h[3], gene, edge).
  radix_u32(occurrences, scratch, offsets,
            [](const SplitOccurrence256& value) { return value.edge; },
            stats.occurrence_passes, stats.occurrence_items);
  radix_u32(occurrences, scratch, offsets,
            [](const SplitOccurrence256& value) { return value.gene; },
            stats.occurrence_passes, stats.occurrence_items);
  for (int word = 3; word >= 0; --word) {
    radix_u64(
        occurrences, scratch, offsets,
        [word](const SplitOccurrence256& value) {
          return value.key.h[static_cast<std::size_t>(word)];
        },
        stats.occurrence_passes, stats.occurrence_items);
  }
  radix_u32(occurrences, scratch, offsets,
            [](const SplitOccurrence256& value) { return value.key.size; },
            stats.occurrence_passes, stats.occurrence_items);

  if (stats.occurrence_passes != kSplitOccurrenceRadixPasses) {
    throw std::logic_error("split-occurrence radix schedule mismatch");
  }
}

}  // namespace

bool split_key_less(const SplitKey256& left, const SplitKey256& right) {
  if (left.size != right.size) return left.size < right.size;
  return left.h < right.h;
}

SplitTaxonFingerprints256 build_split_taxon_fingerprints256(
    int taxon_count,
    const Options& options) {
  if (taxon_count < 0) {
    throw std::length_error("negative split-fingerprint taxon count");
  }
  SplitTaxonFingerprints256 out;
  out.word.resize(static_cast<std::size_t>(taxon_count));
  for (int taxon = 0; taxon < taxon_count; ++taxon) {
    for (int word = 0; word < 4; ++word) {
      const std::uint64_t value = stable_hash_u64(
          static_cast<std::uint64_t>(taxon) + 1,
          static_cast<std::uint64_t>(word) + 0x517cc1b727220a95ULL,
          options.seed ^ 0x9e3779b97f4a7c15ULL);
      out.word[static_cast<std::size_t>(taxon)]
              [static_cast<std::size_t>(word)] = value;
      out.all[static_cast<std::size_t>(word)] += value;
    }
  }
  return out;
}

bool split_occurrence_less(const SplitOccurrence256& left,
                           const SplitOccurrence256& right) {
  if (split_key_less(left.key, right.key)) return true;
  if (split_key_less(right.key, left.key)) return false;
  return std::tie(left.gene, left.edge) < std::tie(right.gene, right.edge);
}

bool split_priority_less_by_key(const SplitAggregate256& left,
                                const SplitAggregate256& right) {
  if (left.count != right.count) return left.count > right.count;
  if (left.first_tree != right.first_tree) {
    return left.first_tree < right.first_tree;
  }
  if (left.first_edge != right.first_edge) {
    return left.first_edge < right.first_edge;
  }
  return split_key_less(left.key, right.key);
}

std::vector<SplitAggregate256> aggregate_split_occurrences256(
    std::vector<SplitOccurrence256>& occurrences,
    SplitRadixStats& stats) {
  if (occurrences.size() >
      static_cast<std::size_t>(std::numeric_limits<long long>::max())) {
    throw std::length_error("split occurrence count exceeds telemetry domain");
  }
  stats.emitted_occurrences = static_cast<long long>(occurrences.size());
  radix_sort_occurrences(occurrences, stats);

  std::vector<SplitAggregate256> aggregates;
  aggregates.reserve(occurrences.size());
  std::size_t begin = 0;
  long long gene_splits = 0;
  while (begin < occurrences.size()) {
    const SplitOccurrence256& first = occurrences[begin];
    if (aggregates.size() >
        static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
      throw std::length_error("dense split id exceeds its 32-bit domain");
    }
    SplitAggregate256 aggregate;
    aggregate.key = first.key;
    aggregate.split_id = static_cast<std::uint32_t>(aggregates.size());
    aggregate.first_tree = first.gene;
    aggregate.first_edge = first.edge;
    aggregate.rep_tree = first.gene;
    aggregate.rep_node = first.node;
    aggregate.one_taxon = first.one_taxon;

    std::uint32_t previous_gene = std::numeric_limits<std::uint32_t>::max();
    std::size_t end = begin;
    while (end < occurrences.size() && occurrences[end].key == first.key) {
      const std::uint32_t gene = occurrences[end].gene;
      if (gene != previous_gene) {
        if (aggregate.count == std::numeric_limits<std::uint32_t>::max()) {
          throw std::length_error("split support exceeds its 32-bit domain");
        }
        ++aggregate.count;
        ++gene_splits;
        previous_gene = gene;
      }
      ++end;
    }
    aggregates.push_back(aggregate);
    begin = end;
  }
  stats.deduplicated_gene_splits = gene_splits;
  stats.distinct_splits = static_cast<long long>(aggregates.size());
  return aggregates;
}

void radix_sort_split_priority256(
    std::vector<SplitAggregate256>& aggregates,
    SplitRadixStats& stats) {
  if (aggregates.empty()) return;
  std::vector<SplitAggregate256> scratch(aggregates.size());
  std::vector<std::size_t> offsets(kRadixBuckets);

  // Stable LSD schedule for
  // (~count, first_tree, first_edge, split_id). Because split_id is the key
  // rank, this is exactly the legacy key-based priority comparator.
  radix_u32(aggregates, scratch, offsets,
            [](const SplitAggregate256& value) { return value.split_id; },
            stats.priority_passes, stats.priority_items);
  radix_u32(aggregates, scratch, offsets,
            [](const SplitAggregate256& value) { return value.first_edge; },
            stats.priority_passes, stats.priority_items);
  radix_u32(aggregates, scratch, offsets,
            [](const SplitAggregate256& value) { return value.first_tree; },
            stats.priority_passes, stats.priority_items);
  radix_u32(
      aggregates, scratch, offsets,
      [](const SplitAggregate256& value) {
        return static_cast<std::uint32_t>(~value.count);
      },
      stats.priority_passes, stats.priority_items);

  if (stats.priority_passes != kSplitPriorityRadixPasses) {
    throw std::logic_error("split-priority radix schedule mismatch");
  }
}

}  // namespace wqfm::internal
