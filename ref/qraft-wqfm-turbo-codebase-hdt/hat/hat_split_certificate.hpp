#pragma once

#include "internal.hpp"

namespace wqfm::internal {

// Immutable, linear-space index of the ORIGINAL, degree-two-suppressed genes.
// A compressed child bank does not retain these graph distances.
struct HatGeneDistanceIndex {
  std::vector<int> depth;
  std::vector<int> head;

  int distance(const GeneIndex& gene, int a, int b) const;
};

struct HatSplitCertificateIndex {
  bool complete_single_copy = false;
  int taxon_count = 0;
  std::vector<HatGeneDistanceIndex> genes;
};

struct HatSplitCertificateResult {
  bool available = false;
  bool accepted = false;
  double gap = 0.0;
  double threshold = 0.0;
};

std::shared_ptr<const HatSplitCertificateIndex> make_hat_split_certificate_index(
    const std::vector<GeneIndex>& genes, int taxon_count);

double hat_split_certificate_threshold(int gene_count, int taxon_count,
                                       double coefficient, double crossover);

HatSplitCertificateResult hat_check_split_certificate(
    const SolverContext& ctx, const std::vector<int>& active_taxa,
    const Bitset& split);

}  // namespace wqfm::internal
