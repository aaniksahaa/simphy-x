#include "wqfm/api.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if defined(_WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif

namespace {

/// Parsed command-line arguments for the CLI.
struct Args {
  std::string genes_path;
  std::string out_path;
  wqfm::Options options;
};

/// Prints CLI usage text.
void usage() {
  std::cerr << "Usage: wqfm-turbo --genes gene_trees.nwk [--out species_tree.nwk] "
                 "[--threads N] "
                 "[--algorithm hat-v1|hat-consensus|hat-consensus-nni-v1] "
                 "[--taxon-order sorted|java-hashset] "
                 "[--seed N] "
                 "[--consensus path] [--benchmark-json path] [--dump-consensus path] "
                 "[--debug-trace path] [--debug-max-depth N] [--debug-stop-after-root] "
                 "[--hat-split raw|consensus-balanced|consensus-polytomy-balanced] "
                 "[--hat-one-third-singleton] "
                 "[--no-hat-split-certificate] "
                 "[--hat-certificate-coefficient X] [--hat-certificate-candidates N] "
                 "[--hat-certificate-crossover K0] "
                 "[--hat-certificate-anchors N] "
                 "[--hat-one-x raw|normalized] "
                 "[--hat-norm-alpha X] "
                 "[--hat-nni-passes N] [--hat-nni-unbounded] "
                 "[--hat-nni-prefix-gate|--no-hat-nni-prefix-gate] "
                 "[--hat-nni-local-rescore|--no-hat-nni-local-rescore] "
                 "[--hat-nni-min-gain-per-gene X] "
                 "[--hat-nni-gate-split-drop-per-gene-v1 X] "
                 "[--hat-nni-pass-tree-dir path] "
                 "[--hat-structural-budget-factor X] "
                 "[--hat-debug-trace path] "
                 "[--progress|--no-progress]\n";
}

/// Parses the top-level algorithm mode.
wqfm::AlgorithmMode parse_algorithm(const std::string& value) {
  if (value == "hat-v1") return wqfm::AlgorithmMode::HatV1;
  if (value == "hat-consensus") return wqfm::AlgorithmMode::HatConsensus;
  if (value == "hat-consensus-nni-v1") {
    return wqfm::AlgorithmMode::HatConsensusNniV1;
  }
  throw std::invalid_argument(
      "invalid --algorithm: " + value +
      " (removed modes include native, attach-v1, and hat-v2; use hat-v1)");
}

/// Parses the taxon id assignment mode.
wqfm::TaxonOrderMode parse_taxon_order(const std::string& value) {
  if (value == "sorted") return wqfm::TaxonOrderMode::Sorted;
  if (value == "java-hashset") return wqfm::TaxonOrderMode::JavaHashSet;
  throw std::invalid_argument("invalid --taxon-order: " + value);
}

/// Parses the HAT active split objective option.
wqfm::HatSplitMode parse_hat_split(const std::string& value) {
  if (value == "raw") return wqfm::HatSplitMode::Raw;
  if (value == "consensus-balanced") {
    return wqfm::HatSplitMode::ConsensusBalanced;
  }
  if (value == "consensus-polytomy-balanced") {
    return wqfm::HatSplitMode::ConsensusPolytomyBalanced;
  }
  throw std::invalid_argument("invalid --hat-split: " + value);
}

/// Parses the HAT one-X objective option.
wqfm::HatOneXMode parse_hat_one_x(const std::string& value) {
  if (value == "raw") return wqfm::HatOneXMode::Raw;
  if (value == "normalized") return wqfm::HatOneXMode::Normalized;
  throw std::invalid_argument("invalid --hat-one-x: " + value);
}

/// Returns true when stderr is connected to an interactive terminal.
bool stderr_is_tty() {
#if defined(_WIN32)
  return _isatty(_fileno(stderr)) != 0;
#else
  return isatty(STDERR_FILENO) != 0;
#endif
}

/// Parses CLI flags and positional paths.
Args parse_args(int argc, char** argv) {
  Args args;
  args.options.progress = wqfm::ProgressMode::Auto;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto need_value = [&](const std::string& name) -> std::string {
      if (i + 1 >= argc) {
        throw std::invalid_argument("missing value for " + name);
      }
      return argv[++i];
    };
    if (arg == "--genes") {
      args.genes_path = need_value(arg);
    } else if (arg == "--out") {
      args.out_path = need_value(arg);
    } else if (arg == "--threads") {
      args.options.threads = std::stoi(need_value(arg));
    } else if (arg == "--algorithm") {
      args.options.algorithm = parse_algorithm(need_value(arg));
    } else if (arg == "--taxon-order") {
      args.options.taxon_order = parse_taxon_order(need_value(arg));
    } else if (arg == "--seed") {
      args.options.seed = std::stoull(need_value(arg));
    } else if (arg == "--consensus") {
      args.options.consensus_path = need_value(arg);
    } else if (arg == "--benchmark-json") {
      args.options.benchmark_json = need_value(arg);
    } else if (arg == "--dump-consensus") {
      args.options.dump_consensus = need_value(arg);
    } else if (arg == "--debug-trace") {
      args.options.debug_trace = need_value(arg);
    } else if (arg == "--debug-max-depth") {
      args.options.debug_max_depth = std::stoi(need_value(arg));
    } else if (arg == "--debug-stop-after-root") {
      args.options.debug_stop_after_root = true;
    } else if (arg == "--hat-split") {
      args.options.hat_split_mode = parse_hat_split(need_value(arg));
    } else if (arg == "--hat-one-third-singleton") {
      args.options.hat_one_third_singleton = true;
    } else if (arg == "--hat-split-certificate") {
      args.options.hat_split_certificate = true;
    } else if (arg == "--no-hat-split-certificate") {
      args.options.hat_split_certificate = false;
    } else if (arg == "--hat-certificate-coefficient") {
      args.options.hat_certificate_coefficient = std::stod(need_value(arg));
    } else if (arg == "--hat-certificate-crossover") {
      args.options.hat_certificate_crossover = std::stod(need_value(arg));
    } else if (arg == "--hat-certificate-candidates") {
      args.options.hat_certificate_candidates = std::stoi(need_value(arg));
    } else if (arg == "--hat-certificate-anchors") {
      args.options.hat_certificate_anchors = std::stoi(need_value(arg));
    } else if (arg == "--hat-one-x") {
      args.options.hat_one_x_mode = parse_hat_one_x(need_value(arg));
    } else if (arg == "--hat-norm-alpha") {
      args.options.hat_norm_alpha = std::stod(need_value(arg));
    } else if (arg == "--hat-nni-passes") {
      args.options.hat_nni_passes = std::stoi(need_value(arg));
      args.options.hat_nni_passes_explicit = true;
      args.options.hat_nni_unbounded_explicit = true;
      if (args.options.hat_nni_passes < 0) {
        args.options.hat_nni_unbounded = true;
      } else {
        args.options.hat_nni_unbounded = false;
      }
      if (!args.options.hat_nni_prefix_gate_explicit) {
        args.options.hat_nni_prefix_gate = false;
      }
    } else if (arg == "--hat-nni-unbounded") {
      args.options.hat_nni_unbounded = true;
      args.options.hat_nni_unbounded_explicit = true;
      if (!args.options.hat_nni_prefix_gate_explicit) {
        args.options.hat_nni_prefix_gate = false;
      }
    } else if (arg == "--no-hat-nni-unbounded") {
      args.options.hat_nni_unbounded = false;
      args.options.hat_nni_unbounded_explicit = true;
    } else if (arg == "--hat-nni-prefix-gate") {
      args.options.hat_nni_prefix_gate = true;
      args.options.hat_nni_prefix_gate_explicit = true;
    } else if (arg == "--no-hat-nni-prefix-gate") {
      args.options.hat_nni_prefix_gate = false;
      args.options.hat_nni_prefix_gate_explicit = true;
    } else if (arg == "--hat-nni-local-rescore") {
      args.options.hat_nni_local_rescore = true;
    } else if (arg == "--no-hat-nni-local-rescore") {
      args.options.hat_nni_local_rescore = false;
    } else if (arg == "--hat-nni-min-gain-per-gene") {
      args.options.hat_nni_min_gain_per_gene = std::stod(need_value(arg));
      if (args.options.hat_nni_min_gain_per_gene < 0.0) {
        throw std::invalid_argument(
            "--hat-nni-min-gain-per-gene must be non-negative");
      }
    } else if (arg == "--hat-nni-gate-low-first-gain-per-gene") {
      args.options.hat_nni_gate_low_first_gain_per_gene =
          std::stod(need_value(arg));
      args.options.hat_nni_prefix_gate = true;
      args.options.hat_nni_prefix_gate_explicit = true;
    } else if (arg == "--hat-nni-gate-high-first-gain-per-gene") {
      args.options.hat_nni_gate_high_first_gain_per_gene =
          std::stod(need_value(arg));
      args.options.hat_nni_prefix_gate = true;
      args.options.hat_nni_prefix_gate_explicit = true;
    } else if (arg == "--hat-nni-gate-high-gain-pass-cap") {
      args.options.hat_nni_gate_high_gain_pass_cap =
          std::stoi(need_value(arg));
      args.options.hat_nni_prefix_gate = true;
      args.options.hat_nni_prefix_gate_explicit = true;
    } else if (arg == "--hat-nni-gate-middle-min-gain-per-gene") {
      args.options.hat_nni_gate_middle_min_gain_per_gene =
          std::stod(need_value(arg));
      args.options.hat_nni_prefix_gate = true;
      args.options.hat_nni_prefix_gate_explicit = true;
    } else if (arg == "--hat-nni-gate-middle-rel-max-gain") {
      args.options.hat_nni_gate_middle_rel_max_gain =
          std::stod(need_value(arg));
      args.options.hat_nni_prefix_gate = true;
      args.options.hat_nni_prefix_gate_explicit = true;
    } else if (arg == "--hat-nni-gate-middle-min-passes") {
      args.options.hat_nni_gate_middle_min_passes =
          std::stoi(need_value(arg));
      args.options.hat_nni_prefix_gate = true;
      args.options.hat_nni_prefix_gate_explicit = true;
    } else if (arg == "--hat-nni-gate-max-passes") {
      args.options.hat_nni_gate_max_passes = std::stoi(need_value(arg));
      args.options.hat_nni_prefix_gate = true;
      args.options.hat_nni_prefix_gate_explicit = true;
    } else if (arg == "--hat-nni-gate-split-drop-v1") {
      args.options.hat_nni_gate_split_drop_v1 = std::stoi(need_value(arg));
      args.options.hat_nni_prefix_gate = true;
      args.options.hat_nni_prefix_gate_explicit = true;
    } else if (arg == "--hat-nni-gate-split-drop-per-gene-v1") {
      args.options.hat_nni_gate_split_drop_per_gene_v1 =
          std::stod(need_value(arg));
      args.options.hat_nni_prefix_gate = true;
      args.options.hat_nni_prefix_gate_explicit = true;
    } else if (arg == "--hat-nni-pass-tree-dir") {
      args.options.hat_nni_pass_tree_dir = need_value(arg);
    } else if (arg == "--hat-structural-budget-factor") {
      args.options.hat_structural_budget_factor = std::stod(need_value(arg));
      if (args.options.hat_structural_budget_factor < 1.0) {
        throw std::invalid_argument(
            "--hat-structural-budget-factor must be at least 1.0");
      }
    } else if (arg == "--hat-debug-trace") {
      args.options.hat_debug_trace = need_value(arg);
    } else if (arg == "--progress") {
      args.options.progress = wqfm::ProgressMode::On;
    } else if (arg == "--no-progress") {
      args.options.progress = wqfm::ProgressMode::Off;
    } else if (arg == "--help" || arg == "-h") {
      usage();
      std::exit(0);
    } else {
      throw std::invalid_argument("unknown argument: " + arg);
    }
  }
  if (args.genes_path.empty()) {
    throw std::invalid_argument("--genes is required");
  }
  if (args.options.progress == wqfm::ProgressMode::Auto) {
    args.options.progress = stderr_is_tty() ? wqfm::ProgressMode::On
                                            : wqfm::ProgressMode::Off;
  }
  return args;
}

/// Tests whether a line contains only whitespace.
bool blank(const std::string& line) {
  for (char ch : line) {
    if (!std::isspace(static_cast<unsigned char>(ch))) {
      return false;
    }
  }
  return true;
}

}  // namespace

/// CLI entry point.
int main(int argc, char** argv) {
  try {
    Args args = parse_args(argc, argv);
    const auto parse_started = std::chrono::steady_clock::now();
    std::ifstream in(args.genes_path);
    if (!in) {
      throw std::runtime_error("failed to open genes file: " + args.genes_path);
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
      if (!blank(line)) {
        lines.push_back(std::move(line));
      }
    }
    in.close();
    if (lines.empty()) {
      throw std::runtime_error("genes file did not contain any Newick trees");
    }
    if (lines.size() >
        static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::length_error("gene count exceeds the integer-indexed domain");
    }
    // Parse Newick trees in parallel: each line is independent and parse_newick
    // is pure, so we fill genes[i] from lines[i] (input order preserved). A
    // local work-stealing pool keeps this independent of the solver executor,
    // which is configured later inside infer_species_tree.
    std::vector<wqfm::GeneTree> genes(lines.size());
    {
      int nw = args.options.threads > 0
                   ? args.options.threads
                   : static_cast<int>(std::thread::hardware_concurrency());
      if (nw < 1) nw = 1;
      nw = std::min<int>(nw, static_cast<int>(lines.size()));
      if (nw <= 1) {
        for (std::size_t i = 0; i < lines.size(); ++i) {
          genes[i] = wqfm::parse_newick(lines[i]);
        }
      } else {
        std::atomic<std::size_t> next{0};
        std::vector<std::exception_ptr> failures(lines.size());
        // jthread owns the unwind path: if a later worker launch throws, every
        // already-created worker is joined before its captured stack state dies.
        std::vector<std::jthread> pool;
        pool.reserve(static_cast<std::size_t>(nw));
        for (int t = 0; t < nw; ++t) {
          pool.emplace_back([&]() {
            std::size_t i;
            while ((i = next.fetch_add(1, std::memory_order_relaxed)) <
                   lines.size()) {
              try {
                genes[i] = wqfm::parse_newick(lines[i]);
              } catch (...) {
                failures[i] = std::current_exception();
              }
            }
          });
        }
        for (std::jthread& th : pool) {
          th.join();
        }
        // Input order is semantic. Report the lowest failing line regardless
        // of which worker happened to complete first.
        for (const std::exception_ptr& failure : failures) {
          if (failure != nullptr) std::rethrow_exception(failure);
        }
      }
    }
    args.options.parse_input_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - parse_started)
            .count();
    std::vector<std::string>().swap(lines);

    const wqfm::SpeciesTree species = wqfm::infer_species_tree(std::move(genes), args.options);
    const std::string output = wqfm::to_newick(species) + "\n";
    if (args.out_path.empty()) {
      std::cout << output;
    } else {
      std::ofstream out(args.out_path);
      if (!out) {
        throw std::runtime_error("failed to open output file: " + args.out_path);
      }
      out << output;
    }
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "wqfm-turbo: " << ex.what() << "\n";
    usage();
    return 1;
  }
}
