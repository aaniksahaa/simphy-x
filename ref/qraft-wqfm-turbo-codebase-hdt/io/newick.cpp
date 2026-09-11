#include "internal.hpp"

#include <cctype>

namespace wqfm {
namespace {

/// A simple recursive descent parser for Newick format.
class Parser {
public:
  /// Creates a parser over immutable Newick text.
  explicit Parser(std::string_view text) : text_(text) {}

  /// Parses an entire Newick document.
  Tree parse() {
    Tree tree;
    skip_ws_and_comments();
    const int root = parse_subtree(tree);
    tree.root = root;
    skip_ws_and_comments();
    if (peek() == ';') {
      ++pos_;
    }
    skip_ws_and_comments();
    if (pos_ != text_.size()) {
      throw std::invalid_argument("unexpected trailing Newick text");
    }
    return tree;
  }

private:
  std::string_view text_;
  std::size_t pos_ = 0;

  /// Returns the current character without consuming it.
  char peek() const { return pos_ < text_.size() ? text_[pos_] : '\0'; }

  /// Consumes and returns the current character.
  char get() {
    if (pos_ >= text_.size()) {
      throw std::invalid_argument("unexpected end of Newick text");
    }
    return text_[pos_++];
  }

  /// Skips whitespace and bracketed Newick comments.
  void skip_ws_and_comments() {
    bool progressed = true;
    while (progressed) {
      progressed = false;
      while (std::isspace(static_cast<unsigned char>(peek()))) {
        ++pos_;
        progressed = true;
      }
      if (peek() == '[') {
        ++pos_;
        int depth = 1;
        while (pos_ < text_.size() && depth > 0) {
          const char ch = text_[pos_++];
          if (ch == '[') {
            ++depth;
          } else if (ch == ']') {
            --depth;
          }
        }
        if (depth != 0) {
          throw std::invalid_argument("unterminated Newick comment");
        }
        progressed = true;
      }
    }
  }

  /// Appends a node to the parsed tree.
  int add_node(Tree &tree, std::string label = {}) {
    const int id = static_cast<int>(tree.nodes.size());
    tree.nodes.push_back(Tree::Node{{}, std::move(label)});
    return id;
  }

  /// Adds an undirected edge in the parsed tree.
  void add_edge(Tree &tree, int a, int b) {
    tree.nodes[a].neighbors.push_back(b);
    tree.nodes[b].neighbors.push_back(a);
  }

  /// Parses an internal subtree or labeled leaf.
  int parse_subtree(Tree &tree) {
    skip_ws_and_comments();
    if (peek() == '(') {
      get();
      const int node = add_node(tree);
      bool need_child = true;
      while (true) {
        skip_ws_and_comments();
        if (!need_child && peek() == ')') {
          break;
        }
        const int child = parse_subtree(tree);
        add_edge(tree, node, child);
        skip_ws_and_comments();
        if (peek() == ',') {
          get();
          need_child = true;
          continue;
        }
        if (peek() == ')') {
          break;
        }
        throw std::invalid_argument("expected ',' or ')' in Newick subtree");
      }
      get();
      parse_label();
      skip_branch_length();
      return node;
    }

    const std::string label = parse_label();
    if (label.empty()) {
      throw std::invalid_argument("leaf without label in Newick tree");
    }
    skip_branch_length();
    return add_node(tree, label);
  }

  /// Parses a quoted or unquoted Newick label.
  std::string parse_label() {
    skip_ws_and_comments();
    if (peek() == '\'') {
      get();
      std::string label;
      while (true) {
        const char ch = get();
        if (ch == '\'') {
          if (peek() == '\'') {
            get();
            label += '\'';
            continue;
          }
          break;
        }
        label += ch;
      }
      skip_ws_and_comments();
      return label;
    }

    std::string label;
    while (pos_ < text_.size()) {
      const char ch = peek();
      if (ch == ':' || ch == ',' || ch == ')' || ch == '(' || ch == ';' ||
          ch == '[' || ch == ']' ||
          std::isspace(static_cast<unsigned char>(ch))) {
        break;
      }
      label += ch;
      ++pos_;
    }
    skip_ws_and_comments();
    return label;
  }

  /// Skips an optional branch length token.
  void skip_branch_length() {
    skip_ws_and_comments();
    if (peek() != ':') {
      return;
    }
    ++pos_;
    while (pos_ < text_.size()) {
      const char ch = peek();
      if (ch == ',' || ch == ')' || ch == ';' || ch == '[') {
        break;
      }
      ++pos_;
    }
    skip_ws_and_comments();
  }
};

} // namespace

/// Parses Newick text into the public tree representation.
Tree parse_newick(const std::string &text) { return Parser(text).parse(); }

} // namespace wqfm
