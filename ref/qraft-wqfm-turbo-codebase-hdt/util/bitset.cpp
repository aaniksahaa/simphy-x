#include "internal.hpp"

#include <bit>

namespace wqfm::internal {

/// Constructs a fixed-size zero-filled bitset.
Bitset::Bitset(int size) : size_(size), words_((size + 63) / 64, 0) {}

/// Returns true when no bit is set.
bool Bitset::empty() const {
  for (uint64_t word : words_) {
    if (word != 0) {
      return false;
    }
  }
  return true;
}

/// Sets one bit after validating the index.
void Bitset::set(int pos) {
  if (pos < 0 || pos >= size_) {
    throw std::out_of_range("bitset set out of range");
  }
  words_[pos / 64] |= (uint64_t{1} << (pos % 64));
}

/// Clears one bit after validating the index.
void Bitset::reset(int pos) {
  if (pos < 0 || pos >= size_) {
    throw std::out_of_range("bitset reset out of range");
  }
  words_[pos / 64] &= ~(uint64_t{1} << (pos % 64));
}

/// Tests one bit, returning false out of range.
bool Bitset::test(int pos) const {
  if (pos < 0 || pos >= size_) {
    return false;
  }
  return (words_[pos / 64] & (uint64_t{1} << (pos % 64))) != 0;
}

/// Counts set bits.
int Bitset::count() const {
  int total = 0;
  for (uint64_t word : words_) {
    total += static_cast<int>(std::popcount(word));
  }
  return total;
}

/// Clears all bits.
void Bitset::clear() {
  std::fill(words_.begin(), words_.end(), uint64_t{0});
}

/// Sets all valid bits.
void Bitset::fill() {
  std::fill(words_.begin(), words_.end(), ~uint64_t{0});
  mask_unused_bits();
}

/// Returns the complement within the fixed bit width.
Bitset Bitset::complement() const {
  Bitset result(size_);
  for (std::size_t i = 0; i < words_.size(); ++i) {
    result.words_[i] = ~words_[i];
  }
  result.mask_unused_bits();
  return result;
}

/// Returns all set bit positions.
std::vector<int> Bitset::positions() const {
  std::vector<int> result;
  result.reserve(count());
  for (int i = 0; i < size_; ++i) {
    if (test(i)) {
      result.push_back(i);
    }
  }
  return result;
}

/// Tests whether two same-sized bitsets intersect.
bool Bitset::intersects(const Bitset& other) const {
  if (size_ != other.size_) {
    throw std::invalid_argument("bitset size mismatch");
  }
  for (std::size_t i = 0; i < words_.size(); ++i) {
    if ((words_[i] & other.words_[i]) != 0) {
      return true;
    }
  }
  return false;
}

/// Tests whether the intersection is empty.
bool Bitset::intersection_empty(const Bitset& other) const {
  return !intersects(other);
}

/// Counts bits shared with another bitset.
int Bitset::intersection_count(const Bitset& other) const {
  if (size_ != other.size_) {
    throw std::invalid_argument("bitset size mismatch");
  }
  int total = 0;
  for (std::size_t i = 0; i < words_.size(); ++i) {
    total += static_cast<int>(std::popcount(words_[i] & other.words_[i]));
  }
  return total;
}

/// Returns the bitwise intersection.
Bitset Bitset::intersection(const Bitset& other) const {
  if (size_ != other.size_) {
    throw std::invalid_argument("bitset size mismatch");
  }
  Bitset result(size_);
  for (std::size_t i = 0; i < words_.size(); ++i) {
    result.words_[i] = words_[i] & other.words_[i];
  }
  return result;
}

/// Returns the bitwise union.
Bitset Bitset::unite(const Bitset& other) const {
  if (size_ != other.size_) {
    throw std::invalid_argument("bitset size mismatch");
  }
  Bitset result(size_);
  for (std::size_t i = 0; i < words_.size(); ++i) {
    result.words_[i] = words_[i] | other.words_[i];
  }
  return result;
}

/// Returns this set minus another set.
Bitset Bitset::difference(const Bitset& other) const {
  if (size_ != other.size_) {
    throw std::invalid_argument("bitset size mismatch");
  }
  Bitset result(size_);
  for (std::size_t i = 0; i < words_.size(); ++i) {
    result.words_[i] = words_[i] & ~other.words_[i];
  }
  result.mask_unused_bits();
  return result;
}

/// Compares bitset size and storage.
bool Bitset::operator==(const Bitset& other) const {
  return size_ == other.size_ && words_ == other.words_;
}

/// Orders bitsets deterministically for tie-breaking.
bool Bitset::lex_less(const Bitset& other) const {
  if (size_ != other.size_) {
    return size_ < other.size_;
  }
  for (int i = 0; i < size_; ++i) {
    const bool a = test(i);
    const bool b = other.test(i);
    if (a != b) {
      return !a && b;
    }
  }
  return false;
}

/// Hashes bitset size and words.
std::size_t Bitset::hash() const {
  uint64_t h = 1469598103934665603ULL ^ static_cast<uint64_t>(size_);
  for (uint64_t word : words_) {
    h ^= word;
    h *= 1099511628211ULL;
  }
  return static_cast<std::size_t>(h);
}

/// Clears padding bits in the last storage word.
void Bitset::mask_unused_bits() {
  if (size_ == 0 || words_.empty()) {
    return;
  }
  const int used = size_ % 64;
  if (used != 0) {
    words_.back() &= ((uint64_t{1} << used) - 1);
  }
}

}  // namespace wqfm::internal
