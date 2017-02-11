/*
 * compiled_shape.h
 *
 *  Created on: Jun 24, 2016
 *      Author: isovic
 */

#ifndef SRC_MINIMIZERINDEX_COMPILED_SHAPE_H_
#define SRC_MINIMIZERINDEX_COMPILED_SHAPE_H_

#include <stdint.h>
#include <string>
#include <sstream>
#include <vector>

namespace is {

struct Mask {
  int32_t start = 0;
  int32_t len = 0;
  uint64_t bits = 0;
  int32_t shift = 0;
};

class CompiledShape {
 public:
  CompiledShape() : num_incl_bits_(0), max_width_(0) { };
  CompiledShape(const std::string new_shape) { Compile(new_shape); };

  /* Creates a compiled shape for a given string description of the shape.
   * The string is composed of characters '0' and '1', where '0' describes
   * a don't care base, and '1' an inclusive base.
   */
  void Compile(const std::string new_shape);

  /* Takes a buffer of bases (max 32 bases in 64-bits), 2bit packed, and extracts those inclusive ones (defined by a shape).
   * @buffer an integer containing 2-bit packed values of the sequence, max. 32 bases from the starting position.
   * @shape a string specifying the shape to be extracted from the buffer. Specified with '1' as inclusive bases and '0' as don't care bases.
   */
  uint64_t CreateSeedFromShape(uint64_t bases2bit) const;

  const std::vector<Mask>& masks() const {
    return masks_;
  }

  int32_t num_incl_bits() const {
    return num_incl_bits_;
  }

  const std::string& shape() const {
    return shape_;
  }

  const int32_t max_width() const {
    return max_width_;
  }

  std::string Verbose() const {
    std::stringstream ss;
    ss << "shape = '" << shape_ << "'" << std::endl;
    ss << "num_incl_bits = " << num_incl_bits_ << std::endl;
    ss << "max_width = " << max_width_ << std::endl;
    ss << "Masks:" << std::endl;
    for (int32_t i=0; i<masks_.size(); i++) {
      ss << "[mask " << i << "] start = " << masks_[i].start << ", len = " << masks_[i].len << ", bits = " << std::hex << masks_[i].bits << std::dec << ", shift = " << masks_[i].shift << std::endl;
    }
    return ss.str();
  }

 private:
  std::vector<Mask> masks_;
  std::string shape_;
  int32_t num_incl_bits_;
  int32_t max_width_;
};

std::vector<CompiledShape> CompileShapes(const std::vector<std::string> &shapes);

/* For a given shape, for every don't cate ('0') base generate all three combinations for the same seed,
 * containing at the position: 0 (deletion), 1 (match/mismatch) and 2 (insertion).
 * E.g. for a shape "1110111", this function will generate three shapes: "11111", "1110111" and "11100111".
 * The shapes vector is not cleared, so the function can be re-used to add multiple shapes.
 */
int CreateLookupShapes(std::string index_shape, std::vector<std::string> &shapes);

/* Converts a base 10 number x to a base N number with n digits.
 */
void Base10ToBaseN(int32_t x, int32_t N, int32_t n, std::vector<int8_t> &digits);

} /* Namespace is. */

#endif /* SRC_MINIMIZERINDEX_COMPILED_SHAPE_H_ */
