/*
 * compiled_shape.cc
 *
 *  Created on: Jun 24, 2016
 *      Author: isovic
 */

#include "compiled_shape.h"
#include <math.h>
#include <sstream>
#include <algorithm>

namespace is {

void CompiledShape::Compile(const std::string new_shape) {
  shape_ = new_shape;
  masks_.clear();
  num_incl_bits_ = 0;

  Mask mask;
  for (int32_t i=0; i<new_shape.size(); i++) {
    if (shape_[i] == '1' && (i == 0 || (i > 0 && shape_[i-1] == '0'))) {  /// This detects the start of a new chain of inclusive bases.
      mask.start = i*2;
      mask.len = 2; /// Counts two bits per base.
    } else if (shape_[i] == '0' && (i > 0 && shape_[i-1] == '1')) {       /// This detects the end of the chain of inclusive bases.
      masks_.push_back(mask);
      mask.start = 0;
      mask.len = 0;
    } else if (shape_[i] == '1') { /// Count the number of inclusive bases.
      mask.len += 2; /// Counts two bits per base.
    }
  }

  /// Leftover mask part.
  if (mask.len > 0) {
    masks_.push_back(mask);
    mask.start = 0;
    mask.len = 0;
  }

  max_width_ = 0;
  for (int32_t i = 0; i < new_shape.size(); i++) {
    max_width_ += ((new_shape[i] == '1') ? 1 : 2);  /// '0' can also mean an insertion, so it can occupy two bases instead of one.
  }

//    shape = "11110111101111";
//    Mask mask;
//    mask.start = 0; mask.len = 4*2; masks.push_back(mask);
//    mask.start = 5*2; mask.len = 4*2; masks.push_back(mask);
//    mask.start = 10*2; mask.len = 4*2; masks.push_back(mask);

  for (int32_t i=0; i<masks_.size(); i++) {
    int64_t shift_mask_raw = 64 - masks_[i].start - masks_[i].len;
    masks_[i].bits = ((uint64_t) (pow(2, masks_[i].len) - 1)) << (shift_mask_raw);
    masks_[i].shift = (i > 0) ? (masks_[i].start - (masks_[i-1].start + masks_[i-1].len) + masks_[i-1].shift) : 0;
    num_incl_bits_ += masks_[i].len;
//      printf ("%s\t\t[%d] start = %d, len = %d, shift = %d, shape_incl_bases = %d\n", ConvertToBinary(masks[i].bits).c_str(), i, masks[i].start, masks[i].len, masks[i].shift, num_incl_bases);
  }
}

uint64_t CompiledShape::CreateSeedFromShape(uint64_t bases2bit) const {
  uint64_t seed = 0;
  for (int32_t i=0; i<masks_.size(); i++) {
    uint64_t buffer_part = bases2bit & masks_[i].bits;
    buffer_part = buffer_part << masks_[i].shift;
    seed |= buffer_part;
  }
  seed = seed >> (64 - num_incl_bits_);
  return seed;
}

std::vector<CompiledShape> CompileShapes(const std::vector<std::string> &shapes) {
  std::vector<std::string> lshapes = shapes;
  std::sort(lshapes.begin(), lshapes.end());
  std::vector<CompiledShape> compiled_shapes;
  for (int32_t i=0; i<lshapes.size(); i++) {
    if (i > 0 && lshapes[i-1] == lshapes[i]) { continue; }      // Remove duplicates.
    CompiledShape compiled_shape(lshapes[i]);
    compiled_shapes.push_back(compiled_shape);
  }
  return compiled_shapes;
}

void Base10ToBaseN(int32_t x, int32_t N, int32_t n, std::vector<int8_t> &digits) {
  digits.clear();

  while (n > 0) {
    digits.push_back(x % N);
    x /= N;
    n -= 1;
  }
  std::reverse(digits.begin(), digits.end());
}

int CreateLookupShapes(std::string index_shape, std::vector<std::string> &shapes) {
  int32_t n = 0;
  for (int32_t i=0; i<index_shape.size(); i++) {
    // Sanity check.
    if (index_shape[i] != '0' && index_shape[i] != '1') { return 1; }
    if (index_shape[i] == '0') { n += 1; }
  }

  int32_t m = pow(3, n);

  for (int32_t i=0; i<m; i++) {
    std::vector<int8_t> digits;
    Base10ToBaseN(i, 3, n, digits);
    int32_t num_digits = digits.size() - 1;

    std::stringstream ss;
    for (int32_t j=0; j<index_shape.size(); j++) {
      if (index_shape[j] == '1') {
        ss << "1";
      } else {
        ss << std::string(digits[num_digits], '0');
        num_digits -= 1;
      }
    }
    shapes.push_back(ss.str());
  }

  return 0;
}

} /* Namespace is. */
