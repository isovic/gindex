#ifndef MINIMIZER_GENERATOR_H_
#define MINIMIZER_GENERATOR_H_

#include <stdint.h>
#include <deque>
#include <memory>

#include "minimizer_index_consts.h"

namespace is {
class MinimizerGenerator;

std::unique_ptr<MinimizerGenerator> createMinimizerGenerator(int64_t window_len);

class MinimizerGenerator {
 public:
  friend std::unique_ptr<MinimizerGenerator> createMinimizerGenerator(int64_t window_len);
  ~MinimizerGenerator();

  /* Pushes the new seed into the queue and returns the minimum previous value
   * within the given window.
   * Returns 1 if the queue has not yet been filled up to window_len size,
   * otherwise returns 0.
   */
  int yield(const uint128_t& seed_in, uint128_t& seed_out);

 private:
  MinimizerGenerator(int64_t window_len);

  int64_t window_len_;
  int64_t num_added_keys_;
  uint128_t last_seed_out_;
  std::deque<uint128_t> q_;
};

}

#endif
