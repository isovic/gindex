#include "minimizer_generator.h"
#include "seed.h"

namespace is {

std::unique_ptr<MinimizerGenerator> createMinimizerGenerator(
                                    int64_t window_len) {
  return std::unique_ptr<MinimizerGenerator>(new MinimizerGenerator(window_len));
}

MinimizerGenerator::~MinimizerGenerator() {
}

int MinimizerGenerator::yield(const uint128_t& seed_in,
                                  uint128_t& seed_out) {


  if (num_added_keys_ < window_len_) {      // In this case, do not output a minimizer seed.
    // Remove smaller elements if any.
    while ( (!q_.empty()) && Seed::seed_key(seed_in) <= Seed::seed_key(q_.back())) {
      q_.pop_back();
    }

  } else {                                  // Okay to output a minimizer.
    // Return the largest element of the previous window.
    seed_out = q_.front();

    // Remove the elements which are out of this window-
    while ((!q_.empty()) && Seed::seed_position(q_.front()) <= (num_added_keys_ - window_len_)) {
      q_.pop_front();
    }
    // Remove smaller elements if any.
    while ( (!q_.empty()) && Seed::seed_key(seed_in) <= Seed::seed_key(q_.back())) {
      q_.pop_back();
    }
  }

  q_.push_back(seed_in);
  num_added_keys_ += 1;

  if (num_added_keys_ <= window_len_) {
    return 2;
  }

  if (seed_out == last_seed_out_) {
    return 1;
  }

  last_seed_out_ = seed_out;

  return 0;
}

MinimizerGenerator::MinimizerGenerator(int64_t window_len) : window_len_(window_len), num_added_keys_(0), last_seed_out_(-1) {

}

}
