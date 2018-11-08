/*
 * seed.h
 *
 *  Created on: Mar 3, 2017
 *      Author: isovic
 */

#ifndef CODEBASE_GINDEX_SRC_MINIMIZER_INDEX_SEED_H_
#define CODEBASE_GINDEX_SRC_MINIMIZER_INDEX_SEED_H_

#include "minimizer_index_consts.h"

namespace is {

class Seed {
 public:
  Seed();
  ~Seed();

  /* Extract the position from a packed seed.
   */
  static inline int32_t seed_position(uint128_t seed) {
    return (int32_t) ((seed & kSeedMaskPos));
  }

  /* Extract the sequence ID from a packed seed.
   */
  static inline int32_t seed_seq_id(uint128_t seed) {
    return (int32_t) ((seed & kSeedMaskSeqId) >> 32);
  }

  /* Extract the key from a packed seed.
   */
  static inline int64_t seed_key(uint128_t seed) {
    return (int64_t) ((seed >> 64));
  }

  /* Create a packed seed from the three components.
   */
  static inline uint128_t make_seed(uint64_t key, uint64_t seq_id, uint64_t position) {
    return (((uint128_t) key) << 64) | (((uint128_t) seq_id) << 32)
        | (((uint128_t) position) << 0);
  }

};

} /* namespace is */

#endif /* CODEBASE_GINDEX_SRC_MINIMIZER_INDEX_SEED_H_ */
