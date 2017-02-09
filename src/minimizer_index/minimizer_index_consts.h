/*
 * minimizer_index_consts.h
 *
 *  Created on: Feb 5, 2017
 *      Author: isovic
 */

#ifndef SRC_INDEX_MINIMIZER_INDEX_CONSTS_H_
#define SRC_INDEX_MINIMIZER_INDEX_CONSTS_H_

#include <stdint.h>

typedef unsigned __int128 uint128_t;

constexpr uint128_t kInvalidSeed = ((uint128_t) 0x0000000000000000FFFFFFFFFFFFFFFF) || (((uint128_t) 0x0000000000000000FFFFFFFFFFFFFFFF) << 64);
constexpr uint128_t kSeedMaskLow64 = ((uint128_t) 0x0000000000000000FFFFFFFFFFFFFFFF);
constexpr uint128_t kSeedMaskHigh64 = (((uint128_t) 0x0000000000000000FFFFFFFFFFFFFFFF) << 64);
constexpr uint128_t kSeedMaskPos = ((uint128_t) 0x000000000000000000000000FFFFFFFF);
constexpr uint128_t kSeedMaskSeqId = ((uint128_t) 0x0000000000000000FFFFFFFF00000000);



static const uint64_t empty_hash_key = 0xFFFFFFFFFFFFFFFF;
const uint64_t kIndexIdReverse128 = ((uint128_t) 1) << 31;
const uint64_t kIndexIdReverse64 = ((uint64_t) 1) << 31;
const uint64_t kIndexMaskStrand64 = ~((uint64_t) kIndexIdReverse64);
const uint64_t kIndexMaskLowerBits64 = 0x00000000FFFFFFFF;
const uint64_t kIndexMaskUpperBits64 = 0xFFFFFFFF00000000;

#define MAKE_HIT(key64, seqid32, pos32) ((((uint128_t) key64) << 64) | (((uint128_t) seqid32) << 32) | ((uint128_t) pos32))
#define GET_KEY_FROM_HIT(x)  ((uint64_t) ((x) >> 64))
/// #define GET_REAL_POS_FROM_HIT(x)  ((uint64_t) (((x) & ((uint128_t) 0x000000000000000000000000EFFFFFFF))))
// & kIndexMaskStrand64))
#define GET_SEQ_ID_FROM_HIT(x)  ((uint64_t) ((((x) & ((uint128_t) 0x0000000000000000FFFFFFFFFFFFFFFF)) >> 32) & kIndexMaskLowerBits64))
#define GET_CODED_SEED_FROM_HIT(x)  ((uint64_t) ((x) &0x0FFFFFFFFFFFFFFFF))
/// #define GET_POS_FROM_HIT_WITH_REV(x)  ((uint64_t) ((x) & kIndexMaskLowerBits64))

// #define GET_REAL_POS_FROM_HIT(x)  (((uint64_t) (x)) & ((uint64_t) 0x00000000EFFFFFFF))
#define GET_REAL_POS_FROM_HIT(x)  ((uint64_t) (x & 0x000000007FFFFFFF))
#define GET_POS_FROM_HIT_WITH_REV(x)  (((uint64_t) (x)) & ((uint64_t) 0x00000000FFFFFFFF))

#define IS_HIT_REVERSE(x) ((((uint128_t) (x) ) & (((uint128_t) 1) << 31)) != 0)



#endif /* SRC_INDEX_MINIMIZER_INDEX_CONSTS_H_ */
