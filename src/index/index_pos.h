/*
 * index_pos.h
 *
 *  Created on: Feb 5, 2017
 *      Author: isovic
 */

#ifndef SRC_INDEX_INDEX_POS_H_
#define SRC_INDEX_INDEX_POS_H_

#include <stdint.h>
#include "minimizer_index_consts.h"

namespace is {

class IndexPos {
 public:
  // This structure holds the info of a seed position in the form of (seq_id, pos, strand) 'tuple'.
  // The strand info is coded as the MSB of the pos. This way, when sorting the hits
  // All hits belonging to the same sequence will be grouped, but the fwd and the rev strand will be separate.
  uint32_t seq_id;     // ID of the originating sequence.
  uint32_t pos;   // Position within the sequence. If >= (1 << 31), the sequence is reverse complemented.

  IndexPos() : seq_id(0), pos(0) { }
  IndexPos(uint32_t nid, uint32_t npos) : seq_id(nid), pos(npos) { }
  IndexPos(uint64_t coded_pos) {
    seq_id = coded_pos >> 32;
    pos = coded_pos & (0x00000000FFFFFFFF);
  }
  bool is_rev() { return (pos & (kIndexIdReverse64)) != 0; }
  inline uint32_t get_pos() { return (pos & kIndexMaskStrand64); }
  inline uint32_t get_seq_id() { return seq_id; }
  static uint64_t CodePos(uint32_t ref_id, uint32_t pos, bool reverse) {
    return (((uint64_t) ref_id) << 32) | ((uint64_t) pos);
  }
};

} /* namespace is */

#endif /* SRC_INDEX_INDEX_POS_H_ */
