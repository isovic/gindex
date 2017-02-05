/*
 * index_gapped_minimizer.h
 *
 *  Created on: Jun 24, 2016
 *      Author: isovic
 */

#ifndef SRC_OWLER2_INDEX_GAPPED_MINIMIZER_H_
#define SRC_OWLER2_INDEX_GAPPED_MINIMIZER_H_

#include "sparsehash/dense_hash_map"

//typedef dense_hash_map<int64_t, DenseType2, std::hash<int64_t> > DenseType;
//std::vector<DenseType> bins_map1;
//bins_map1[i].set_empty_key(-1);
//struct DenseType2 {
//  int32_t timestamp = 0;
//  float count = 0.0;
//};
//DenseType2 &hit = temp_map[position_bin];
//if (hit.timestamp == (i + 1)) { continue; }
//hit.count += 1.0f;
//hit.timestamp = (i + 1);

#include <vector>
#include <map>
#include <stdint.h>
#include <stdint.h>
#include <memory>
#include "sequences/sequence_file.h"
#include "utility/utility_general.h"
#include "compiled_shape.h"
#include "index_pos.h"
#include "minimizer_index_consts.h"

namespace is {

using google::dense_hash_map;      // namespace where class lives by default

struct SeedHashValue {
  int64_t start = 0, num = 0;
};
typedef dense_hash_map<uint64_t, SeedHashValue, std::hash<uint64_t> > SeedHashType;     // SeedHashType encodes the following: key is the hash key of a seed, and value is the position in the list of seeds where the seed starts as well as the number of positions with the same key.

class IndexPos;
class MinimizerIndex;

std::shared_ptr<MinimizerIndex> createMinimizerIndex();

class MinimizerIndex {
 public:
  friend std::shared_ptr<MinimizerIndex> createMinimizerIndex();
  ~MinimizerIndex();

  int Create(const SequenceFile &seqs, const std::vector<std::string> &shapes,
                             float min_avg_seed_qv, bool index_reverse_strand, bool use_minimizers,
                             int32_t minimizer_window_len, int32_t num_threads);

  int GetHits(uint64_t seed_key, uint128_t **seeds, int64_t *num_seeds) const;

  static void CollectMinimizers(const int8_t *seqdata, const int8_t *seqqual, int64_t seqlen,
                                float min_avg_seed_qv, uint64_t seq_id_fwd, uint64_t seq_id_rev,
                                const std::vector<CompiledShape> &compiled_shapes, int64_t minimizer_window_len,
                                std::vector<uint128_t> &seed_list, int64_t *num_minimizers);

  /* A debug function which writes the hash to a file on disk.
   * @out_path Path to a file on disk where data will be written to.
   * @num_bases The number of bases that a seed contains. Needed because seeds are stored as uint64 which can hold up to 32 bases. This number refers only to the inclusive bases, and not don't cares.
   */
  void DumpHash(std::string out_path, int32_t num_bases);
  void DumpHashSortedByCount(std::string out_path, int32_t num_bases);
  void DumpHashSortedByName(std::string out_path, int32_t num_bases);
  void DumpSeeds(std::string out_path, int32_t num_bases);

 private:
  std::vector<uint128_t> seeds_;      // A seed is encoded with: upper (MSB) 64 bits are the seed key, and lower (LSB) 64 bits are the ID of the sequence and the position (1-based, and encoded as in the IndexPos class). The position is 1-based to allow for undefined values.
  SeedHashType hash_;                 // A lookup for seeds by their key.
  std::vector<std::vector<int8_t> > data_;      // Actual sequences which have been indexed. The outer vector contains both fwd and revcmp sequences (fwd come first).

  int64_t num_sequences_;
  int64_t num_sequences_forward_;
  std::vector<int64_t> reference_lengths_;
  std::vector<std::string> headers_;
//  // Seeds are extracted for each sequence separately, but are stored in a giant array. Each sequence 'i' is designated to belong to a part of that array, starting with seq_seed_starts_[i] position in seed_list_.
//  std::vector<int64_t> seq_seed_starts_;
//  std::vector<int64_t> seq_seed_counts_;
  double count_cutoff_;
  std::vector<CompiledShape> lookup_shapes_;

  MinimizerIndex();
  MinimizerIndex(const MinimizerIndex&) = delete;
  const MinimizerIndex& operator=(const MinimizerIndex&) = delete;

  void Clear_();
  void AssignData_(const SequenceFile &seqs, bool index_reverse_strand);
  void AllocateSpaceForSeeds_(const SequenceFile &seqs, bool index_reverse_strand, int64_t num_shapes, int64_t max_seed_len, int64_t num_fwd_seqs, std::vector<int64_t> &seed_starts_for_seq, int64_t *total_num_seeds);

  static int CollectSeedsForSeq_(const int8_t *seqdata, const int8_t *seqqual, int64_t seqlen, float min_avg_seed_qv, uint64_t seq_id_fwd, uint64_t seq_id_rev, const std::vector<CompiledShape> &compiled_shapes, uint128_t* seed_list_fwd, uint128_t* seed_list_rev);
  static inline uint64_t SeedHashFunction_(uint64_t seed);
  static inline uint64_t ReverseComplementSeed_(uint64_t seed, int32_t num_bases);
  static int MakeMinimizers_(uint128_t *seed_list, int64_t num_seeds, int64_t num_seeds_per_base, int32_t window_len);
  // Removes empty values, and shifts other seeds towards the top of the list.
  static int MakeSeedListDense_(uint128_t *seed_list, int64_t num_seeds);

  int FlagDuplicates_(uint128_t *seed_list, int64_t num_seeds) const;
  int OccurrenceStatistics_(double percentil, int32_t num_threads, double* ret_avg, double* ret_stddev, double *ret_percentil_val) const;
};

// Helper functions for debugging.

inline std::string SeedToString(uint64_t seed, int32_t num_bases) {
  std::stringstream ss;
  for (int32_t i=0; i<num_bases; i++) {
    ss << "ACGT"[(seed & 0x03)];
    seed >>= 2;
  }
  std::string seed_string = ss.str();
  std::reverse(seed_string.begin(), seed_string.end());
  return seed_string;
}

inline uint128_t make_seed(uint64_t key, uint64_t seq_id, uint64_t position) {
  return (((uint128_t) key) << 64) | (((uint128_t) seq_id) << 32) | (((uint128_t) position) << 0);
}

} /* Namespace is. */

/** TODO:
 * 1. Make minimizers uniformly sampled by - resetting the starting position for minimizer search at the position of the last found minimizer.
 * 2. Check if there is any usage of the 'default' value of 0 for a seed (e.g. such as in FlagDuplicates_). Positions are now 0-based, and this is not valid. MakeMinimizers_ uses it as well, line 413.
 * 3. Implement a different hash function.
 */

#endif /* SRC_OWLER2_INDEX_GAPPED_MINIMIZER_H_ */
