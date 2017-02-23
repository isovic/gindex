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

#include <stdint.h>
#include <map>
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include "sequences/sequence_file.h"
#include "utility/utility_general.h"
#include "compiled_shape.h"
#include "index_pos.h"
#include "minimizer_index_consts.h"

namespace is {

class IndexPos;
class MinimizerIndex;

struct SeedHashValue {
  int64_t start = 0, num = 0;
};

using google::dense_hash_map;
typedef dense_hash_map<uint64_t, SeedHashValue, std::hash<uint64_t> > SeedHashType;  // SeedHashType encodes the following: key is the hash key of a seed, and value is the position in the list of seeds where the seed starts as well as the number of positions with the same key.

/** A factory function for Minimizer Index. Uses a default seed hash function.
 */
std::shared_ptr<MinimizerIndex> createMinimizerIndex(const std::vector<std::string> &shapes, double freq_percentil);

/** A factory function for Minimizer Index. Loads the index from a given path.
 */
std::shared_ptr<MinimizerIndex> createMinimizerIndex(const std::string& path, double freq_percentil);

class MinimizerIndex {
 public:
  friend std::shared_ptr<MinimizerIndex> createMinimizerIndex(const std::vector<std::string> &shape, double freq_percentils);
  friend std::shared_ptr<MinimizerIndex> createMinimizerIndex(const std::string& path, double freq_percentil);

  ~MinimizerIndex();

  int Create(const SequenceFile &seqs,
             float min_avg_seed_qv, bool index_reverse_strand,
             bool use_minimizers, int32_t minimizer_window_len,
             int32_t num_threads, bool verbose=false);

  void CollectIndexSeeds(const int8_t *seqdata, const int8_t *seqqual, int64_t seqlen,
                    float min_avg_seed_q, bool index_reverse_strand,
                    bool use_minimizers, int32_t minimizer_window_len,
                    std::vector<uint128_t> &seed_list) const;

  void CollectLookupSeeds(const int8_t *seqdata, const int8_t *seqqual, int64_t seqlen,
                    float min_avg_seed_q, bool index_reverse_strand,
                    bool use_minimizers, int32_t minimizer_window_len,
                    std::vector<uint128_t> &seed_list) const;

  /* Writes the contents of the index to a file in binary format.
   * @path Path to a file where the index will be written to.
   * @return C-style return, 0 if everything went fine, 1 otherwisel
   */
  int Store(const std::string& path);

  /* Loads the index from a file.
   * @path Path to a file containing the index to load.
   * @return C-style return, 0 if everything went fine, 1 otherwisel
   */
  int Load(const std::string& path);

  /* For a given C-style seed, all lookup keys are compiled from lookup shapes, and queried.
   * Results are given as a vector of pointers to the corresponding hash buckets. This is done
   * in order to prevent potential copying of a large number of hits (for performance purposes).
   */
  int Find(const int8_t* seed, int32_t seed_len, bool threshold_hits, std::vector<const uint128_t *> &hits, std::vector<int64_t> &num_hits) const;

  /* For a given C-style seed, all lookup keys are compiled from lookup shapes, and queried.
   * All hits are copied to a new vector (unlike the alternative Find function) and self contained.
   */
  int FindAndJoin(const int8_t* seed, int32_t seed_len, bool threshold_hits, std::vector<uint128_t> &hits) const;

  void ConstructHashUnderCutoff();

  /* For a given seed, calculate keys for all lookup shapes.
   */
  void CalcKeysFromSeed(const int8_t *seed, int32_t seed_len, std::vector<uint64_t> &keys) const;

  /* Looks up all hits for a specific, calculated keys.
   */
  int KeyLookup(uint64_t seed_key, const uint128_t **seeds, int64_t *num_seeds) const;

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

  /* A debug function which writes the hash to a file on disk.
   * @out_path Path to a file on disk where data will be written to.
   * @num_bases The number of bases that a seed contains. Needed because seeds are stored as uint64 which can hold up to 32 bases.
   *            This number refers only to the inclusive bases, and not don't cares.
   */
  void DumpHash(std::string out_path, int32_t num_bases);
  void DumpHashSortedByCount(std::string out_path, int32_t num_bases);
  void DumpHashSortedByName(std::string out_path, int32_t num_bases);
  void DumpSeeds(std::string out_path, int32_t num_bases);


  /* Returns a pointer to the serialized data of all sequences making
   * up the index.
   */
  const std::vector<int8_t>& get_data() const {
    return data_serialized_;
  }

  const std::vector<std::string>& get_headers() const {
    return headers_;
  }

  const std::vector<int64_t>& get_reference_lengths() const {
    return reference_lengths_;
  }

  int64_t get_num_sequences_forward() const {
    return num_sequences_forward_;
  }

  int64_t get_num_sequences() const {
    return num_sequences_;
  }

  int64_t get_data_length() const {
    return data_length_;
  }

  int64_t get_data_length_forward() const {
    return data_length_forward_;
  }

  int64_t get_shape_max_width() const {
    return shape_max_width_;
  }

  const std::vector<int64_t>& get_reference_starting_pos() const {
    return reference_starting_pos_;
  }

  const std::vector<CompiledShape>& get_index_shapes() const {
    return index_shapes_;
  }

  const std::vector<CompiledShape>& get_lookup_shapes() const {
    return lookup_shapes_;
  }

  double count_cutoff() const {
    return count_cutoff_;
  }

  double avg_seed_occurrence() const {
    return avg_seed_occurrence_;
  }

  double max_seed_occurrence() const {
    return max_seed_occurrence_;
  }

  int32_t max_incl_bases() const {
    return max_incl_bases_;
  }

  ///////////////////////////////////////////////
  /// Legacy support. Needs updating.
  /// Also, needs implementing.
  int64_t RawPositionConverterWithRefId(int64_t raw_position, int64_t reference_index, int64_t query_length, int64_t *ret_absolute_position=NULL, int64_t *ret_relative_position=NULL, SeqOrientation *ret_orientation=NULL, int64_t *ret_reference_index_with_reverse=NULL) const;
  ///////////////////////////////////////////////

 private:
  explicit MinimizerIndex(const std::vector<std::string> &shapes, double freq_percentil);           // Private constructor, prevent memory leaks. Initializes an empty index for a given set of index shapes.
  explicit MinimizerIndex(const std::string& path, double freq_percentil);                          // Private constructor, prevent memory leaks. Loads the index from file.
  MinimizerIndex(const MinimizerIndex&) = delete;                   // No copy constructor.
  const MinimizerIndex& operator=(const MinimizerIndex&) = delete;  // No assignment operator.

  void Clear_();             // Clears internal storage when creating the index.
  void InitializeShapes_(const std::vector<std::string> &shapes);
  void ConstructHash_();

  int Serialize_(FILE *fp);

  int Deserialize_(FILE *fp);

  /* Calculates some info on input indexing shapes, such as the maximum seed len, maximum number of inclusive bits
   * (2-bit packed format) and the maximum width of a shape (should all don't care bases be replaced with the
   * potential maximum of 2 bases).
   */
  void CalcShapeStats_(const std::vector<CompiledShape>& index_shapes, int64_t &max_seed_len, int32_t &max_incl_bits, int32_t &shape_max_width, int32_t &max_inclusive_bases) const;

  /* Create shapes which will be used for lookup.
   * For every '0' base create 3 combinations: (mis)match, insertion, deletion.
   */
  std::vector<CompiledShape> CreateCompiledLookupShapes_(const std::vector<CompiledShape> &index_shapes);

  void AssignData_(const SequenceFile &seqs, bool index_reverse_strand);
  void AllocateSpaceForSeeds_(const SequenceFile &seqs,
                              bool index_reverse_strand, int64_t num_shapes,
                              int64_t max_seed_len, int64_t num_fwd_seqs,
                              std::vector<int64_t> &seed_starts_for_seq,
                              int64_t *total_num_seeds);

  static int CollectAllSeedsForSeq_(
      const int8_t *seqdata, const int8_t *seqqual, int64_t seqlen,
      float min_avg_seed_qv, bool index_reverse_strand,
      uint64_t seq_id_fwd, uint64_t seq_id_rev,
      const std::vector<CompiledShape> &compiled_shapes,
      uint128_t* seed_list_fwd, uint128_t* seed_list_rev);

  static int CollectAllSeedsForSeqWithQual_(
      const int8_t *seqdata, const int8_t *seqqual, int64_t seqlen,
      float min_avg_seed_qv, bool index_reverse_strand,
      uint64_t seq_id_fwd, uint64_t seq_id_rev,
      const std::vector<CompiledShape> &compiled_shapes,
      uint128_t* seed_list_fwd, uint128_t* seed_list_rev);

  static int CollectAllSeedsForSeqWithOutQual_(
      const int8_t *seqdata, int64_t seqlen,
      bool index_reverse_strand,
      uint64_t seq_id_fwd, uint64_t seq_id_rev,
      const std::vector<CompiledShape> &compiled_shapes,
      uint128_t* seed_list_fwd, uint128_t* seed_list_rev);

  static inline uint64_t SeedHashFunctionDefault_(uint64_t seed, int32_t k);

  static inline uint64_t ReverseComplementSeed_(uint64_t seed,
                                                int32_t num_bases);
  static int MakeMinimizers_(uint128_t *seed_list, int64_t num_seeds,
                             int64_t num_seeds_per_base, int32_t window_len);
  // Removes empty values, and shifts other seeds towards the top of the list.
  static int64_t MakeSeedListDense_(uint128_t *seed_list, int64_t num_seeds);

  int FlagDuplicates_(uint128_t *seed_list, int64_t num_seeds) const;
  int OccurrenceStatistics_(double percentil, int32_t num_threads,
                            double* ret_avg, double *ret_max, double* ret_stddev,
                            double *ret_percentil_val) const;

  /* Creates a vector of seeds for a given sequence. The seeds are packed in
   * the uint128_t format, where the MSB 64 bits specify the key, next 32 bits
   * the sequence ID and the 32 LSB bits the position of the seed on the sequence.
   * If specified, the seeds will be minimized as well.
   * @seqdata A C-style pointer to the raw sequence.
   * @seqqual A pointer to the quality values of the sequence, if available. Can be NULL to ommit.
   * @seqlen Length of the seqdata and seqqual arrays.
   * @min_avg_seed_qv A seed will be skipped if the QV's of its bases are below this threshold.
   * @seq_id_fwd The ID of the sequence in the forward direction. For example. Will be encoded in the seed.
   * @seq_id_rev When indexing a sequence, fwd and rev are separately indexed. Fwd are usually indexed first, and rev after that, which means that rev will usually have ID of (seq_id_fwd + num_sequences_forward_.
   * @seed_list A vector of all compiled seeds.
   */
  void CollectSeeds_(const int8_t* seqdata, const int8_t* seqqual, int64_t seqlen,
                     bool index_reverse_strand, uint64_t seq_id_fwd, uint64_t seq_id_rev,
                     int64_t max_seed_len, float min_avg_seed_qv, bool use_minimizers, int32_t minimizer_window_len,
                     const std::vector<CompiledShape> &index_shapes, std::vector<uint128_t> &seed_list) const;



  ///////////////////////////////////////////////
  /// Legacy support. Needs updating.
  /// This is the generic part of the index,
  /// and index independent.
  std::vector<int8_t> data_serialized_;             // TODO: remove data_ and replace with this!
  std::vector<int64_t> reference_lengths_;          // Initialized by AssignData_(...). Values for reverse complements are initialized as well.
  std::vector<int64_t> reference_starting_pos_;     // Initialized by AssignData_(...). Starting positions of each sequence in the serialized data_ array (absolute coordinates).
  std::vector<std::string> headers_;                // Initialized by AssignData_(...). Headers for rev complements are initialized as well.
  int64_t num_sequences_;                           // Initialized by AssignData_(...).
  int64_t num_sequences_forward_;                   // Initialized by AssignData_(...).
  int64_t data_length_;                             // Initialized by AssignData_(...).
  int64_t data_length_forward_;                     // Initialized by AssignData_(...).
  ///////////////////////////////////////////////

  ///////////////////////////////////////////////
  /// These are data for the index specific
  /// things
  std::vector<uint128_t> seeds_;                // A seed is encoded with: upper (MSB) 64 bits are the seed key, and lower (LSB) 64 bits are the ID of the sequence and the position (1-based, and encoded as in the IndexPos class). The position is 1-based to allow for undefined values.
  SeedHashType hash_;                           // A lookup for seeds by their key.

  // These are initialized.
  bool use_minimizers_;                         // Initialized by Create(...).
  int32_t minimizer_window_len_;                // Initialized by Create(...).
  double percentil_;                            // A value in range [0.0, 1.0] specifying the percentil of all seed occurrences to use as a cutoff threshold. E.g. percentil_ of 0.5 is a median.

  // These are calculated.
  double count_cutoff_;                         // Initialized by Create(...).
  double avg_seed_occurrence_;                  // Initialized by Create(...).
  double stddev_seed_occurrence_;
  double max_seed_occurrence_;

  // The following shape related values can all
  // be derived from index shapes in string form
  // by calling InitializeShapes_.
  std::vector<CompiledShape> index_shapes_;     // Initialized by the constructor.
  std::vector<CompiledShape> lookup_shapes_;    // Initialized by the constructor.
  int64_t max_seed_len_;                        // Initialized by the constructor.
  int32_t max_incl_bits_;                       // Initialized by the constructor.
  int32_t shape_max_width_;                     // Initialized by the constructor.
  int32_t max_incl_bases_;                      // Initialized by the constructor.

  ///////////////////////////////////////////////
  /// Constant values
  static const int64_t version_ = 12;
};

/** Helper functions for debugging.
 */
inline std::string SeedToString(uint64_t seed, int32_t num_bases) {
  std::stringstream ss;
  for (int32_t i = 0; i < num_bases; i++) {
    ss << "ACGT"[(seed & 0x03)];
    seed >>= 2;
  }
  std::string seed_string = ss.str();
  std::reverse(seed_string.begin(), seed_string.end());
  return seed_string;
}

inline bool invalid_seed(const uint128_t& seed) {
  return (seed == kInvalidSeed);
}

inline void set_invalid_seed(uint128_t& seed) {
  seed = kInvalidSeed;
}

std::string PrintSeed(uint128_t seed);

} /* Namespace is. */

/** TODO:
 * 1. Make minimizers uniformly sampled by - resetting the starting position for minimizer search at the position of the last found minimizer.
 * 2. Check if there is any usage of the 'default' value of 0 for a seed (e.g. such as in FlagDuplicates_). Positions are now 0-based, and this is not valid. MakeMinimizers_ uses it as well, line 413.
 * 3. Implement a different hash function.
 */

#endif /* SRC_OWLER2_INDEX_GAPPED_MINIMIZER_H_ */
