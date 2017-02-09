/*
 * index_gapped_minimizer.cc
 *
 *  Created on: Jun 24, 2016
 *      Author: isovic
 */

#include "minimizer_index.h"

#include "log_system/log_system.h"
#include "omp_sort.hpp"
#include <omp.h>
#include <algorithm>
#include <tuple>
#include <deque>

namespace is {

std::shared_ptr<MinimizerIndex> createMinimizerIndex(const std::vector<std::string> &shapes) {
  return std::shared_ptr<MinimizerIndex>(new MinimizerIndex(shapes));
}

MinimizerIndex::MinimizerIndex(const std::vector<std::string> &shapes) : percentil_(0.999) {
  hash_.set_empty_key(empty_hash_key);

  Clear_();
  // Compile the given shapes.
  index_shapes_ = CompileShapes(shapes);
  lookup_shapes_ = CreateCompiledLookupShapes_(index_shapes_);
  CalcShapeStats_(index_shapes_, max_seed_len_, max_incl_bits_, shape_max_width_);
}

MinimizerIndex::~MinimizerIndex() {
  Clear_();
}

void MinimizerIndex::Clear_() {
  count_cutoff_ = 0.0;
  num_sequences_ = 0;
  num_sequences_forward_ = 0;
  seeds_.clear();
  hash_.clear();
//  data_.clear();
  data_serialized_.clear();
  reference_lengths_.clear();
  headers_.clear();

  // These should not be cleared, they are defined by
  // the constructor, and no other parameter should be
  // able to change these values.
  //  lookup_shapes_.clear();
  //  max_seed_len_ = 0;
  //  max_incl_bits_ = 0;
  //  shape_max_width_ = 0;

  use_minimizers_ = false;
  minimizer_window_len_ = 1;
}

std::vector<CompiledShape> MinimizerIndex::CreateCompiledLookupShapes_(const std::vector<CompiledShape> &index_shapes) {
  std::vector<std::string> ls;
  for (int32_t i=0; i<index_shapes.size(); i++) {
    CreateLookupShapes(index_shapes[i].shape(), ls);
  }

  std::vector<CompiledShape> lookup_shapes;
  lookup_shapes = CompileShapes(ls);
  return lookup_shapes;
}

void MinimizerIndex::CalcShapeStats_(const std::vector<CompiledShape>& index_shapes, int64_t &max_seed_len, int32_t &max_incl_bits, int32_t &shape_max_width) const {
  // Determine the maximum length of all shapes, to limit the last kmer being checked.
  max_seed_len = 0;
  max_incl_bits = 0;
  shape_max_width = 0;
  for (int32_t i=0; i<index_shapes_.size(); i++) {
    max_seed_len = std::max(max_seed_len, (int64_t) index_shapes[i].shape().size());
    max_incl_bits = std::max(max_incl_bits, index_shapes[i].num_incl_bits());
    shape_max_width  = std::max(shape_max_width, index_shapes[i].max_width());
  }
}

int MinimizerIndex::Create(const SequenceFile& seqs, float min_avg_seed_qv, bool index_reverse_strand, bool use_minimizers, int32_t minimizer_window_len, int32_t num_threads, bool verbose) {
  Clear_();

  clock_t absolute_time = clock();
  clock_t diff_time = clock();

  use_minimizers_ = use_minimizers;
  minimizer_window_len_ = minimizer_window_len;
  int64_t n_seqs = seqs.get_sequences().size();

  AssignData_(seqs, index_reverse_strand);

  int64_t total_num_seeds = 0;
  std::vector<int64_t> seed_starts_for_seq;

  AllocateSpaceForSeeds_(seqs, index_reverse_strand, index_shapes_.size(), max_seed_len_, num_sequences_forward_, seed_starts_for_seq, &total_num_seeds);

  if (verbose) {
    LOG_ALL("Allocated memory for a list of %ld seeds (128 bits each) (%.5f sec, diff: %.5f sec).\n", seeds_.size(), (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC));
    LOG_ALL("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str());
    LOG_ALL("Collecting seeds with %ld threads.\n", num_threads);
    if (use_minimizers) {
      LOG_ALL("Minimizer seeds will be used. Minimizer window is %ld.\n", minimizer_window_len);
    }
  }

  diff_time = clock();

  // Process all reads in parallel.
  #pragma omp parallel for num_threads(num_threads) shared(seqs) schedule(dynamic, 1)
  for (int64_t i=0; i<seqs.get_sequences().size(); i++) {
    uint32_t thread_id = omp_get_thread_num();
    if (thread_id == 0 && verbose) {
      LOG_ALL("\rRead %ld/%ld", (i + 1), seqs.get_sequences().size());
    }

    int8_t *seqdata = (int8_t *) seqs.get_sequences()[i]->get_data();
    int8_t *seqqual = (int8_t *) seqs.get_sequences()[i]->get_quality();
    int64_t seqlen = seqs.get_sequences()[i]->get_data_length();

    uint64_t seq_id = seqs.get_sequences()[i]->get_sequence_absolute_id() + 0; // The '+ 0' is actually short for this: + ((orientation == kReverse) ? num_sequences_forward_ : 0);

    uint128_t *seeds_fwd_ptr = &(seeds_[seed_starts_for_seq[i]]);
    uint128_t *seeds_rev_ptr = NULL;

    if (index_reverse_strand) {
      seeds_rev_ptr = &(seeds_[seed_starts_for_seq[n_seqs + i]]);
    }

    // Collect all seeds.
    int64_t num_seeds_processed = CollectAllSeedsForSeq_(seqdata, seqqual, seqlen, min_avg_seed_qv,
                                                         index_reverse_strand,
                                                         seq_id, seq_id + n_seqs,
                                                         index_shapes_,
                                                         seeds_fwd_ptr,
                                                         seeds_rev_ptr);

    // Create minimizers if needed.
    if (use_minimizers) {
//      MakeMinimizers_(&(seeds_[seed_starts_for_seq[i]]), num_seeds_processed, 2, minimizer_window_len);     // 2 only refers to the fwd and rev complement (2 seeds per base).
      MakeMinimizers_(seeds_fwd_ptr, num_seeds_processed, 1, minimizer_window_len);     // 2 only refers to the fwd and rev complement (2 seeds per base).
      if (index_reverse_strand) {
        MakeMinimizers_(seeds_rev_ptr, num_seeds_processed, 1, minimizer_window_len);     // 2 only refers to the fwd and rev complement (2 seeds per base).
      }
    }
  }

  if (verbose) { LOG_NEWLINE; }

//  DumpSeeds("temp/seeds.sparse.minimizers.csv", max_incl_bits/2);

  if (use_minimizers) {
    // Remove all excess seeds so that sorting will be faster.
    if (verbose) {
      LOG_ALL("Removing excess seeds.\n");
    }
    int64_t num_dense_seeds = MakeSeedListDense_(&(seeds_[0]), seeds_.size());
    seeds_.resize(num_dense_seeds);
  }

//  DumpSeeds("temp/seeds.dense.minimizers.csv", max_incl_bits/2);

  if (verbose) {
    LOG_ALL("Sorting the seeds (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC));
  }

  diff_time = clock();
  pquickSort(&(seeds_[0]), seeds_.size(), num_threads);
  //    FlagDuplicates_(&(seeds_[seed_starts_for_seq[i]]), num_seeds_processed);
//  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Making unique (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
//  diff_time = clock();

  if (verbose) {
    LOG_ALL("Generating the hash table (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC));
  }

  diff_time = clock();

  // This part generates the lookup table for each seed key.
  uint64_t prev_key = 0;
  bool init = false;
  SeedHashValue new_hash_val;
  for (int64_t i=0; i<seeds_.size(); i++) {
    uint64_t key = (uint64_t) (seeds_[i] >> 64);
    if (init == false) {
      new_hash_val.start = i;
      new_hash_val.num = 1;
      init = true;
    } else if (key != prev_key && new_hash_val.num > 0) {
      hash_[prev_key] = new_hash_val;
      new_hash_val.start = i;
      new_hash_val.num = 1;
    } else {
      new_hash_val.num += 1;
    }
    prev_key = key;
  }
  // Store the last hash key in the list.
  if (new_hash_val.num > 0) {
    hash_[prev_key] = new_hash_val;
    new_hash_val.start = 0;
    new_hash_val.num = 0;
  }

  // Calculate a cutoff threshold as a percentil of the occurrence of a seed.
  double avg = 0.0f, stddev = 0.0f, percentil_value = 0.0f;
  if (verbose) {
    LOG_ALL("Calculating the distribution statistics for key counts (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC));
  }
  OccurrenceStatistics_(percentil_, num_threads, &avg, &stddev, &percentil_value);
  count_cutoff_ = percentil_value;
  if (verbose) {
    LOG_ALL("Index statistics: average key count = %f, std dev = %f, percentil(%.2f%%) = %f\n", avg, stddev, percentil_*100.0, percentil_value);
  }

//  DumpSeeds("temp/seeds.csv", max_incl_bits_/2);
//  DumpHash("temp/hash.csv", max_incl_bits_/2);
//  DumpHashSortedByCount("temp/hash.minimizers.sorted.csv", max_incl_bits_/2);

  return 0;
}

void MinimizerIndex::AssignData_(const SequenceFile& seqs, bool index_reverse_strand) {
  /// The array seq_seeds_starts will hold starting positions for seeds comming from
  /// individual reads. I.e. each read will get a part of seq_seeds_ array which it will
  /// modify, and this array marks the part which is designated for this particular read.
  num_sequences_ = (index_reverse_strand == false) ? seqs.get_sequences().size() : seqs.get_sequences().size()*2;
  num_sequences_forward_ = seqs.get_sequences().size();
  reference_lengths_.resize(num_sequences_);
  reference_starting_pos_.resize(num_sequences_);
  headers_.resize(num_sequences_);

  // Calculate the total length of all reference sequences.
  int64_t total_seq_len = 0;
  for (int64_t i=0; i<num_sequences_forward_; i++) {
    total_seq_len += seqs.get_sequences()[i]->get_sequence_length();
  }
  if (index_reverse_strand) {
    total_seq_len *= 2;
  }

  // Preallocate space for data.
  data_serialized_.clear();
  data_serialized_.reserve(total_seq_len);

  for (int64_t i=0; i<num_sequences_forward_; i++) {
//    if (seqs.get_sequences()[i]->get_data_format() != kDataFormat2BitSparse) {
//      ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Data format is not 2-bit sparse!\n");
//    }
    reference_starting_pos_[i] = data_serialized_.size();
    auto rdata = seqs.get_sequences()[i]->get_data();
    auto rlen = seqs.get_sequences()[i]->get_data_length();
    data_serialized_.insert(data_serialized_.end(), rdata, rdata + rlen);
    headers_[i] = seqs.get_sequences()[i]->get_header();
    reference_lengths_[i] = rlen;
  }
  data_length_forward_ = data_serialized_.size();
  data_length_ = data_serialized_.size();

  if (index_reverse_strand == true) {
    for (int64_t i=0; i<num_sequences_forward_; i++) {
      int8_t *revcmp = seqs.get_sequences()[i]->GetReverseComplement();
      if (revcmp == NULL) {
        FATAL_REPORT(ERR_MEMORY, "Could not allocate memory for reverse complement!\n");
      }

      auto rlen = seqs.get_sequences()[i]->get_data_length();
      int64_t data_id = i + num_sequences_forward_;

      reference_starting_pos_[data_id] = data_serialized_.size();
      data_serialized_.insert(data_serialized_.end(), revcmp, revcmp + rlen);
      headers_[data_id] = seqs.get_sequences()[i]->get_header();
      reference_lengths_[data_id] = rlen;

      if (revcmp) {
        delete[] revcmp;
        revcmp = NULL;
      }
    }
    data_length_ = data_serialized_.size();
  }
}

void MinimizerIndex::AllocateSpaceForSeeds_(const SequenceFile& seqs, bool index_reverse_strand, int64_t num_shapes, int64_t max_seed_len, int64_t num_fwd_seqs, std::vector<int64_t>& seed_starts_for_seq, int64_t* total_num_seeds) {
  seeds_.clear();
  seed_starts_for_seq.clear();

  if (index_reverse_strand == true) {
    seed_starts_for_seq.resize(seqs.get_sequences().size() * 2);
  } else {
    seed_starts_for_seq.resize(seqs.get_sequences().size() * 1);
  }

  // Distribute bins for forward strand.
  *total_num_seeds = 0;
  for (int64_t i=0; i<seqs.get_sequences().size(); i++) {
    seed_starts_for_seq[i] = *total_num_seeds;
    *total_num_seeds += (seqs.get_sequences()[i]->get_sequence_length() - max_seed_len) * num_shapes;
  }

  // Distribute bins for reverse strand if required.
  if (index_reverse_strand == true) {
    for (int64_t i=0; i<seqs.get_sequences().size(); i++) {
      seed_starts_for_seq[seqs.get_sequences().size() + i] = *total_num_seeds;
      *total_num_seeds += (seqs.get_sequences()[i]->get_sequence_length() - max_seed_len) * num_shapes;
    }
  }

  /// The seed_list_will contain all seeds that were obtained from the input sequences.
  seeds_.resize(*total_num_seeds);
}

int MinimizerIndex::CollectAllSeedsForSeq_( const int8_t* seqdata, const int8_t* seqqual, int64_t seqlen,
                                            float min_avg_seed_qv, bool index_reverse_strand,
                                            uint64_t seq_id_fwd, uint64_t seq_id_rev,
                                            const std::vector<CompiledShape>& compiled_shapes,
                                            uint128_t* seed_list_fwd, uint128_t* seed_list_rev) {
  int64_t max_seed_len = 0;
  for (int32_t i=0; i<compiled_shapes.size(); i++) {
    max_seed_len = std::max(max_seed_len, (int64_t) compiled_shapes[i].shape().size());
  }

  /// The seqdata will be split in parts separated by N bases (similar to DALIGNER).
  /// Parts shorter than seed_len will be skipped.
  std::vector<int64_t> split_start;
  std::vector<int64_t> split_len;
  int64_t start = 0;
  for (int64_t i=0; i<seqlen; i++) {
    if (kBaseToBwa[seqdata[i]] > 3 || (i+1) >= seqlen) {
      int64_t current_len = ((i+1) >= seqlen) ? ((i - start) + 1) : (i - start);
      if (current_len >= max_seed_len) {
        split_start.push_back(start);
        split_len.push_back(current_len);
      }
      start = i + 1;
    }
  }

  int64_t seed_id = 0;

  double qv_sum = 0.0f;
  if (seqqual != NULL) {
    for (int32_t j=0; j<max_seed_len && j<seqlen; j++) {
      qv_sum += (seqqual[j] - 33.0);
    }
  }

  /// Iterate through all split regions.
  for (int64_t i=0; i<split_start.size(); i++) {
    /// We will use a buffer of 2-bit encoded bases, which will keep up to 32 bases at once.
    /// Fill the buffer for gapped spaced seed making. Buffer holds the next 8 bytes of data (max. 32 bases).
    /// Only seed_len bases are used, so the rest of the buffer is just to have the data ready.
    uint64_t buffer = 0;
    int64_t num_bases_ahead = sizeof(buffer) * 8 / 2;

    for (uint64_t pos=0; pos<(split_len[i] - max_seed_len); pos++) {
      /// Prepare the base-buffer. This would be an equivalent of a full-seed, including the
      /// don't care bases. This buffer will be used to extract the gapped spaced seed.
      /// Initialize the buffer.
      if (pos == 0) {
        for (int32_t j=0; j<num_bases_ahead && j<split_len[i]; j++) {
          int8_t seqbase = kBaseToBwa[seqdata[j + split_start[i]]];
          buffer |= (((uint64_t) seqbase) << (sizeof(buffer)*8 - j*2 - 2));
        }
        if (seqqual != NULL) {
          qv_sum = 0.0;
          for (int32_t j=0; j<max_seed_len && j<split_len[i]; j++) {
            qv_sum += (seqqual[j] - 33.0);
          }
        }
      } else {  /// Else, shift the buffer and re-fill.
        buffer = buffer << 2;
        if ((pos + num_bases_ahead) < split_len[i]) {
          int8_t seqbase = kBaseToBwa[seqdata[pos + split_start[i] + num_bases_ahead - 1]];
          buffer |= (((uint64_t) seqbase) << (0));
        }
        if (seqqual != NULL) {
          qv_sum = qv_sum - seqqual[pos-1] + seqqual[pos+max_seed_len-1];
        }
      }

      if (min_avg_seed_qv > 0 && (qv_sum / max_seed_len) < min_avg_seed_qv) {
        continue;
      }

      /// Extract gapped spaced seeds from the buffer.
      for (int64_t shape_id=0; shape_id<compiled_shapes.size(); shape_id++) {
        auto num_inclusive_bases = compiled_shapes[shape_id].num_incl_bits()/2;

        uint64_t seed = compiled_shapes[shape_id].CreateSeedFromShape(buffer);
        uint64_t key = SeedHashFunctionDefault_(seed, num_inclusive_bases);
        uint64_t position = pos + split_start[i]; // + 1;     // Make the position 1-based.
        seed_list_fwd[seed_id] = make_seed(key, seq_id_fwd, position);
//        seed_id += 1;

        if (index_reverse_strand) {
          uint64_t rev_seed = ReverseComplementSeed_(seed, num_inclusive_bases);
          uint64_t rev_key = SeedHashFunctionDefault_(rev_seed, num_inclusive_bases);
          uint64_t rev_position = (seqlen - position - 1); // | kIndexIdReverse64;
          seed_list_rev[seed_id] = make_seed(rev_key, seq_id_rev, rev_position);
        }

        seed_id += 1;
      }

    }
  }

  return seed_id;
}

inline uint64_t MinimizerIndex::SeedHashFunctionDefault_(uint64_t seed, int32_t k) {
  return seed;
}

inline uint64_t MinimizerIndex::ReverseComplementSeed_(uint64_t seed, int32_t num_bases) {
  uint64_t rev_seed = 0;
  uint64_t complement_seed = ((1 << (num_bases * 2)) - 1) - seed;      // Create a complement of the seed.
  for (int32_t i=0; i<num_bases; i++) {
    rev_seed |= (complement_seed & 0x03) << ((num_bases - i - 1) * 2);  // Reverse the bases.
    complement_seed >>= 2;
  }
  return rev_seed;
}

// Parameter window_len specifies the length of the window in the number of bases. For each base, there may be more than one seed
// (e.g. rev. complement, multiple indexes, etc.), and for this reason the parameter num_seeds_per_base is given.
// For example, if the array consists of a linear chain of only the forward strand and is obtained using only one hash function
// (gapped spaced index), then num_seeds_per_base = 1.
int MinimizerIndex::MakeMinimizers_(uint128_t* seed_list, int64_t num_seeds, int64_t num_seeds_per_base, int32_t window_len) {
  if (window_len > num_seeds) { return 1; }

  std::deque<int> q(num_seeds_per_base*window_len);

  // Setup the initial deque.
  for (int64_t i=0; i<window_len*num_seeds_per_base && i<num_seeds; i++) {
    // Remove smaller elements if any.
    while ( (!q.empty()) && GET_KEY_FROM_HIT(seed_list[i]) <= GET_KEY_FROM_HIT(seed_list[q.back()])) {
      q.pop_back();
    }
    q.push_back(i);
  }

  std::vector<int64_t> minimizer_indices;
  minimizer_indices.reserve(num_seeds);

  // If num_seeds_per_base is e.g. equal to 2, then
  // every other seed is potentially the reverse complement of the previous one.
  // Thus the sliding window will skip 2 instead of 1 seed.
  for (int64_t i=window_len*num_seeds_per_base; i<num_seeds; i+=num_seeds_per_base) {
    // Store the largest element of the previous window.
    minimizer_indices.push_back(q.front());
    // Remove the elements which are out of this window-
    while ((!q.empty()) && q.front() <= (i - window_len*num_seeds_per_base)) {
      q.pop_front();
    }
    // Remove smaller elements if any.
    while ( (!q.empty()) && GET_KEY_FROM_HIT(seed_list[i]) <= GET_KEY_FROM_HIT(seed_list[q.back()])) {
      q.pop_back();
    }

    q.push_back(i);
  }

  minimizer_indices.push_back(q.front());

  // Remove excess seeds.
  std::vector<bool> keep;
  keep.resize(num_seeds, false);
  for (int64_t i=0; i<minimizer_indices.size(); i++) {
    keep[minimizer_indices[i]] = true;
  }
  for (int64_t i=0; i<num_seeds; i++) {
    if (keep[i] == false) {
      set_invalid_seed(seed_list[i]);
    }
  }

  return 0;
}

int MinimizerIndex::FlagDuplicates_(uint128_t* seed_list, int64_t num_seeds) const {
  for (int64_t i=1; i<num_seeds; i++) {
    if (seed_list[i-1] == seed_list[i]) {
      set_invalid_seed(seed_list[i-1]);
//      seed_list[i-1] = 0;
    }
  }
  return 0;
}

void MinimizerIndex::DumpHash(std::string out_path, int32_t num_bases) {
  FILE *fp = fopen(out_path.c_str(), "w");
  fprintf (fp, "Key\tStart\tNum\tPositions...\n");
  if (fp != NULL) {
    for (auto it=hash_.begin(); it!=hash_.end(); it++) {
      uint64_t key = it->first;
      SeedHashValue shv = it->second;
  //    fprintf (fp, "%6X\t%ld\t%ld\n", key, shv.start, shv.num);
      std::string key_string = SeedToString(key, num_bases);
      fprintf (fp, "%s\t%ld\t%ld\t", key_string.c_str(), shv.start, shv.num);

      for (int64_t j=0; j<shv.num; j++) {
        fprintf (fp, "%ld ", GET_REAL_POS_FROM_HIT(seeds_[shv.start + j]));
      }
      fprintf (fp, "\n");
    }
    fclose(fp);
  }
}

void MinimizerIndex::DumpHashSortedByCount(std::string out_path, int32_t num_bases) {
  std::vector<std::tuple<std::string, int64_t, int64_t> > sorted_hash;
  sorted_hash.reserve(hash_.size());
  for (auto it=hash_.begin(); it!=hash_.end(); it++) {
    uint64_t key = it->first;
    SeedHashValue shv = it->second;
    std::string key_string = SeedToString(key, num_bases);
    sorted_hash.push_back(std::make_tuple(key_string, shv.start, shv.num));
  }
  std::sort(sorted_hash.begin(), sorted_hash.end(),
            [](const std::tuple<std::string, int64_t, int64_t>& a, const std::tuple<std::string, int64_t, int64_t>& b) { return ((std::get<2>(a)) > (std::get<2>(b))); });

  FILE *fp = fopen(out_path.c_str(), "w");
  fprintf (fp, "Key\tStart\tNum\n");
  for (int64_t i=0; i<sorted_hash.size(); i++) {
    fprintf (fp, "%s\t%ld\t%ld\n", std::get<0>(sorted_hash[i]).c_str(), std::get<1>(sorted_hash[i]), std::get<2>(sorted_hash[i]));
  }
  fclose(fp);
}

void MinimizerIndex::DumpHashSortedByName(std::string out_path, int32_t num_bases) {
  std::vector<std::tuple<std::string, int64_t, int64_t> > sorted_hash;
  sorted_hash.reserve(hash_.size());
  for (auto it=hash_.begin(); it!=hash_.end(); it++) {
    uint64_t key = it->first;
    SeedHashValue shv = it->second;
    std::string key_string = SeedToString(key, num_bases);
    sorted_hash.push_back(std::make_tuple(key_string, shv.start, shv.num));
  }
  std::sort(sorted_hash.begin(), sorted_hash.end(), [](const std::tuple<std::string, int64_t, int64_t>& a, const std::tuple<std::string, int64_t, int64_t>& b) { return ((std::get<0>(a)) < (std::get<0>(b))); });

  FILE *fp = fopen(out_path.c_str(), "w");
  fprintf (fp, "Key\tStart\tNum\n");
  for (int64_t i=0; i<sorted_hash.size(); i++) {
    fprintf (fp, "%s\t%ld\t%ld\n", std::get<0>(sorted_hash[i]).c_str(), std::get<1>(sorted_hash[i]), std::get<2>(sorted_hash[i]));
  }
  fclose(fp);
}

void MinimizerIndex::DumpSeeds(std::string out_path, int32_t num_bases) {
  FILE *fp = fopen(out_path.c_str(), "w");
  if (fp != NULL) {
    fprintf (fp, "Key\tPosition\tSeqID\tis_rev\tkey\n");
    for (int64_t i=0; i<seeds_.size(); i++) {
      uint64_t key = (uint64_t) (seeds_[i] >> 64);
      uint64_t coded_pos = seeds_[i] & 0x0000000000000000FFFFFFFFFFFFFFFF;
      uint64_t ref_id = GET_SEQ_ID_FROM_HIT(seeds_[i]);
      IndexPos ipos(coded_pos);
  //    fprintf (fp, "%6X %ld\n", key, ipos.get_pos());
      std::string key_string = SeedToString(key, num_bases);
      fprintf (fp, "%s %ld %ld %s %15ld\n", key_string.c_str(), ipos.get_pos(), ipos.seq_id, ((ipos.is_rev() == true) ? "r" : "f"), key);
    }
    fclose(fp);
  }
}

int MinimizerIndex::OccurrenceStatistics_(double percentil, int32_t num_threads, double* ret_avg, double* ret_stddev, double *ret_percentil_val) const {
  if (percentil < 0.0 || percentil > 1.0) { return 1; }

  std::vector<int32_t> key_counts;
  key_counts.resize(hash_.size(), 0);
  int64_t currkey = 0;
  double avg = 0.0, stddev = 0.0, sum = 0.0;
  // Initialize the array of counts. Needed for percentil calculation.
  // Also, calculate the avg on the fly.
  for (auto it = hash_.begin(); it != hash_.end(); it++) {
    key_counts[currkey++] = it->second.num;
    avg += it->second.num;
    sum += it->second.num;
  }
  if (key_counts.size() > 0) { avg /= key_counts.size(); }

  // Calculate the standard deviation.
  for (int64_t i=0; i<key_counts.size(); i++) {
    stddev += (key_counts[i] - avg) * (key_counts[i] - avg);
  }
  if (key_counts.size() > 1) { stddev /= (key_counts.size() - 1); } // Unbiased estimator.
  stddev = sqrt(stddev);

  // Calculate the percentil.
  pquickSort(&(key_counts[0]), key_counts.size(), num_threads);
  double perc_val = key_counts[percentil * (key_counts.size() - 1)];
//  LOG_ALL("key_counts.size() = %ld\n", key_counts.size());
//  LOG_ALL("Simple percentil (percentage of the total array): %f\n", perc_val);

//  double perc_sum = 0.0;
//  for (int64_t i=0; i<key_counts.size(); i++) {
//    if ((perc_sum + key_counts[i]) >= (sum * percentil)) {
//      perc_val = key_counts[i];
//      break;
//    }
//    perc_sum += key_counts[i];
//  }
//  LOG_ALL("Weighted percentil (percentage of the total array): %f\n", perc_val);

  *ret_avg = avg;
  *ret_stddev = stddev;
  *ret_percentil_val = perc_val;

  return 0;
}

int MinimizerIndex::MakeSeedListDense_(uint128_t* seed_list, int64_t num_seeds) {
  int64_t offset = 0;
  int64_t i = 0;
  for (i=0; (i+offset)<num_seeds; i++) {
    while ((i+offset) < num_seeds && invalid_seed(seed_list[i+offset]) == true) {
      offset += 1;
    }
    if ((i+offset) < num_seeds) {
      seed_list[i] = seed_list[i+offset];
    }
  }
  return i;
}

void MinimizerIndex::CollectSeeds_(const int8_t* seqdata, const int8_t* seqqual, int64_t seqlen,
                                   bool index_reverse_strand, uint64_t seq_id_fwd, uint64_t seq_id_rev,
                                   int64_t max_seed_len, float min_avg_seed_qv, bool use_minimizers, int32_t minimizer_window_len,
                                   const std::vector<CompiledShape> &index_shapes, std::vector<uint128_t> &seed_list) const {

  int64_t num_seeds_fwd = (seqlen - max_seed_len + 1) * index_shapes.size();

  seed_list.clear();
  seed_list.resize(num_seeds_fwd * 2);     // Allocate maximum space for the seeds. Factor 2 is for the reverse complement.

  int64_t num_seeds_processed = CollectAllSeedsForSeq_(seqdata, seqqual, seqlen, min_avg_seed_qv,
                                                       index_reverse_strand, seq_id_fwd, seq_id_rev, index_shapes, &(seed_list[0]), &(seed_list[num_seeds_fwd]));

  if (use_minimizers == true) {
  //  MakeMinimizers_(&(seed_list[0]), num_seeds_processed, 2*compiled_shapes.size(), minimizer_window_len);     // 2 only refers to the fwd and rev complement (2 seeds per base).
    MakeMinimizers_(&(seed_list[0]), num_seeds_fwd, index_shapes.size(), minimizer_window_len);                  // Forward strand.
    if (index_reverse_strand) {
      MakeMinimizers_(&(seed_list[num_seeds_fwd]), num_seeds_fwd, index_shapes.size(), minimizer_window_len);      // Reverse complement.
    }

    int64_t num_dense_seeds = MakeSeedListDense_(&(seed_list[0]), seed_list.size());
    seed_list.resize(num_dense_seeds);
  }
}

int64_t MinimizerIndex::RawPositionConverterWithRefId(int64_t raw_position, int64_t reference_index, int64_t query_length, int64_t *ret_absolute_position, int64_t *ret_relative_position, SeqOrientation *ret_orientation, int64_t *ret_reference_index_with_reverse) const {
  if (raw_position < 0 || raw_position >= data_length_)
    return -2;

  if (reference_index < 0) {
    LogSystem::GetInstance().Error(SEVERITY_INT_ERROR, __FUNCTION__, LogSystem::GetInstance().GenerateErrorMessage(ERR_UNEXPECTED_VALUE, "Offending variable: reference_index. Values: reference_index = %ld, raw_position = %ld, data_length = %ld.", reference_index, raw_position, data_length_));
    return reference_index;
  }

  int64_t relative_pos = (raw_position - reference_starting_pos_[(uint64_t) reference_index]);
  SeqOrientation orientation = kForward;

  if (ret_reference_index_with_reverse != NULL)
    *ret_reference_index_with_reverse = reference_index;

  if (((uint64_t) reference_index) >= num_sequences_forward_) {
    // Relative position has to be changed, because, from the outside, it is expected that we have reverse-complemented the seed and not the reference sequence.
//    relative_pos = reference_lengths_[(uint64_t) reference_index] - relative_pos - query_length - 1 - reference_index;      // The '-1' is to compensate for the '!' character added at the end of every sequence in the data array.
    relative_pos = reference_lengths_[(uint64_t) reference_index] - relative_pos - query_length - 1;

    // Unlike BWA, we haven't reversed the order of sequences when their reverse complements were added to the index. That is why we only need to subtract the number of forward sequences, and not do (2*num_forward_sequences - gene_idx - 1).
    reference_index = reference_index - ((int64_t) num_sequences_forward_);
    orientation = kReverse;
  }

  if (ret_absolute_position != NULL)
    *ret_absolute_position = raw_position;

  if (ret_relative_position != NULL)
    *ret_relative_position = relative_pos;

  if (ret_orientation != NULL)
    *ret_orientation = orientation;

  return reference_index;
}

void MinimizerIndex::CollectSeeds(const int8_t* seqdata, const int8_t* seqqual, int64_t seqlen,
                                  float min_avg_seed_qv, bool index_reverse_strand,
                                  std::vector<uint128_t>& seed_list) const {
  CollectSeeds_(seqdata, seqqual, seqlen, index_reverse_strand, 0, 0, max_seed_len_, min_avg_seed_qv, use_minimizers_, minimizer_window_len_, index_shapes_, seed_list);
}

std::string VerboseSeed(uint128_t seed) {
  std::stringstream ss;
  ss << "Key = " << MinimizerIndex::seed_key(seed) << " seq_id = " << MinimizerIndex::seed_seq_id(seed) << " pos = " << MinimizerIndex::seed_position(seed);
  return ss.str();
}

int MinimizerIndex::KeyLookup(uint64_t seed_key, const uint128_t** seeds,
                                  int64_t *num_seeds) const {
  *seeds = NULL;
  *num_seeds = 0;

  auto it = hash_.find(seed_key);
  if (it == hash_.end()) { return 1; }

  *seeds = &(seeds_[it->second.start]);
  *num_seeds = it->second.num;

  return 0;
}

int MinimizerIndex::Find(const int8_t* seed, int32_t seed_len,
                         std::vector<const uint128_t*>& hits,
                         std::vector<int64_t> &num_hits) const {

  std::vector<uint64_t> keys;
  CalcKeysFromSeed(seed, seed_len, keys);

  hits.clear();
  num_hits.clear();

  for (auto& key: keys) {
    const uint128_t *seeds = NULL;
    int64_t num_seeds = 0;

    int lookup_ret = KeyLookup(key, &seeds, &num_seeds);
    if (!lookup_ret) {
      hits.push_back(seeds);
      num_hits.push_back(num_seeds);
    }
  }

  if (hits.size() == 0) {
    return 1;
  }

  return 0;
}

int MinimizerIndex::FindAndJoin(const int8_t* seed, int32_t seed_len, std::vector<uint128_t> &hits) const {
  std::vector<const uint128_t*> hits_ptr;
  std::vector<int64_t> num_hits_ptr;
  Find(seed, seed_len, hits_ptr, num_hits_ptr);

  hits.clear();
  int64_t total_hits = 0;
  for (int32_t i=0; i<hits_ptr.size(); i++) {
    total_hits += num_hits_ptr[i];
  }
  hits.reserve(total_hits);

//  int64_t c = 0;
  for (int32_t i=0; i<hits_ptr.size(); i++) {
    hits.insert(hits.end(), hits_ptr[i], hits_ptr[i] + num_hits_ptr[i]);
//    memcpy(&(hits[c]), hits_ptr[i], num_hits_ptr[i]);
//    c += num_hits_ptr[i];
  }

  if (hits.size() == 0) {
    return 1;
  }

  return 0;
}

void MinimizerIndex::CalcKeysFromSeed(const int8_t* seed, int32_t seed_len,
                                       std::vector<uint64_t> &keys) const {
  uint64_t buffer = 0;
  for (int32_t j=0; j<seed_len && j<32; j++) {
    int8_t seqbase = kBaseToBwa[seed[j]];
    buffer |= (((uint64_t) seqbase) << (sizeof(buffer)*8 - j*2 - 2)); // -2 is for a 62bit lshift.
  }

  keys.clear();
  for (int32_t i=0; i<lookup_shapes_.size(); i++) {
    uint64_t seed = lookup_shapes_[i].CreateSeedFromShape(buffer);
    keys.push_back(seed);
  }
}

} /* Namespace is. */
