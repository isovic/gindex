/*
 * experimental.cc
 *
 *  Created on: Feb 26, 2017
 *      Author: isovic
 */

#include "minimizer_index.h"

#include "log_system/log_system.h"
#include "omp_sort.hpp"
#include <omp.h>
#include <algorithm>
#include <tuple>
#include <deque>
#include "utility/utility_general.h"
#include "utility/tictoc.h"
#include "minimizer_generator.h"

namespace is {

int MinimizerIndex::Create(const SequenceFile &seqs,
           float min_avg_seed_qv, bool index_reverse_strand,
           bool use_minimizers, int32_t minimizer_window_len,
           int32_t num_threads, bool verbose) {
  Clear_();

  use_minimizers_ = use_minimizers;
  minimizer_window_len_ = minimizer_window_len;
  int64_t n_seqs = seqs.get_sequences().size();

  TicToc tt_all;
  tt_all.start();

  AssignData_(seqs, index_reverse_strand);

  // The seed_list_will contain all seeds that were obtained from the input sequences.
  TicToc tt_alloc;
  tt_alloc.start();

  seeds_.clear();
  int64_t approx_num_seeds = (data_length_ / (std::max(1, minimizer_window_len - 1))) + 1;     // An estimate. Uses minimizer_window_len - 1 to allocate more space.
  seeds_.reserve(approx_num_seeds);

  tt_alloc.stop();

  if (verbose) {
    LOG_ALL("Allocated memory for a list of %ld seeds (128 bits each) (%.5f sec, diff: %.5f sec).\n", seeds_.capacity(), tt_alloc.get_secs(), tt_all.get_secs_current());
    LOG_ALL("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str());
    LOG_ALL("Collecting seeds.\n");
    if (use_minimizers) {
      LOG_ALL("Minimizer seeds will be used. Minimizer window is %ld.\n", minimizer_window_len);
    }
  }

  TicToc tt_minimizers;
  tt_minimizers.start();

  for (int64_t i=0; i<num_sequences_; i++) {
    if (verbose) {
      LOG_ALL("\r%s Sequence: %ld/%ld, len: %ld, name: '%s'", FormatMemoryConsumptionAsString().c_str(), (i + 1), num_sequences_, reference_lengths_[i], TrimToFirstSpace(headers_[i]).c_str());
    }

    int8_t *seqdata = (int8_t *) (&data_serialized_[0] + reference_starting_pos_[i]);
    int64_t seqlen = reference_lengths_[i];
    uint64_t seq_id = i;

    // Collect all seeds.
    int64_t num_seeds_processed = AddSeeds_(seqdata, seqlen, seq_id,
                                            use_minimizers, minimizer_window_len, index_shapes_, seeds_);
  }

  if (verbose) {
    LOG_NEWLINE;
    LOG_ALL("Final memory allocation after collecting seeds: %s\n", FormatMemoryConsumptionAsString().c_str());
  }

//  if (verbose) { LOG_NEWLINE; }

//  DumpSeeds("temp/seeds.dense.minimizers.csv", max_incl_bits_/2);

  if (verbose) {
    LOG_ALL("Sorting the seeds using %ld threads.\n", num_threads);
  }

  TicToc tt_sort;
  tt_sort.start();

  pquickSort(&(seeds_[0]), seeds_.size(), num_threads);

  tt_sort.stop();
  //    FlagDuplicates_(&(seeds_[seed_starts_for_seq[i]]), num_seeds_processed);
//  LogSystem::GetInstance().Log(VERBOSE_LEVEL_ALL, true, FormatString("Making unique (%.5f sec, diff: %.5f sec).\n", tt_sort.get_secs(), (((float) (clock() - diff_time))/CLOCKS_PER_SEC)), std::string(__FUNCTION__));
//  diff_time = clock();

  if (verbose) {
//    LOG_ALL("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str());
    LOG_ALL("Generating the hash table.\n");
  }

  TicToc tt_hash;
  tt_hash.start();
  ConstructHash_();
  tt_hash.stop();

  // Calculate a cutoff threshold as a percentil of the occurrence of a seed.
  if (verbose) {
//    LOG_ALL("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str());
    LOG_ALL("Calculating the distribution statistics for key counts.\n");
  }

  OccurrenceStatistics_(percentil_, num_threads, &avg_seed_occurrence_, &max_seed_occurrence_, &stddev_seed_occurrence_, &count_cutoff_);

  if (verbose) {
//    LOG_ALL("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str());
    LOG_ALL("Index statistics: average key count = %f, max key count = %f, std dev = %f, percentil (%.2f%%) (count cutoff) = %f\n",
            avg_seed_occurrence_, max_seed_occurrence_, stddev_seed_occurrence_, percentil_*100.0, count_cutoff_);
  }

  ConstructHashUnderCutoff();

  if (verbose) {
    LOG_ALL("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str());
  }

//  DumpSeeds("temp/seeds.csv", max_incl_bits_/2);
//  DumpHash("temp/hash.csv", max_incl_bits_/2);
//  DumpHashSortedByCount("temp/hash.minimizers.sorted.csv", max_incl_bits_/2);

  return 0;
}

int64_t MinimizerIndex::AddSeeds_(const int8_t *seqdata, int64_t seqlen, int64_t seq_id, bool use_minimizers, int32_t minimizer_window_len,
                                  const std::vector<CompiledShape> &compiled_shapes, std::vector<uint128_t>& seeds) const {
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

  /// Iterate through all split regions.
  for (int64_t i=0; i<split_start.size(); i++) {
    // Start generating minimzers from scratch.
    auto mg = createMinimizerGenerator(minimizer_window_len);

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
      } else {  /// Else, shift the buffer and re-fill.
        buffer = buffer << 2;
        if ((pos + num_bases_ahead) < split_len[i]) {
          int8_t seqbase = kBaseToBwa[seqdata[pos + split_start[i] + num_bases_ahead - 1]];
          buffer |= (((uint64_t) seqbase) << (0));
        }
      }

      /// Extract gapped spaced seeds from the buffer.
      for (int64_t shape_id=0; shape_id<compiled_shapes.size(); shape_id++) {
        auto num_inclusive_bases = compiled_shapes[shape_id].num_incl_bits()/2;

        uint64_t seed = compiled_shapes[shape_id].CreateSeedFromShape(buffer);
        uint64_t key = MinimizerIndex::SeedHashFunctionDefault_(seed, num_inclusive_bases);
        uint64_t position = pos + split_start[i]; // + 1;     // Make the position 1-based.

        uint128_t seed_in = Seed::make_seed(key, seq_id, position);

        uint128_t seed_out = 0;
        int ret_val_mg = mg->yield(seed_in, seed_out);

        if (!ret_val_mg) {
          seeds.emplace_back(seed_out);
        }

        seed_id += 1;
      }

    }
  }

  return seed_id;
}

void ReverseComplement(const int8_t *seqdata, int64_t seqlen, std::vector<int8_t>& ret) {
  ret.resize(seqlen);
  for (int64_t i=0; i<seqlen; i++) {
    ret[i] = kBaseComplement[seqdata[seqlen-i-1]];
  }
}

void MinimizerIndex::CollectIndexSeeds(const int8_t *seqdata, const int8_t *seqqual, int64_t seqlen,
                  float min_avg_seed_qv, bool index_reverse_strand,
                  bool use_minimizers, int32_t minimizer_window_len,
                  std::vector<uint128_t> &seed_list) const {
//  CollectSeeds_(seqdata, seqqual, seqlen, index_reverse_strand, 0, 0, max_seed_len_, min_avg_seed_qv, use_minimizers, minimizer_window_len, index_shapes_, seed_list);
  seed_list.clear();
  AddSeeds_(seqdata, seqlen, 0, use_minimizers, minimizer_window_len, index_shapes_, seed_list);

  if (index_reverse_strand) {
    std::vector<int8_t> rev;
    ReverseComplement(seqdata, seqlen, rev);
    AddSeeds_(&rev[0], rev.size(), 0, use_minimizers, minimizer_window_len, index_shapes_, seed_list);
  }
}

void MinimizerIndex::CollectLookupSeeds(const int8_t *seqdata, const int8_t *seqqual, int64_t seqlen,
                  float min_avg_seed_qv, bool index_reverse_strand,
                  bool use_minimizers, int32_t minimizer_window_len,
                  std::vector<uint128_t> &seed_list) const {
//  CollectSeeds_(seqdata, seqqual, seqlen, index_reverse_strand, 0, 0, max_seed_len_, min_avg_seed_qv, use_minimizers, minimizer_window_len, lookup_shapes_, seed_list);
  seed_list.clear();
  AddSeeds_(seqdata, seqlen, 0, use_minimizers, minimizer_window_len, lookup_shapes_, seed_list);

  if (index_reverse_strand) {
    std::vector<int8_t> rev;
    ReverseComplement(seqdata, seqlen, rev);
    AddSeeds_(&rev[0], rev.size(), 0, use_minimizers, minimizer_window_len, lookup_shapes_, seed_list);
  }
}

}
