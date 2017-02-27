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

namespace is {

int MinimizerIndex::Create(const SequenceFile &seqs,
           float min_avg_seed_qv, bool index_reverse_strand,
           bool use_minimizers, int32_t minimizer_window_len,
           int32_t num_threads, bool verbose=false) {
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
  int64_t approx_num_seeds = data_length_ / minimizer_window_len + 1;     // An estimate.
  seeds_.reserve(approx_num_seeds);
  tt_alloc.stop();

  if (verbose) {
    LOG_ALL("Allocated memory for a list of %ld seeds (128 bits each) (%.5f sec, diff: %.5f sec).\n", seeds_.capacity(), tt_alloc.get_secs(), tt_all.get_secs_current());
    LOG_ALL("Memory consumption: %s\n", FormatMemoryConsumptionAsString().c_str());
    LOG_ALL("Collecting seeds with %ld threads.\n", num_threads);
    if (use_minimizers) {
      LOG_ALL("Minimizer seeds will be used. Minimizer window is %ld.\n", minimizer_window_len);
    }
  }

  TicToc tt_minimizers;
  tt_minimizers.start();

  for (int64_t i=0; i<seqs.get_sequences().size(); i++) {
    uint32_t thread_id = omp_get_thread_num();
    if (thread_id == 0 && verbose) {
      LOG_ALL("\rSequence %ld/%ld", (i + 1), seqs.get_sequences().size());
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
    if (verbose) { LOG_ALL("Removing excess seeds.\n"); }
    int64_t num_dense_seeds = MakeSeedListDense_(&(seeds_[0]), seeds_.size());
    seeds_.resize(num_dense_seeds);
  }

//  DumpSeeds("temp/seeds.dense.minimizers.csv", max_incl_bits_/2);

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

  ConstructHash_();

  // Calculate a cutoff threshold as a percentil of the occurrence of a seed.
  if (verbose) {
    LOG_ALL("Calculating the distribution statistics for key counts (%.5f sec, diff: %.5f sec).\n", (((float) (clock() - absolute_time))/CLOCKS_PER_SEC), (((float) (clock() - diff_time))/CLOCKS_PER_SEC));
  }

  OccurrenceStatistics_(percentil_, num_threads, &avg_seed_occurrence_, &max_seed_occurrence_, &stddev_seed_occurrence_, &count_cutoff_);

  if (verbose) {
    LOG_ALL("Index statistics: average key count = %f, max key count = %f, std dev = %f, percentil (%.2f%%) (count cutoff) = %f\n", avg_seed_occurrence_, max_seed_occurrence_, stddev_seed_occurrence_, percentil_*100.0, count_cutoff_);
  }

  ConstructHashUnderCutoff();

//  DumpSeeds("temp/seeds.csv", max_incl_bits_/2);
//  DumpHash("temp/hash.csv", max_incl_bits_/2);
//  DumpHashSortedByCount("temp/hash.minimizers.sorted.csv", max_incl_bits_/2);

  return 0;

}

void MinimizerIndex::CollectIndexSeeds(const int8_t *seqdata, const int8_t *seqqual, int64_t seqlen,
                  float min_avg_seed_qv, bool index_reverse_strand,
                  bool use_minimizers, int32_t minimizer_window_len,
                  std::vector<uint128_t> &seed_list) const {
  CollectSeeds_(seqdata, seqqual, seqlen, index_reverse_strand, 0, 0, max_seed_len_, min_avg_seed_qv, use_minimizers, minimizer_window_len, index_shapes_, seed_list);
}

void MinimizerIndex::CollectLookupSeeds(const int8_t *seqdata, const int8_t *seqqual, int64_t seqlen,
                  float min_avg_seed_qv, bool index_reverse_strand,
                  bool use_minimizers, int32_t minimizer_window_len,
                  std::vector<uint128_t> &seed_list) const {
  CollectSeeds_(seqdata, seqqual, seqlen, index_reverse_strand, 0, 0, max_seed_len_, min_avg_seed_qv, use_minimizers, minimizer_window_len, lookup_shapes_, seed_list);

}

}
