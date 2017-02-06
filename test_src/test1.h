#ifndef TEST1_
#define TEST1_

#include "../src/index/minimizer_index.h"
#include "gtest/gtest.h"

TEST(GenerateIndex, Normal) {

//  std::string ref = "test-data/lambda/NC_001416.fa";
  std::string ref = "test-data/simple/test1.fa";

  SequenceFile ref_seqs(ref);
  ref_seqs.ConvertDataFormat(kDataFormat2BitSparse);

  std::vector<std::string> index_shapes = {"1111110111111"};

  auto index = is::createMinimizerIndex();
//  index->Create(ref_seqs, index_shapes, 0.0f, true, false, 0, 4);
  index->Create(ref_seqs, index_shapes, 0.0f, true, true, 5, 4);

  index->DumpHash("dump-hash.txt", 12);
  index->DumpHashSortedByName("dump-sorted-hash.txt", 12);
  index->DumpSeeds("dump-seeds.txt", 12);

  // std::string cigar_string = is::CigarToString(al->getResults()[0].cigar);
  // std::cout << cigar_string << "\n";

  // EXPECT_EQ(cigar_string, std::string("8964="));
}

#endif
