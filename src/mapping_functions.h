#ifndef MAPPING_FUNCTIONS_H
#define MAPPING_FUNCTIONS_H

# include "parser.h"
# include "index.h"
# include <zlib.h>
# include <htslib/kseq.h>

#ifndef KSEQ_GZ
#define KSEQ_GZ
KSEQ_INIT(gzFile, gzread)
#endif  // KSEQ_GZ

struct barcodeData{
  std::string barcode;
  std::vector<Dna5String> Reads;

  barcodeData(){
    Reads.reserve(1000);
  }

  void append(kseq_t * seq1, kseq_t * seq2){
    Reads.push_back(seq1->seq.s);
    Reads.push_back(seq2->seq.s);
  }
};

struct ReadData{
  Dna5String read;
  CharString id;
  CharString qual;
};

int mapLinkedZipped(openAddressingKmerHashtable & Index, mapOptions & options, int argc, char const ** argv);
void MapKmerListLinked(std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, uint8_t & window_count, std::string barcode, std::vector<result_t> & results, std::vector<std::string> & lookChrom, std::vector<uint32_t> & histogram, mapOptions & options);
void trimmWindowLinked(std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrstart, std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrk, std::tuple<double,uint8_t,uint32_t,uint32_t> & candidate, std::vector<float> & lookQual);
void AppendPos(std::vector<std::tuple <uint8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, const int64_t & hash, openAddressingKmerHashtable & Index, uint8_t & minimizer_active_bases, mapOptions & options);
uint8_t getBarcodeLength(std::string & readfile1, std::streampos & readfile1_size);

// ################################################################################

int mapLongUnzipped(openAddressingKmerHashtable & Index, longmapOptions & options, int argc, char const ** argv);
void readBatch(kseq_t * seq1, std::vector<ReadData> & Batch, int & batchSize);
int mapLongZipped(openAddressingKmerHashtable & Index, longmapOptions & options, int argc, char const ** argv);
void moveFileToStart(SeqFileIn & file1, std::streampos & startpos, int & t);
int MapKmerListLong(std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>> & kmer_list, BamAlignmentRecord & result, std::vector<uint32_t> & histogram, uint32_t readLength, longmapOptions & options);
void trimmWindowLong(std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>> & kmer_list, std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrstart, std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrk, std::tuple<double,uint8_t,uint32_t,uint32_t,bool> & candidate, std::vector<float> & lookQual);
void AppendPosLong(std::vector<std::tuple <uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>> & kmer_list, const int64_t & hash, openAddressingKmerHashtable & Index, uint8_t & minimizer_active_bases, uint32_t order, longmapOptions & options);

#endif
