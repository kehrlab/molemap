#ifndef PARSER_H
#define PARSER_H

# include <seqan/arg_parse.h>
# include <seqan/misc/interval_tree.h>


struct indexOptions{
  std::string reference_file;
  std::string kmer_index_name;
  std::string preset;
  uint32_t k;
  uint32_t m;//int64_t bucket_count;

  indexOptions():
  k(31), m(61)// bucket_count(3221225472)
  {}
};

struct mapOptions{
  std::string readfile1;
  std::string readfile2;
  std::string kmer_index_name;
  std::string read_index_name;
  uint32_t k;
  uint32_t mini_window_size;
  uint32_t max_window_size;
  uint32_t max_gap_size;
  uint32_t max_abundance;
  std::string output_file;
  unsigned s;
  unsigned l;
  unsigned threads;
  bool Sort;
  bool CoverageAnalysis;
  uint32_t batchSize;

  mapOptions() :
  kmer_index_name("Index"), read_index_name("ReadIndex"), k(31), mini_window_size(61), max_window_size(300000), max_gap_size(20000), max_abundance(20),output_file("barcode_windows.bed"),l(10000) , s(0), threads(16), Sort(false), CoverageAnalysis(false), batchSize(5000)
  {}
};

struct region{
  std::string chrom;
  uint32_t start;
  uint32_t end;
};

struct longmapOptions{
  std::string readfile1;
  std::string kmer_index_name;
  uint32_t k;
  uint32_t mini_window_size;
  uint32_t max_gap_size;
  uint32_t max_abundance;
  std::string output_file;
  std::string output_format;
  unsigned s;
  unsigned l;
  unsigned threads;
  std::string readGroup;
  std::string readGroupId;
  bool regionDefined;
  // std::vector<region> regions;
  std::map<std::string,seqan::IntervalTree<uint32_t, bool>> regions;
  uint32_t batchSize;

  longmapOptions() :
  kmer_index_name("Index"), k(31), mini_window_size(61), max_gap_size(20000), max_abundance(20), output_file("stdout"), output_format("sam"),l(500) , s(0), threads(16), readGroup(""), readGroupId(""), regionDefined(false), batchSize(5000)
  {}
};

struct getOptions{
  std::string readfile1;
  std::string readfile2;
  std::string read_index_name;
  std::string barcodes;
  std::string output_file;

  getOptions() :
  output_file("molemapGetOut"), read_index_name("BarcodeIndex")
  {}
};



seqan::ArgumentParser::ParseResult parseCommandLine_index(indexOptions & options, int argc, char const ** argv);
void printParseResults_index(indexOptions & options);
seqan::ArgumentParser::ParseResult parseCommandLine_map(mapOptions & options, int argc, char const ** argv);
void printParseResults_map(mapOptions & options);
seqan::ArgumentParser::ParseResult parseCommandLine_long_map(longmapOptions & options, int argc, char const ** argv);
void printParseResults_long_map(longmapOptions & options);
seqan::ArgumentParser::ParseResult parseCommandLine_get(getOptions & options, int argc, char const ** argv);
void printParseResults_get(getOptions & options);
void loadIndexParameters(uint32_t & k, uint32_t & m, std::string & IndexName);

#endif
