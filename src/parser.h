#ifndef PARSER_H
#define PARSER_H

# include <seqan/arg_parse.h>

struct indexOptions{
  std::string reference_file;
  std::string kmer_index_name;
  unsigned k;
  unsigned m;//int64_t bucket_count;

  indexOptions():
  k(31), m(61)// bucket_count(3221225472)
  {}
};

struct mapOptions{
  std::string readfile1;
  std::string readfile2;
  std::string kmer_index_name;
  std::string read_index_name;
  unsigned k;
  unsigned mini_window_size;
  unsigned max_window_size;
  unsigned max_gap_size;
  std::string output_file;
  unsigned s;
  unsigned l;
  unsigned threads;
  bool Sort;
  bool CoverageAnalysis;

  mapOptions() :
  kmer_index_name("Index"), read_index_name("ReadIndex"), k(31), mini_window_size(61), max_window_size(300000), max_gap_size(20000),output_file("barcode_windows.bed"),l(10000) , s(0), threads(16), Sort(false), CoverageAnalysis(false)
  {}
};

struct getOptions{
  std::string readfile1;
  std::string readfile2;
  std::string read_index_name;
  std::string barcodes;
  std::string output_file;

  getOptions() :
  output_file("bcmapGetOut"), read_index_name("BarcodeIndex")
  {}
};

seqan::ArgumentParser::ParseResult parseCommandLine_index(indexOptions & options, int argc, char const ** argv);
void printParseResults_index(indexOptions & options);
seqan::ArgumentParser::ParseResult parseCommandLine_map(mapOptions & options, int argc, char const ** argv);
void printParseResults_map(mapOptions & options);
seqan::ArgumentParser::ParseResult parseCommandLine_get(getOptions & options, int argc, char const ** argv);
void printParseResults_get(getOptions & options);


#endif
