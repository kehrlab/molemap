#ifndef PARSER_H
#define PARSER_H

# include <seqan/arg_parse.h>

struct bcmapOptions{
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

  bcmapOptions() :
  kmer_index_name("Index"), read_index_name("ReadIndex"), k(31), mini_window_size(61), max_window_size(300000), max_gap_size(20000),output_file("barcode_windows.bed"),l(10000) , s(0), threads(16), Sort(false), CoverageAnalysis(false)
  {}
};

seqan::ArgumentParser::ParseResult parseCommandLine(bcmapOptions & options, int argc, char const ** argv);
void printParseResults(bcmapOptions options);

#endif
