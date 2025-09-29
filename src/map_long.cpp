# include <seqan/seq_io.h>
# include <seqan/bam_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <iostream>
# include <fstream>
# include <time.h>
# include "map_long.h"
# include "parser.h"
# include "functions.h"
# include "IO_functions.h"
# include "index.h"
# include "mapping_functions.h"
using namespace seqan;

void AppendPosLong(std::vector<std::tuple <uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>> & kmer_list, const int64_t & hash, openAddressingKmerHashtable & Index,uint8_t & minimizer_active_bases, uint32_t order);
void moveFileToStart(SeqFileIn & file1, std::streampos & startpos, int & t);

int maplong(int argc, char const ** argv){

  // parsing command line arguments
  longmapOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine_long_map(options, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK) return res;
  printParseResults_long_map(options);

  if(options.k > 21){
    std::cerr << "WARNING: Large value for k = " << options.k << ". For reads with low base calling accuracy choose an index with smaller k for optimal results.\n";
  }

  std::cerr << "Reading in the k-mer index";

  openAddressingKmerHashtable Index;
  readKmerIndex(Index, options.kmer_index_name);

  std::cerr <<"..done.\n";

  std::string filenameExtension = options.readfile1.substr(options.readfile1.rfind('.')+1);

  if(filenameExtension == "fq" || filenameExtension == "fastq"){
    if(checkReadfile(options.readfile1)){return 1;};
    if(mapLongUnzipped(Index, options, argc, argv)){
      return 1;
    }
  }else if(filenameExtension == "gz"){
    if(mapLongZipped(Index, options, argc, argv)){
      return 1;
    }
  }else{
    std::cerr << "\nUnknown readfile filename extension '" << filenameExtension <<  "'. Supported file extensions: fq, fastq, fq.gz or fastq.gz\n";
    return 1;
  }

  std::cerr << "Reads mapped sucessfully!\n\n";

  return 0;
}
