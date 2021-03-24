# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <map>
# include <iostream>
# include <time.h>
# include <fstream>
# include "./src/functions.h"
using namespace seqan;

/*
g++ countKmers.cpp -o countK
*/

struct countKOptions{
  std::string reference_file;
  std::string index_name;
  unsigned k;
  long unsigned bucket_count;

  countKOptions():
  k(31), bucket_count(3221225472)
  {}
  };

seqan::ArgumentParser::ParseResult parseCommandLine(countKOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("countK");

    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "reference_file"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "index_name"));

    // Define Options
    addOption(parser, seqan::ArgParseOption(
        "k", "kmer_length", "length of k-mers in kmer indexed",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "k", "31");
    addOption(parser, seqan::ArgParseOption(
        "b", "bucket_count", "number of buckets in index",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "b", "3221225472");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;}

    // Extract option values.
    getOptionValue(options.k, parser, "k");
    getOptionValue(options.bucket_count, parser, "b");

    getArgumentValue(options.reference_file, parser, 0);
    getArgumentValue(options.index_name, parser, 1);

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char *argv[]){

  // parsing command line arguments
  countKOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

  // If parsing was not successful then exit with code 1 if there were errors.
  // Otherwise, exit with code 0 (e.g. help was printed).
  if (res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;
  std::cout <<'\n'
            << "k                \t" << options.k << '\n'
            << "bucket_count     \t" << options.bucket_count << '\n'
            << "reference        \t" << options.reference_file << '\n'
            << "index_name       \t" << options.index_name << "\n\n";

  uint_fast8_t k = options.k;
  uint_fast32_t bucket_number=options.bucket_count; // should depend on k and the length of the indexed sequence

  // auto begin = std::chrono::high_resolution_clock::now();

// reading the FastQ file

  StringSet<CharString> ids;
  StringSet<Dna5String> seqs;

  try {
    SeqFileIn file(toCString(options.reference_file));

    readRecords(ids, seqs, file);

    close(file);
  }
  catch (ParseError const & e){
    std::cerr << "ERROR: input record is badly formatted. " << e.what() << std::endl;
  }
  catch (IOError const & e){
    std::cerr << "ERROR: input file can not be opened. " << e.what() << std::endl;
  }

  std::cerr << "Genome read in. \n";

  int64_t maxhash;
  for (uint_fast8_t i=0;i<k;i++){
    maxhash= maxhash << 2 | 3;
  }

  std::srand(0);
  int64_t random_seed=0;
  for (uint_fast8_t i=0;i<k;++i){
    random_seed= random_seed << 2 | (int64_t)(std::rand()%3);
  }


  // building index storage

  String<uint32_t> dir;
  resize(dir,bucket_number+1,0);
  String<std::pair <uint_least8_t,uint32_t>> pos;
  resize(pos,length(concat(seqs)));   // may be re
  String<int64_t> C;
  resize(C,bucket_number,-1);

  typedef Iterator<String<uint32_t>>::Type Titrs;

  uint64_t c;
  uint_least8_t CHROM = 0;


  std::cerr << "Index prepared. \n";
  // iterating over the stringSet (Chromosomes)
  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  for (TStringSetIterator seq = begin(seqs); seq != end(seqs); ++seq){
    std::cerr << "Chrom: " << (int)CHROM << "\n";
    // counting k-mers

    std::pair<int64_t, int64_t> hash=hashkMer(infix(*seq,0,k),k);    // calculation of the hash value for the first k-mer

    for (uint64_t i = 0;i<length(*seq)-k;++i){
      c=ReqBkt(ReturnSmaller(hash.first,hash.second,random_seed),C,bucket_number);     // indexing the hashed k-mers
      dir[c+1]+=1;
      if ((*seq)[i+k]!='N'){                                             // calculation of the hash value for the next k-mer
        rollinghashkMer(hash.first,hash.second,(*seq)[i+k],k,maxhash);
      }
      else {                                                          // reinitialization of the hashvalue after encountering an "N"
        i+=k+1;
        hash=hashkMer(infix(*seq,i,i+k),k);
      }
    }
    c=ReqBkt(ReturnSmaller(hash.first,hash.second,random_seed),C,bucket_number);       // indexing of the last element
    dir[c+1]+=1;
    CHROM++;
  }

  std::cerr << "Index initially filled. \n";


  // cumulative sum

  uint64_t sum=length(concat(seqs))-k+1;

  for (Titrs itrs=end(dir)-1;itrs!=begin(dir)-1;--itrs){
    if (*itrs!=0){   //tracking k-mer abundances
      sum-=(uint64_t)*itrs;
    }
    // abundance.push_back(*itrs);} //tracking k-mer abundances
    *itrs=sum;
  }

  std::cerr << "cumulated sum culculation finished. \n";

  // iterating over the stringSet (Chromosomes)
  CHROM=0;
  for (TStringSetIterator seq = begin(seqs); seq != end(seqs); ++seq){

    // filling pos

    std::pair<int64_t, int64_t> hash=hashkMer(infix(*seq,0,k),k);                                // calculation of the hash value for the first k-mer

    for (uint64_t i = 0;i<length(*seq)-k;++i){
      c=GetBkt(ReturnSmaller(hash.first,hash.second,random_seed),C,bucket_number);   // filling of the position table
      pos[dir[c+1]]=std::make_pair(CHROM,i);
      dir[c+1]++;
      if ((*seq)[i+k]!='N'){                                           // calculation of the hash value for the next k-mer
        rollinghashkMer(hash.first,hash.second,(*seq)[i+k],k,maxhash);
      }
      else {                                                        // reinitialization of the hashvalue after encountering an "N"
        i+=k+1;
        hash=hashkMer(infix(*seq,i,i+k),k);
      }
    }
    c=GetBkt(ReturnSmaller(hash.first,hash.second,random_seed),C,bucket_number);     // filling the position table for the last element
    pos[dir[c+1]]=std::make_pair(CHROM,length(*seq)-k);
    dir[c+1]++;
    CHROM++;
  }

  std::cerr << "Index build. \n";

  //write index to file

  std::string IndPos=options.index_name;
  IndPos.append("_pos.txt");
  std::string IndDir=options.index_name;
  IndDir.append("_dir.txt");
  std::string IndC=options.index_name;
  IndC.append("_C.txt");


  String<std::pair <uint_least8_t,uint32_t>, External<ExternalConfigLarge<>> > extpos;
  if (!open(extpos, IndPos.c_str(), OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index counts file." );
  }
  assign(extpos, pos, Exact());
  close(extpos);

  String<uint32_t, External<> > extdir;
  if (!open(extdir, IndDir.c_str(), OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index counts file." );
  }
  assign(extdir, dir, Exact());
  close(extdir);

  String<int64_t, External<> > extC;
  if (!open(extC, IndC.c_str(), OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index counts file." );
  }
  assign(extC, C, Exact());
  close(extC);

  std::cerr << "Index writen to file.\n";

}
