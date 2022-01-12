# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <map>
# include <iostream>
# include <time.h>
# include <fstream>
# include "functions.h"
# include "index.h"

using namespace seqan;

/*
g++ countKmers.cpp -o countK
*/

struct countKOptions{
  std::string reference_file;
  std::string kmer_index_name;
  unsigned k;
  unsigned m;//int64_t bucket_count;

  countKOptions():
  k(31), m(61)// bucket_count(3221225472)
  {}
};

seqan::ArgumentParser::ParseResult parseCommandLine(countKOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bcmap index");

    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "reference(.fastq/.fasta)"));
    // addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Index_name[OUT]"));

    // Define Options
    addOption(parser, seqan::ArgParseOption(
        "o", "kmer_index_name", "Name of the folder in which the kmer index is stored.",
        seqan::ArgParseArgument::STRING, "kmer_index_name[OUT]"));
    setDefaultValue(parser, "o", "Index");
    addOption(parser, seqan::ArgParseOption(
        "k", "kmer_length", "Length of kmers in index.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "k", "31");
    setMinValue(parser, "k", "8");
    setMaxValue(parser, "k", "31");
    addOption(parser, seqan::ArgParseOption(
        "m", "minimizer_window", "Length of window a minimizer is chosen from.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "m", "61");
    // addOption(parser, seqan::ArgParseOption(
    //     "b", "bucket_count", "Number of buckets in index.",
    //     seqan::ArgParseArgument::INT64, "unsigned"));
    // setDefaultValue(parser, "b", "3221225472");
    seqan::addUsageLine(parser,"reference.fq [OPTIONS]");
    setShortDescription(parser, "Build an index of a reference genome.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addDescription(parser,"Builds an open adressing k-mer index for the given reference genome(fastq/fasta).");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;}

    // Extract option values.
    getOptionValue(options.k, parser, "k");
    getOptionValue(options.m, parser, "m");
    // getOptionValue(options.bucket_count, parser, "b");
    getOptionValue(options.kmer_index_name, parser, "o");

    getArgumentValue(options.reference_file, parser, 0);

    return seqan::ArgumentParser::PARSE_OK;
}

int index(int argc, char const **argv){

  // parsing command line arguments
  countKOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

  if (res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;
  std::cout <<'\n'
            << "reference        \t" << options.reference_file << '\n'
            << "kmer_index_name  \t" << options.kmer_index_name << '\n'
            << "k                \t" << options.k << '\n'
            << "minimizer_window \t" << options.m << "\n\n";
            // << "bucket_count     \t" << options.bucket_count << "\n\n";

  uint_fast8_t k = options.k;
  int k_2 = k+1;
  uint_fast8_t m = options.m;
  // options.bucket_count; // should depend on k and the length of the indexed sequence

  // auto begin = std::chrono::high_resolution_clock::now();

// reading the FastQ file

  StringSet<CharString> ids;
  StringSet<Dna5String> seqs;
  std::cerr << "Loading reference genome...";
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

  std::cerr << "..done.\n";

  uint_fast32_t bucket_number=round((double)length(concat(seqs))/((double)(m-k+1)/2));
  if (bucket_number>3221225472){
    bucket_number=3221225472;
  }

  std::cerr << "Bucket number set to: " << bucket_number << "\n";

  std::cerr << "Loading ref.fai...";

  FaiIndex faiIndex;
  if (!open(faiIndex, toCString(options.reference_file))){
    std::cerr << ".........failed.\nBuilding ref.fai..";
    if (!build(faiIndex, toCString(options.reference_file))){
        std::cerr << "\nERROR: FASTA index could not be loaded or built.\n";
      return 1;
    }
    if (!save(faiIndex)) // Name is stored from when reading.
      {
        std::cerr << "\nWARNING: FASTA index could not be written to disk.\n";
      }
  }

  if (mkdir(toCString(options.kmer_index_name), 0777) == -1){
    // std::cerr << "Error for index target:  " << strerror(errno) << "\n";
  }

  std::fstream output;
  std::string IndFai=options.kmer_index_name;
  IndFai.append("/fai.txt");
  output.open(toCString(IndFai),std::ios::out);
  for (int id=0; id<numSeqs(faiIndex); id++){
    output << sequenceName(faiIndex, id) << "\n";
  }
  output.close();

  std::cerr << "...........done.\n";
  std::cerr << "Preparing index...";

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
  String<uint32_t> pos;
  String<uint_fast8_t> ref;
  String<int32_t> C;
  resize(dir,bucket_number+1,0);
  std::cerr << "....";
  // resize(pos,length(concat(seqs)));
  // resize(ref,length(concat(seqs)));
  resize(C,bucket_number,-1);
  std::cerr << "....";


  typedef Iterator<String<uint32_t>>::Type Titrs;

  // uint32_t c;

  std::cerr << "...done.\nFilling index initially...";
  // iterating over the stringSet (Chromosomes)
  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  TStringSetIterator seqG = begin(seqs);
  // uint_fast8_t CHROM=0;
  // int laenge=length(seqs);
  uint64_t sum=0;
  #pragma omp parallel
  {
    #pragma omp for schedule(dynamic)
    for (int j=0; j<(int)length(seqs); j++){ // iterating over chromosomes/contigs
      uint32_t c;
      uint64_t local_sum=0;
      TStringSetIterator seq=seqG+j;
      // counting k-mers
      minimizer mini;
      minimizedSequence miniSeq(*seq,k,m,random_seed,maxhash);

      while(!miniSeq.at_end){ // iterating through the sequence
        mini=miniSeq.pop();
          c=ReqBkt(mini.value^random_seed,C,bucket_number,k_2);     // indexing the hashed k-mers
          #pragma omp atomic
          dir[c+1]+=1;
          local_sum++;
      }
      #pragma omp atomic
      sum+=local_sum;
    }
  }

  std::cerr << "...done. \n";
  std::cerr << "Calculating cumulated sum...";

  // cumulative sum

  resize(pos,sum);
  resize(ref,sum);
  // uint64_t sum=length(concat(seqs))-k+1;

  for (Titrs itrs=end(dir)-1;itrs!=begin(dir)-1;--itrs){
    if (*itrs!=0){   //tracking k-mer abundances
      sum-=(uint64_t)*itrs;
    }
    *itrs=sum;
  }

  std::cerr << ".done.\n";

  // iterating over the stringSet (Chromosomes)
  std::cerr << "Writing positions to index...";
  seqG = begin(seqs);
  // #pragma omp parallel
  // {
    // #pragma omp for schedule(dynamic)
    for (int j=0; j<(int)length(seqs); j++){
      uint32_t c;
      TStringSetIterator seq=seqG+j;
      uint_fast8_t Chromosome=j;
      // filling pos
      minimizer mini;
      minimizedSequence miniSeq(*seq,k,m,random_seed,maxhash);
      while(!miniSeq.at_end){
        mini=miniSeq.pop();
        c=GetBkt(mini.value^random_seed,C,bucket_number,k_2);     // indexing the hashed k-mers

        // #pragma omp critical
        // {
        pos[dir[c+1]]=mini.position;
        ref[dir[c+1]]=Chromosome;
        dir[c+1]++;
        // }
      }
    }
  // }
  std::cerr << "done. \n";

  //write index to file

  std::string IndPos=options.kmer_index_name;
  IndPos.append("/pos.txt");
  std::string IndRef=options.kmer_index_name;
  IndRef.append("/ref.txt");
  std::string IndDir=options.kmer_index_name;
  IndDir.append("/dir.txt");
  std::string IndC=options.kmer_index_name;
  IndC.append("/C.txt");

  std::cerr << "Writing index to file...";

  String<uint32_t, External<ExternalConfigLarge<>> > extpos;
  if (!open(extpos, IndPos.c_str(), OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index positions file." );
  }
  assign(extpos, pos, Exact());
  close(extpos);

  String<uint_fast8_t, External<ExternalConfigLarge<>> > extref;
  if (!open(extref, IndRef.c_str(), OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index reference file." );
  }
  assign(extref, ref, Exact());
  close(extref);

  String<uint32_t, External<> > extdir;
  if (!open(extdir, IndDir.c_str(), OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index directory file." );
  }
  assign(extdir, dir, Exact());
  close(extdir);

  String<int32_t, External<> > extC;
  if (!open(extC, IndC.c_str(), OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index counts file." );
  }
  assign(extC, C, Exact());
  close(extC);

  std::cerr << ".....done.\n";
  std::cerr << "Index finished!\n";
  return 0;
}
