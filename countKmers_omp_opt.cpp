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
  int64_t bucket_count;
  uint32_t thread_count;

  countKOptions():
  k(31), bucket_count(3221225472), thread_count(32)
  {}
  };

seqan::ArgumentParser::ParseResult parseCommandLine(countKOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("countK");

    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "Path to reference.(fastq/fasta)"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Index_name[OUT]"));

    // Define Options
    addOption(parser, seqan::ArgParseOption(
        "k", "kmer_length", "Length of kmers in index.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "k", "31");
    setMinValue(parser, "k", "8");
    setMaxValue(parser, "k", "31");
    addOption(parser, seqan::ArgParseOption(
        "b", "bucket_count", "Number of buckets in index.",
        seqan::ArgParseArgument::INT64, "unsigned"));
    setDefaultValue(parser, "b", "3221225472");
    addOption(parser, seqan::ArgParseOption(
        "t", "thread_count", "Number of threads available.",
        seqan::ArgParseArgument::INT64, "unsigned"));
    setDefaultValue(parser, "t", "32");

    setShortDescription(parser, "Build index of reference genome.");
    setVersion(parser, "0.1");
    setDate(parser, "March 24 2021");
    addDescription(parser,"Builds an open adressing k-mer index for the given reference genome(fastq/fasta).");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;}

    // Extract option values.
    getOptionValue(options.k, parser, "k");
    getOptionValue(options.bucket_count, parser, "b");
    getOptionValue(options.thread_count, parser, "t");

    // Extract argument calues.
    getArgumentValue(options.reference_file, parser, 0);
    getArgumentValue(options.index_name, parser, 1);

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const **argv){

  // parsing command line arguments
  countKOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

  if (res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;
  std::cout <<'\n'
            << "reference   \t" << options.reference_file << '\n'
            << "index_name  \t" << options.index_name << '\n'
            << "k           \t" << options.k << '\n'
            << "bucket_count\t" << options.bucket_count <<'\n'
            << "thread_count\t" << options.thread_count << "\n\n";

  uint_fast8_t k = options.k;
  int k_2 = k+1;
  uint_fast32_t bucket_number=options.bucket_count; // should depend on k and the length of the indexed sequence

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
  std::cerr << "Preparing index...";

  // split large chromosomes into smaller contigs
  std::vector<uint32_t> Chromtable={}; // vector that translates place in seqs into chromosome identifier

  std::cerr << "\n";
  for (int i=0; i<length(seqs); i++){
    Chromtable.push_back(i);
    std::cerr << "chrom: " << i << "\tsize: " << length(seqs[i]) << "\n";
  }


  int maxseqlen=50000000;
  StringSet<Dna5String> seqs2;


  Chromtable={};
  for(int i=0; i<length(seqs); i++){
    if (length(seqs[i])>maxseqlen){
      for (int j=0; j<(length(seqs[i])%maxseqlen);j++){
        appendValue(seqs2,infix(seqs[i],j*maxseqlen,(j+1)*maxseqlen));
        Chromtable.push_back(i);
      }
      appendValue(seqs2,suffix(seqs[i],(length(seqs[i])%maxseqlen)*maxseqlen));
    }
  }


  //
  // if(options.thread_count>1){
  //     int seqpos=0;
  //     while (seqpos<length(seqs)){
  //       if (length(seqs[seqpos])>maxseqlen){
  //         insertValue(seqs,seqpos,suffix(seqs[seqpos],maxseqlen));
  //         seqs[seqpos]=prefix(seqs[seqpos],maxseqlen);
  //         Chromtable.insert(Chromtable.begin()+seqpos,Chromtable[seqpos]);
  //         seqpos++;
  //       }else{
  //         seqpos++;
  //       }
  //     }
  // }

  std::cerr << "\n";
  for (int i=0; i<length(seqs); i++){
    Chromtable.push_back(i);
    std::cerr << "chrom: " << Chromtable[i] << "\tsize: " << length(seqs[i]) << "\n";
  }



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
  std::cerr << "..";
  resize(pos,length(concat(seqs)));
  std::cerr << "..";
  resize(ref,length(concat(seqs)));
  std::cerr << "..";
  resize(C,bucket_number,-1);
  std::cerr << "..";


  typedef Iterator<String<uint32_t>>::Type Titrs;

  uint32_t c;
  uint_fast8_t CHROMG = 0;

  std::cerr << "...done.\nFilling index initially:";
  // iterating over the stringSet (Chromosomes)
  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  TStringSetIterator seqG = begin(seqs);
  uint_fast8_t CHROM=0;
  // int laenge=length(seqs);
  #pragma omp parallel
  {
    #pragma omp for schedule(dynamic)
    for (int i=0; i<(int)length(seqs); i++){
      TStringSetIterator seq=seqG+i;
      // counting k-mers
      std::pair<int64_t, int64_t> hash=hashkMer(infix(*seq,0,k),k);    // calculation of the hash value for the first k-mer

      for (uint64_t i = 0;i<length(*seq)-k;++i){
        c=ReqBkt(ReturnSmaller(hash.first,hash.second,random_seed),C,bucket_number,k_2);     // indexing the hashed k-mers
        #pragma omp atomic
        dir[c+1]+=1;
        if ((*seq)[i+k]!='N'){                                             // calculation of the hash value for the next k-mer
          rollinghashkMer(hash.first,hash.second,(*seq)[i+k],k,maxhash);
        }
        else {                                                          // reinitialization of the hashvalue after encountering an "N"
          i+=k+1;
          hash=hashkMer(infix(*seq,i,i+k),k);
        }
      }

      c=ReqBkt(ReturnSmaller(hash.first,hash.second,random_seed),C,bucket_number,k_2);       // indexing of the last element
      #pragma omp atomic
      dir[c+1]+=1;
      #pragma omp atomic
      CHROM++;
      std::cerr << "." ;
      if ((CHROM-5)%29==0) {std::cerr << "\n";}
    }
  }

  std::cerr << "done. \n";
  std::cerr << "Calculating cumulated sum...";

  // cumulative sum

  uint64_t sum=length(concat(seqs))-k+1;

  for (Titrs itrs=end(dir)-1;itrs!=begin(dir)-1;--itrs){
    if (*itrs!=0){   //tracking k-mer abundances
      sum-=(uint64_t)*itrs;
    }
    // abundance.push_back(*itrs);} //tracking k-mer abundances
    *itrs=sum;
  }

  std::cerr << ".done.\n";

  // iterating over the stringSet (Chromosomes)
  std::cerr << "Writing positions to index:";
  seqG = begin(seqs);
  #pragma omp parallel
  {
    #pragma omp for schedule(dynamic)
    for (int i=0; i<(int)length(seqs); i++){
      TStringSetIterator seq=seqG+i;
      CHROM=i;
      // filling pos

      std::pair<int64_t, int64_t> hash=hashkMer(infix(*seq,0,k),k);                                // calculation of the hash value for the first k-mer

      for (uint64_t i = 0;i<length(*seq)-k;++i){
        c=GetBkt(ReturnSmaller(hash.first,hash.second,random_seed),C,bucket_number,k_2);   // filling of the position table
        pos[dir[c+1]]=i;
        ref[dir[c+1]]=Chromtable[i];
        dir[c+1]++;
        if ((*seq)[i+k]!='N'){                                           // calculation of the hash value for the next k-mer
          rollinghashkMer(hash.first,hash.second,(*seq)[i+k],k,maxhash);
        }
        else {                                                        // reinitialization of the hashvalue after encountering an "N"
          i+=k+1;
          hash=hashkMer(infix(*seq,i,i+k),k);
        }
      }
      c=GetBkt(ReturnSmaller(hash.first,hash.second,random_seed),C,bucket_number,k_2);     // filling the position table for the last element
      pos[dir[c+1]]=length(*seq)-k;
      ref[dir[c+1]]=Chromtable[i];
      dir[c+1]++;

      std::cerr << ".";
      if ((CHROM-2)%29==0) {std::cerr << "\n";}
    }
  }
  std::cerr << "done. \n";

  //write index to file

  std::string IndPos=options.index_name;
  IndPos.append("_pos.txt");
  std::string IndRef=options.index_name;
  IndRef.append("_ref.txt");
  std::string IndDir=options.index_name;
  IndDir.append("_dir.txt");
  std::string IndC=options.index_name;
  IndC.append("_C.txt");

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
