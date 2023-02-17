# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <map>
# include <iostream>
# include <time.h>
# include <fstream>
# include "functions.h"
# include "parser.h"
# include "index.h"

int index(int argc, char const **argv){

  // parsing command line arguments
  indexOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine_index(options, argc, argv);

  if (res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;

  printParseResults_index(options);

  //define key parameters

  uint_fast8_t k = options.k;
  int k_2;
  if (k>16){
    k_2=(k-15)*2;
  }else{
    k_2=0;
  }

  uint_fast8_t m = options.m;


  // reading the FastQ file

  StringSet<CharString> ids;
  StringSet<IupacString> seqsIn;
  std::cerr << "Loading reference genome...";
  try {
    SeqFileIn file(toCString(options.reference_file));

    readRecords(ids, seqsIn, file);

    close(file);
  }
  catch (ParseError const & e){
    std::cerr << "ERROR: input record is badly formatted. " << e.what() << std::endl;
    return 0;
  }
  catch (IOError const & e){
    std::cerr << "ERROR: input file can not be opened. " << e.what() << std::endl;
    return 0;
  }

  StringSet<Dna5String> seqs=seqsIn; //convert Iupac to Dna5
  clear(seqsIn); //delete Iupac version

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
    std::cerr << "Error for index target:  " << strerror(errno) << "\n";
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

  int64_t maxhash=0;
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
  resize(C,bucket_number,-1);
  std::cerr << "....";


  typedef Iterator<String<uint32_t>>::Type Titrs;

  std::cerr << "...done.\nFilling index initially...";
  // iterating over the stringSet (Chromosomes)

  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  TStringSetIterator seqG = begin(seqs);
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
      pos[dir[c+1]]=mini.position;
      ref[dir[c+1]]=Chromosome;
      dir[c+1]++;
    }
  }
  std::cerr << "done. \n";

  std::cerr << "Writing index to file...";

  std::string IndPos=options.kmer_index_name;
  IndPos.append("/pos.txt");
  std::string IndRef=options.kmer_index_name;
  IndRef.append("/ref.txt");
  std::string IndDir=options.kmer_index_name;
  IndDir.append("/dir.txt");
  std::string IndC=options.kmer_index_name;
  IndC.append("/C.txt");

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
