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
# include "IO_functions.h"
using namespace seqan;


int writeIndex(std::string & kmer_index_name, openAddressingKmerHashtable & Index);

float getLoadFactor(String<int32_t> & C){
  uint64_t count=0;
  #pragma omp parallel
  {
    uint64_t localCount = 0;
    #pragma omp for
    for(int i=0;i<length(C);i++){
      localCount+=1*(C[i]>=0);
    }
    #pragma omp atomic
    count+=localCount;
  }
  float loadFactor=(float)count/(float)length(C);
  return loadFactor;
}

int index(int argc, char const **argv){

  // parsing command line arguments
  indexOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine_index(options, argc, argv);

  if (res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;

  printParseResults_index(options);

  //define key parameters

  std::cerr << "Loading reference genome...";

  StringSet<Dna5String> seqs;
  if(loadReference(seqs, options)){return 1;}

  std::cerr << "..done.\n";

  uint32_t bucket_number=round((double)length(concat(seqs))/((double)(options.m-options.k+1)/2));
  if (bucket_number>3221225472){
    bucket_number=3221225472;
  }
  // std::cerr << "Bucket number set to: " << bucket_number << "\n";

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
    std::cerr << "\nError for index target:  " << strerror(errno) << "\n";
    std::cerr << "Old Index will be overwriten\n";
  }

  std::fstream output;
  std::string IndFai=options.kmer_index_name;
  IndFai.append("/fai.txt");
  output.open(toCString(IndFai),std::ios::out);
  for (int id=0; id<numSeqs(faiIndex); id++){
    output << sequenceName(faiIndex, id) << "\t" << sequenceLength(faiIndex, id) <<"\n";
  }
  output.close();

  std::cerr << "...........done.\n";
  std::cerr << "Preparing index...";

  // building index storage
  openAddressingKmerHashtable Index(options.k, options.m, bucket_number);

  std::cerr << "...........done.\nFilling index initially...";
  // iterating over the stringSet (Chromosomes)

  uint64_t random_seed = getRandSeed(options.k);
  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  TStringSetIterator seqG = begin(seqs);
  uint64_t sum=0;
  while(true){
    #pragma omp parallel
    {
      #pragma omp for schedule(dynamic)
      for (int j=0; j<(int)length(seqs); j++){ // iterating over chromosomes/contigs
        uint32_t c;
        uint64_t local_sum=0;
        TStringSetIterator seq=seqG+j;
        // counting k-mers
        minimizedSequence miniSeq(options.k, options.m, random_seed);
        miniSeq.init(*seq);
        minimizer mini;

        while(!miniSeq.at_end){ // iterating through the sequence
          mini=miniSeq.pop();
          c=Index.ReqBkt(mini.value^miniSeq.random_seed);     // indexing the hashed k-mers
          #pragma omp atomic
          Index.dir[c+1]+=1;
          local_sum++;
        }
        #pragma omp atomic
        sum+=local_sum;
      } // for
    } // pragma omp parallel

    // check load factor
    float loadFactor=getLoadFactor(Index.C);
    // std::cerr << "\nLoadFactor: " << loadFactor << "\n";
    if(loadFactor < 0.8){
      sum=0;
      seqG = begin(seqs);
      uint32_t new_bucket_number=round((double)Index.bucket_number*(double)loadFactor/0.85);
      Index.clear();
      Index.resize(new_bucket_number);
      // std::cerr << "\nload_factor recalculated!\n";
    }else{
      break;
    }

  } // while(true)

  std::cerr << "...done. \n";
  std::cerr << "Calculating cumulated sum...";

  // calculate newSum of entrys in dir which are smaller than "a"
  // only substract from newSum if *itrs < "a"
  // only insert kmer positions if dir[c+1]-dir[c] > 0


// ##################################################################

  resize(Index.pos,sum);
  resize(Index.ref,sum);
  typedef Iterator<String<uint32_t>>::Type Titrs;
  for (Titrs itrs=seqan::end(Index.dir)-1; itrs!=seqan::begin(Index.dir)-1; --itrs){
    if (*itrs!=0){   //tracking k-mer abundances
      sum-=(uint64_t)*itrs;
    }
    *itrs=sum;
  }

  std::cerr << ".done.\n";
  std::cerr << "Writing positions to index...";

  seqG = begin(seqs);

  for (int j=0; j<(int)length(seqs); j++){
    uint32_t c;
    TStringSetIterator seq=seqG+j;
    uint8_t Chromosome=j;
    // filling pos
    minimizedSequence miniSeq(options.k, options.m, random_seed);
    miniSeq.init(*seq);
    minimizer mini;

    while(!miniSeq.at_end){
      mini=miniSeq.pop();
      c=Index.GetBkt(mini.value^miniSeq.random_seed);     // indexing the hashed k-mers
      Index.pos[Index.dir[c+1]]=mini.position;
      Index.ref[Index.dir[c+1]]=Chromosome;
      Index.dir[c+1]++;
    }
  }
  std::cerr << "done. \n";

  std::cerr << "Writing index to file...";

  if(writeIndex(options.kmer_index_name, Index)){
    return 1;
  }

  std::cerr << ".....done.\n";
  std::cerr << "Index finished!\n";
  return 0;
}

int writeIndex(std::string & kmer_index_name, openAddressingKmerHashtable & Index){
  std::string IndPos=kmer_index_name;
  IndPos.append("/pos.txt");
  std::string IndRef=kmer_index_name;
  IndRef.append("/ref.txt");
  std::string IndDir=kmer_index_name;
  IndDir.append("/dir.txt");
  std::string IndC=kmer_index_name;
  IndC.append("/C.txt");
  std::string InfoFileName=kmer_index_name;
  InfoFileName.append("/Info.txt");

  try{
    String<uint32_t, External<ExternalConfigLarge<>> > extpos;
    if (!open(extpos, IndPos.c_str(), OPEN_WRONLY | OPEN_CREATE)){
      throw std::runtime_error("Could not open index positions file." );
    }
    assign(extpos, Index.pos, Exact());
    close(extpos);

    String<uint8_t, External<ExternalConfigLarge<>> > extref;
    if (!open(extref, IndRef.c_str(), OPEN_WRONLY | OPEN_CREATE)){
      throw std::runtime_error("Could not open index reference file." );
    }
    assign(extref, Index.ref, Exact());
    close(extref);

    String<uint32_t, External<> > extdir;
    if (!open(extdir, IndDir.c_str(), OPEN_WRONLY | OPEN_CREATE)){
      throw std::runtime_error("Could not open index directory file." );
    }
    assign(extdir, Index.dir, Exact());
    close(extdir);

    String<int32_t, External<> > extC;
    if (!open(extC, IndC.c_str(), OPEN_WRONLY | OPEN_CREATE)){
      throw std::runtime_error("Could not open index counts file." );
    }
    assign(extC, Index.C, Exact());
    close(extC);
  }
  catch(std::runtime_error e){
    std::cerr << e.what();
  }

  std::ofstream InfoFile(InfoFileName);
  InfoFile << Index.k << " " << Index.m;
  InfoFile.close();

  return 0;
}
