# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <map>
# include <iostream>
# include <time.h>
# include <fstream>
# include "./src/functions.h"
using namespace seqan;

/*
g++ countKmers.cpp -o countK
*/

/*
k can be 31 at max because of the hash function
*/

int main(int argc, char *argv[]){

  // auto begin = std::chrono::high_resolution_clock::now();

  // for (int a=0;a<100;a++){

  if(argc!=5){
    std::cerr << "Usage: ./countK InputFILE k bucket_number Index_name \n\n";
    exit(-1);
  }

// reading the FastQ file

  StringSet<CharString> ids;
  StringSet<Dna5String> seqs;

  try {
    SeqFileIn file(argv[1]);

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
  // defining key parameters

  uint_fast8_t k=std::stoi(argv[2]); // length of k-mer

  // int64_t maxhash=pow(2,k*2)-1;

  int64_t maxhash;
  for (uint_fast8_t i=0;i<k;i++){
    maxhash= maxhash << 2 | 3;
  }

  std::srand(0);
  int64_t random_seed=0;
  for (uint_fast8_t i=0;i<k;++i){
    random_seed= random_seed << 2 | (int64_t)(std::rand()%3);
  }

  uint_fast32_t bucket_number=std::stoul(argv[3]); // should depend on k and the length of the indexed sequence

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

  std::string IndPos=argv[4];
  IndPos.append("_pos.txt");
  std::string IndDir=argv[4];
  IndDir.append("_dir.txt");
  std::string IndC=argv[4];
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
