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

  unsigned k=std::stoi(argv[2]); // length of k-mer

  // long long int maxhash=pow(2,k*2)-1;

  long long int maxhash;
  for (unsigned i=0;i<k;i++){
    maxhash= maxhash << 2 | 3;
  }

  std::srand(0);
  long long int random_seed=0;
  for (unsigned i=0;i<k;++i){
    random_seed= random_seed << 2 | (long long int)(std::rand()%3);
  }

  unsigned bucket_number=std::stoul(argv[3]); // should depend on k and the length of the indexed sequence
  // unsigned maxfreq=std::stoll(argv[4]);
  // choosing Chromosome

  // Dna5String seq=seqs[36];
  //
  // std::cerr << "ID: " << ids[36] << "\n";
  //
  // std::cerr << "Chromosome lengths: " << length(seq) << "\n";

  // building index storage

  String<unsigned> dir;
  resize(dir,bucket_number+1,0);
  String<std::pair <unsigned char,unsigned>> pos;
  resize(pos,length(concat(seqs)));   // may be re
  String<int long long> C;
  resize(C,bucket_number,-1);

  typedef Iterator<String<unsigned>>::Type Titrs;

  unsigned long long c;
  unsigned char CHROM =0;


  std::cerr << "Index prepared. \n";
  // iterating over the stringSet (Chromosomes)
  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  for (TStringSetIterator seq = begin(seqs); seq != end(seqs); ++seq){
    std::cerr << "Chrom: " << (int)CHROM << "\n";
    // counting k-mers

    std::pair<long long int, long long int> hash=hashkMer(infix(*seq,0,k),k);    // calculation of the hash value for the first k-mer

    for (long long unsigned i = 0;i<length(*seq)-k;++i){
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

  }

  std::cerr << "Index initially filled. \n";


  // cumulative sum

  long long unsigned sum=length(concat(seqs))-k+1;
  // std::vector<unsigned> abundance; //tracking k-mer abundances

  for (Titrs itrs=end(dir)-1;itrs!=begin(dir)-1;--itrs){
    if (*itrs!=0){   //tracking k-mer abundances
      sum-=(long long unsigned)*itrs;
    }
    // abundance.push_back(*itrs);} //tracking k-mer abundances
    *itrs=sum;
  }

  std::cerr << "cumulated sum culculation finished. \n";

  // iterating over the stringSet (Chromosomes)

  for (TStringSetIterator seq = begin(seqs); seq != end(seqs); ++seq){

    // filling pos

    std::pair<long long int, long long int> hash=hashkMer(infix(*seq,0,k),k);                                // calculation of the hash value for the first k-mer

    for (long long unsigned i = 0;i<length(*seq)-k;++i){
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


  String<std::pair <unsigned char,unsigned>, External<ExternalConfigLarge<>> > extpos;
  if (!open(extpos, IndPos.c_str(), OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index counts file." );
  }
  assign(extpos, pos, Exact());
  close(extpos);

  String<unsigned, External<> > extdir;
  if (!open(extdir, IndDir.c_str(), OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index counts file." );
  }
  assign(extdir, dir, Exact());
  close(extdir);

  String<long long int, External<> > extC;
  if (!open(extC, IndC.c_str(), OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index counts file." );
  }
  assign(extC, C, Exact());
  close(extC);

  std::cerr << "Index writen to file.\n";

}
