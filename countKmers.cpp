# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <map>
# include <iostream>
# include <time.h>
# include <fstream>
using namespace seqan;

/*
g++ countKmers.cpp -o countK
*/
std::pair <unsigned,unsigned> hashkMer(const DnaString & kmer, const unsigned k);
std::pair <unsigned,unsigned> rollinghashkMer(unsigned & oldHash, unsigned & oldHash2, const Dna & newnuc, const unsigned k);
unsigned  GetBkt(const unsigned & hash, const String<unsigned> & C, const unsigned long long bucket_number);
unsigned  ReqBkt(const unsigned & hash, String<unsigned> & C, const unsigned long long bucket_number);
std::vector<std::pair <unsigned,unsigned>> RetPos(const DnaString & kmer, const String<unsigned> & C,const String<unsigned> & dir,const String<std::pair <unsigned,unsigned>> & pos, const unsigned long long bucket_number);

int main(int argc, char *argv[]){

  // auto begin = std::chrono::high_resolution_clock::now();

  // for (int a=0;a<100;a++){

  if(argc!=5){
    std::cerr << "Usage: ./countK InputFILE k bucket_number \n\n";
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

  unsigned long long bucket_number=std::stoll(argv[3]); // should depend on k and the length of the indexed sequence

  // choosing Chromosome

  // Dna5String seq=seqs[36];
  //
  // std::cerr << "ID: " << ids[36] << "\n";
  //
  // std::cerr << "Chromosome lengths: " << length(seq) << "\n";

  // building index storage

  String<unsigned> dir;
  resize(dir,bucket_number,0);
  String<std::pair <unsigned,unsigned>> pos;
  resize(pos,length(concat(seqs)));
  String<unsigned> C;
  resize(C,bucket_number,-1);

  typedef Iterator<String<unsigned>>::Type Titrs;

  std::vector<unsigned>::iterator itrv;
  std::vector<unsigned>::reverse_iterator itrvr;

  unsigned long long c;
  unsigned CHROM =0;

  std::cerr << "Index prepared. \n";

  // iterating over the stringSet (Chromosomes)
  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  for (TStringSetIterator seq = begin(seqs); seq != end(seqs); ++seq){

    // counting k-mers

    std::pair<unsigned, unsigned> hash=hashkMer(infix(*seq,0,k),k);    // calculation of the hash value for the first k-mer


    for (unsigned i = 0;i<length(*seq)-k;++i){
      c=ReqBkt(std::min(hash.first,hash.second),C,bucket_number);     // indexing the hashed k-mers
      dir[c+1]+=1;
      if ((*seq)[i+k]!='N'){                                             // calculation of the hash value for the next k-mer
        hash=rollinghashkMer(hash.first,hash.second,(*seq)[i+k],k);
      }
      else {                                                          // reinitialization of the hashvalue after encountering an "N"
        i+=k+1;
        hash=hashkMer(infix(*seq,i,i+k),k);
      }
    }
    c=ReqBkt(std::min(hash.first,hash.second),C,bucket_number);       // indexing of the last element
    dir[c+1]+=1;

  }

  std::cerr << "Index initially filled. \n";

  // cumulative sum

  unsigned sum=length(concat(seqs))-k+1;
  std::vector<unsigned> abundance; //tracking k-mer abundances

  for (Titrs itrs=end(dir)-1;itrs!=begin(dir)-1;--itrs){
    if (*itrs!=0){   //tracking k-mer abundances
      sum-=*itrs;
    abundance.push_back(*itrs);} //tracking k-mer abundances
    *itrs=sum;
  }

  std::cerr << "cumulated sum culculation finished. \n";

  // iterating over the stringSet (Chromosomes)

  for (TStringSetIterator seq = begin(seqs); seq != end(seqs); ++seq){

    // filling pos

    std::pair<unsigned, unsigned> hash=hashkMer(infix(*seq,0,k),k);                                // calculation of the hash value for the first k-mer

    for (unsigned i = 0;i<length(*seq)-k;++i){
      c=GetBkt(std::min(hash.first,hash.second),C,bucket_number);   // filling of the position table
      pos[dir[c+1]]=std::make_pair(CHROM,i);
      dir[c+1]++;
      if ((*seq)[i+k]!='N'){                                           // calculation of the hash value for the next k-mer
        hash=rollinghashkMer(hash.first,hash.second,(*seq)[i+k],k);
      }
      else {                                                        // reinitialization of the hashvalue after encountering an "N"
        i+=k+1;
        hash=hashkMer(infix(*seq,i,i+k),k);
      }
    }
    c=GetBkt(std::min(hash.first,hash.second),C,bucket_number);     // filling the position table for the last element
    pos[dir[c+1]]=std::make_pair(CHROM,length(*seq)-k);
    dir[c+1]++;
    CHROM++;
  }

  std::cerr << "Index build. \n";

  // write index to file

  String<std::pair <unsigned,unsigned>, External<ExternalConfigLarge<>> > extpos;
  if (!open(extpos, "index_pos.txt", OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index counts file." );
  }
  assign(extpos, pos, Exact());
  close(extpos);

  String<unsigned, External<> > extdir;
  if (!open(extdir, "index_dir.txt", OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index counts file." );
  }
  assign(extdir, dir, Exact());
  close(extdir);

  String<unsigned, External<> > extC;
  if (!open(extC, "index_C.txt", OPEN_WRONLY | OPEN_CREATE)){
    throw std::runtime_error("Could not open index counts file." );
  }
  assign(extC, C, Exact());
  close(extC);

  std::cerr << "Index writen to file.\n";

  // Kontrollausgabe


  DnaString testDNA=argv[4];
  std::vector<std::pair <unsigned,unsigned>> positions=RetPos(testDNA,C,dir,pos,bucket_number);
  std::vector<std::pair <unsigned,unsigned>>::iterator itrpv;
  for (itrpv=positions.begin();itrpv!=positions.end();itrpv++){
    std::cout << (*itrpv).first << " " << (*itrpv).second <<"\n";
  }
  std::cout << "\n";

// }

// auto end = std::chrono::high_resolution_clock::now();
// std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();// << "ns" << std::endl;

}

// return k-mer positions
std::vector<std::pair <unsigned,unsigned>> RetPos(const DnaString & kmer, const String<unsigned> & C,const String<unsigned> & dir,const String<std::pair <unsigned,unsigned>> & pos, const unsigned long long bucket_number){
      std::vector<std::pair <unsigned,unsigned>> positions;
      std::pair <unsigned,unsigned> hash=hashkMer(kmer,length(kmer));
      int c=GetBkt(std::min(hash.first,hash.second),C,bucket_number);
      for (unsigned i = dir[c];i!=dir[c+1];i++){
        positions.push_back(pos[i]);
      }
      return positions;
}

// Find correct Bucket
unsigned  GetBkt(const unsigned & hash, const String<unsigned> & C, const unsigned long long bucket_number){
  std::srand(hash);
  unsigned long long i=std::rand()%bucket_number;
  unsigned d=0;
  // unsigned counter=0;
  while(C[i]!=hash and C[i]!=-1){
    // counter+=1;
    i=(i+2*d+1)%bucket_number;
    d++;
    // if (counter > 1000){   // error if bucket_number not high enough
      // std::cerr<<"\nERROR: Bucket number to small.\n";
      // break;}
  }
  return i;
}

// Request a Bucket
unsigned  ReqBkt(const unsigned & hash, String<unsigned> & C, const unsigned long long bucket_number){
  unsigned long long i = GetBkt(hash,C,bucket_number);
  C[i]=hash;
  return i;
}

//  Hashfunction for k-mer
std::pair <unsigned,unsigned> hashkMer(const DnaString & kmer, const unsigned k){
  unsigned hash=0;
  unsigned hash2=0;
  for (int i=0;i<k;++i){
    hash= hash << 2 | ordValue(kmer[i]);
    hash2= hash2 << 2 | (3-ordValue(kmer[k-1-i]));
  }
  return std::make_pair(hash,hash2);
}

// Rolling hashfunction for k-mer
std::pair <unsigned,unsigned> rollinghashkMer(unsigned & oldHash, unsigned & oldHash2, const Dna & newnuc, const unsigned k){
  oldHash=((oldHash << 2) | ordValue(newnuc)) % ((unsigned long)2 << (k*2-1));
  oldHash2=(oldHash2 >> 2) | (3-ordValue(newnuc)) << (k*2-2);
  return std::make_pair(oldHash,oldHash2);
}
