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
std::pair <unsigned,unsigned> hashkMer(const DnaString & kmer, const unsigned & k);
std::pair <unsigned,unsigned> rollinghashkMer(unsigned & oldHash, unsigned & oldHash2, const Dna & newnuc, const unsigned & k);
unsigned  GetBkt(const unsigned & hash, const std::vector<unsigned> & C, const unsigned bucket_number);
unsigned  ReqBkt(const unsigned & hash, std::vector<unsigned> & C, const unsigned bucket_number);
std::vector<unsigned> RetPos(const DnaString & kmer, const std::vector<unsigned> & C,const std::vector<unsigned> & dir,const std::vector<unsigned> & pos, const unsigned bucket_number);

int main(int argc, char *argv[]){

  // auto begin = std::chrono::high_resolution_clock::now();

  // for (int a=0;a<100;a++){

  if(argc!=4){
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

  // defining key parameters

  unsigned k=std::stoi(argv[2]); // length of k-mer
  unsigned bucket_number=std::stoi(argv[3]); // should depend on k and the length of the indexed sequence

  // concatination of all sequences

  Dna5String seq=concat(seqs);

  // building index storage

  std::vector<unsigned> dir(bucket_number+1,0);       // pow(4,k) depending on k-mer size
  std::vector<unsigned> pos(length(seq),0);         // length(seq)-k+1 runns into error
  std::vector<unsigned> C(bucket_number+1,-1);
  std::vector<unsigned>::iterator itrv;
  std::vector<unsigned>::reverse_iterator itrvr;

  // counting k-mers

  std::pair<unsigned, unsigned> hash=hashkMer(infix(seq,0,k),k);
  unsigned c;

  for (unsigned i = 0;i<length(seq)-k+1;++i){
    c=ReqBkt(std::min(hash.first,hash.second),C,bucket_number);
    dir[c+1]+=1;
    if (seq[i+k]!='N'){
      hash=rollinghashkMer(hash.first,hash.second,seq[i+k],k);
    }
    else {
      i+=k+1;
      hash=hashkMer(infix(seq,i,i+k),k);

    }
  }
  unsigned sum=length(seq)-k+1;

  // cumulative sum

  for (itrvr=dir.rbegin();itrvr!=dir.rend();++itrvr){
    sum-=*itrvr;
    *itrvr=sum;
  }

  // filling pos

  hash=hashkMer(infix(seq,0,k),k);

  for (unsigned i = 0;i<length(seq)-k+1;++i){
    c=GetBkt(std::min(hash.first,hash.second),C,bucket_number);
    pos[dir[c+1]]=i;
    dir[c+1]++;
    if (seq[i+k]!='N'){
      hash=rollinghashkMer(hash.first,hash.second,seq[i+k],k);
    }
    else {
      i+=k+1;
      hash=hashkMer(infix(seq,i,i+k),k);

    }
  }

  // write index to file
  std::ofstream index;
  index.open("index.txt");

  index << bucket_number << " " << k << " " << length(seq) << "\n";

  index << "\n";
  for (itrv=C.begin();itrv!=C.end();itrv++){
    index << *itrv << " ";
  }
  index << "\n";

  index << "\n";
  for (itrv=dir.begin();itrv!=dir.end();itrv++){
    index << *itrv << " ";
  }
  index << "\n";

  index << "\n";
  for (itrv=pos.begin();itrv!=pos.end();itrv++){
    index << *itrv << " ";
  }
  index << "\n";
  index.close();


  // Kontrollausgabe


  DnaString testDNA="TGAAGTGT";
  std::vector<unsigned> positions=RetPos(testDNA,C,dir,pos,bucket_number);
  for (itrv=positions.begin();itrv!=positions.end();itrv++){
    std::cout << *itrv <<" ";
  }
  std::cout << "\n";

// }

// auto end = std::chrono::high_resolution_clock::now();
// std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();// << "ns" << std::endl;

}

// return k-mer positions
std::vector<unsigned> RetPos(const DnaString & kmer, const std::vector<unsigned> & C,const std::vector<unsigned> & dir,const std::vector<unsigned> & pos, const unsigned bucket_number){
      std::vector<unsigned> positions;
      std::pair <unsigned,unsigned> hash=hashkMer(kmer,length(kmer));
      int c=GetBkt(std::min(hash.first,hash.second),C,bucket_number);
      for (unsigned i = dir[c];i!=dir[c+1];i++){
        positions.push_back(pos[i]);
      }
      return positions;
}

// Find correct Bucket
unsigned  GetBkt(const unsigned & hash, const std::vector<unsigned> & C, const unsigned bucket_number){
  int testo=0;
  unsigned i=0;
  try{
  std::srand(hash);
  i=std::rand()%bucket_number;
  unsigned d=0;
  unsigned counter=0;
  while(C[i]!=hash and C[i]!=-1){
    counter+=1;
    i=(i+2*d+1)%bucket_number;
    d++;
    if (counter > bucket_number){   // error if bucket_number not high enough
      std::cerr<<"\nERROR: Bucket number to small.\n";
      break;}
  }
  testo=1;

}
catch(const int exception){} //try bracket
if (testo==0){std::cerr<<"error in GetBkt  ";}
  return i;
}

// Request a Bucket
unsigned  ReqBkt(const unsigned & hash, std::vector<unsigned> & C, const unsigned bucket_number){
  unsigned i = GetBkt(hash,C,bucket_number);
  C[i]=hash;
  return i;
}

//  Hashfunction for k-mer
std::pair <unsigned,unsigned> hashkMer(const DnaString & kmer, const unsigned & k){
  unsigned hash=0;
  unsigned hash2=0;
  for (int i=0;i<k;++i){
    hash= hash << 2 | ordValue(kmer[i]);
    hash2= hash2 << 2 | (3-ordValue(kmer[k-1-i]));
  }
  return std::make_pair(hash,hash2);
}

// Rolling hashfunction for k-mer
std::pair <unsigned,unsigned> rollinghashkMer(unsigned & oldHash, unsigned & oldHash2, const Dna & newnuc, const unsigned & k){
  oldHash=((oldHash << 2) | ordValue(newnuc)) % ((unsigned long)2 << (k*2-1));
  oldHash2=(oldHash2 >> 2) | (3-ordValue(newnuc)) << (k*2-2);
  return std::make_pair(oldHash,oldHash2);
}
