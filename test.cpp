# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <map>
# include <iostream>
using namespace seqan;

/*
g++ test.cpp -o testit
*/
std::pair <unsigned,unsigned> hashkMer(const DnaString & kmer, const unsigned & k);
std::pair <unsigned,unsigned> rollinghashkMer(unsigned & oldHash, unsigned & oldHash2, const Dna & newnuc, const unsigned & k);


int main(int argc, char *argv[]){


  if(argc!=2){
    std::cerr << "Usage: ./testit k-mer \n\n";
    exit(-1);
  }

DnaString DNA=argv[1];
Dna newnuc='C';
unsigned k = length(DNA);
std::pair <unsigned,unsigned> hash=hashkMer(DNA,k);

std::cout << "hash1: " << hash.first << "\n";
std::cout << "hash2: " << hash.second << "\n";

hash=rollinghashkMer(hash.first,hash.second,newnuc,k);

std::cout << "hash1: " << hash.first << "\n";
std::cout << "hash2: " << hash.second << "\n";
}


//  Hashfunction for 8-mer
std::pair <unsigned,unsigned> hashkMer(const DnaString & kmer, const unsigned & k){
  unsigned hash=0;
  unsigned hash2=0;
  for (int i=0;i<k;++i){
    hash= hash << 2 | ordValue(kmer[i]);
  }
  for (int i=k-1;i>-1;--i){
    hash2= hash2 << 2 | (3-ordValue(kmer[i]));
  }
  return std::make_pair(hash,hash2);
}

// Rolling hashfunction for 8-mer
std::pair <unsigned,unsigned> rollinghashkMer(unsigned & oldHash, unsigned & oldHash2, const Dna & newnuc, const unsigned & k){
  oldHash=((oldHash << 2) | ordValue(newnuc)) % ((unsigned long)2 << (k*2-1));
  // oldHash2=(oldHash2 >> 2) |
  oldHash2=(3-ordValue(newnuc)) << (k*2-1);
  return std::make_pair(oldHash,oldHash2);
}
