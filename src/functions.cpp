# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <iostream>
#include "functions.h"

using namespace seqan;

//Insert k-mer positions into vector in sorted order
void AppendPos(std::vector<std::tuple <unsigned,unsigned,unsigned>> & kmer_list, const long long int & hash, const String<int long long> & C,const String<unsigned long long> & dir,const String<std::pair <unsigned,unsigned>> & pos, const unsigned long long bucket_number){
      int c=GetBkt(hash,C,bucket_number);
      unsigned abundance=dir[c+1]-dir[c];
      kmer_list.reserve(kmer_list.size()+abundance);
      for (unsigned i = dir[c];i!=dir[c+1];i++){
        kmer_list.push_back(std::make_tuple(pos[i].first,pos[i].second,abundance));
      }
      return;
}
// return k-mer positions
std::vector<std::pair <unsigned,unsigned>> RetPos(const long long int & hash, const String<int long long> & C,const String<unsigned long long> & dir,const String<std::pair <unsigned,unsigned>> & pos, const unsigned long long bucket_number){
      std::vector<std::pair <unsigned,unsigned>> positions;
      int c=GetBkt(hash,C,bucket_number);
      for (unsigned i = dir[c];i!=dir[c+1];i++){
        positions.push_back(pos[i]);
      }
      return positions;
}

// Find correct Bucket
unsigned long long  GetBkt(const long long int & hash, const String<int long long> & C, const unsigned long long bucket_number){
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
unsigned  long long ReqBkt(const long long int & hash, String<int long long> & C, const unsigned long long bucket_number){
  unsigned long long i = GetBkt(hash,C,bucket_number);
  C[i]=hash;
  return i;
}

//  Hashfunction for k-mer
std::pair <long long int, long long int> hashkMer(const DnaString & kmer, const unsigned k){
  long long int hash=0;
  long long int hash2=0;
  for (unsigned i=0;i<k;++i){
    hash= hash << 2 | ordValue(kmer[i]);
    hash2= hash2 << 2 | (3-ordValue(kmer[k-1-i]));
  }
  return std::make_pair(hash,hash2);
}

// Rolling hashfunction for k-mer
void rollinghashkMer(long long int & oldHash, long long int & oldHash2, const Dna & newnuc, const unsigned k, const long long int & maxhash){
  oldHash=((oldHash << 2) | ordValue(newnuc)) & maxhash;
  oldHash2=(oldHash2 >> 2) | (3-ordValue(newnuc)) << (k*2-2);
  return;
}
