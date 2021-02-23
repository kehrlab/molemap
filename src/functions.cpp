# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <iostream>
#include "functions.h"

using namespace seqan;

// checks if candidate should be inserted into best_windows and inserts it at the correct palce
void report_window(std::vector<std::tuple<double,unsigned,unsigned,unsigned>> & best_windows, std::tuple<double,unsigned,unsigned,unsigned> & candidate){
  std::vector<std::tuple<double,unsigned,unsigned,unsigned>>::iterator itrbw;
  unsigned inserted=0;
  if (std::get<0>(candidate)>std::get<0>(best_windows.front())) {  // if current window better than worst window:
    for (itrbw=best_windows.begin();itrbw!=best_windows.end();itrbw++){
      if(std::get<0>(candidate) < std::get<0>(*itrbw)){                            // if (as soon as) quality is worse than quality in best_windows
        best_windows.insert(itrbw,candidate);                                      // insert new window there
        best_windows.erase(best_windows.begin());                                  // and delete worst window
        candidate=std::make_tuple(0,0,0,4294967295);
        return;
      }
    }
    best_windows.push_back(candidate);
    best_windows.erase(best_windows.begin());
    candidate=std::make_tuple(0,0,0,4294967295);
    return;
  }
}

// randomizes the hashvalues order
long long int ReturnSmaller(const long long int hash1,const long long int hash2,const long long int random_seed){
  if ((hash1^random_seed) < (hash2^random_seed)){
    return hash1^random_seed;
  } else {
    return hash2^random_seed;
  }
}

// initializes the minimizer
long long int InitMini(const DnaString & string, const unsigned k, std::pair <long long int, long long int> & hash, const long long int & maxhash,const long long int random_seed, long long int & minimizer_position){
  long long int minimizer=ReturnSmaller(hash.first,hash.second,random_seed);
  // long long int minimizer=std::min(hash.first,hash.second);
  // std::cerr << "minimizer: "<< minimizer << "min:" << std::min(hash.first,hash.second) << "\n";
  long long int minimizer_pos=0;
  for (unsigned i=1;i<length(string)-k+1;i++){
      rollinghashkMer(hash.first,hash.second,string[i+k-1],k,maxhash);
      if (minimizer > (hash.first^random_seed)){
      // if (std::min(minimizer,hash.first)!=minimizer){
        minimizer=hash.first;
        minimizer_pos=i;
      }
      if (minimizer > (hash.second^random_seed)){
      // if (std::min(minimizer,hash.second)!=minimizer){
        minimizer=hash.second;
        minimizer_pos=i;
      }
  }
  minimizer_position+=minimizer_pos;
  return minimizer;
}

// calculates following minimizer and reports if it replaces the old minimizer
// int RollMini(long long int & minimizer, std::pair <long long int, long long int> & hash, const Dna5 & newnuc, const unsigned k, const long long int & maxhash,const long long int random_seed,const unsigned long long bucket_number){
//   rollinghashkMer(hash.first,hash.second,newnuc,k,maxhash); // inline?!
//   if (minimizer > ReturnSmaller(hash.first,hash.second,random_seed)){
//     AppendPos(kmer_list, minimizer, C, dir, pos, bucket_number);
//     minimizer=ReturnSmaller(hash.first,hash.second,random_seed);
//     return 1;
//   }
//   return 0;
// }

//Insert k-mer positions into vector in sorted order
void AppendPos(std::vector<std::tuple <unsigned,unsigned,unsigned,unsigned>> & kmer_list, const long long int & hash, const String<int long long> & C,const String<unsigned long long> & dir,const String<std::pair <unsigned,unsigned>> & pos, const unsigned long long bucket_number,unsigned & minimizer_active_bases){
      // std::cerr <<"\nhash: " << hash << "\n";
      unsigned long long c=GetBkt(hash,C,bucket_number);
      // std::cerr << "c: " << c << "\n";
      unsigned long long abundance=dir[c+1]-dir[c];
      // std::cerr << "dir[c+1]: " << dir[c+1] << " dir[c]: " << dir[c] << "\n";
      // std::cerr << "abundance: " << abundance << "\n";
      kmer_list.reserve(kmer_list.size()+abundance);
      if (abundance<=10){
        for (unsigned long long i = dir[c];i!=dir[c+1];i++){
          kmer_list.push_back(std::make_tuple(pos[i].first,pos[i].second,abundance,minimizer_active_bases));
        }
      }
      return;
}

// return k-mer positions
std::vector<std::pair <unsigned,unsigned>> RetPos(const long long int & hash, const String<int long long> & C,const String<unsigned long long> & dir,const String<std::pair <unsigned,unsigned>> & pos, const unsigned long long bucket_number){
      std::vector<std::pair <unsigned,unsigned>> positions;
      unsigned long long c=GetBkt(hash,C,bucket_number);
      for (unsigned long long i = dir[c];i!=dir[c+1];i++){
        positions.push_back(pos[i]);
      }
      return positions;
}

// Find correct Bucket
unsigned long long  GetBkt(const long long int & hash, const String<int long long> & C, const unsigned long long bucket_number){
  // std::srand(hash);
  // unsigned long long i=std::rand()%bucket_number;
  long long int i=hash%(long long int)bucket_number;
  long long int d=0;
  unsigned counter=0;
  while(C[i]!=hash and C[i]!=-1){
    // std::cerr <<counter <<"\n";
    // std::cerr << "i before: " << i << "\n";
    i=(i^(hash>>((d*16)%31)));
    i=(i+2*d+1)%(long long int)bucket_number;
    // std::cerr << "i after: " << i << "\n";

    counter+=1;
    d++;
    if (counter > 100){   // error if bucket_number not high enough
      if (counter==101) {std::cerr<<"\nERROR: Bucket number to small.\n";
      std::cerr << " hash: " << hash << " seq: " << hash2kmer(hash,31) <<"\n";
      }
      // if (counter > 1000) {break;}
    }
  }
  // std::cerr << "tries: " << counter << "\n";
  return i;
}

// unsigned long long  GetBkt(const long long int & hash, const String<int long long> & C, const unsigned long long bucket_number){
//   // std::srand(hash);
//   // unsigned long long i=std::rand()%bucket_number;
//   long long int i=hash%(long long int)bucket_number;
//   long long int d=0;
//   unsigned counter=0;
//   while(C[i]!=hash and C[i]!=-1){
//
//     counter+=1;
//     i=(i+2*d+1)%(long long int)bucket_number;
//     d++;
//     if (counter > 100){   // error if bucket_number not high enough
//       if (counter==101) {std::cerr<<"\nERROR: Bucket number to small.\n";}
//       // if (counter > 1000) {break;}
//     }
//   }
//   return i;
// }


// Request a Bucket
unsigned long long ReqBkt(const long long int & hash, String<int long long> & C, const unsigned long long bucket_number){
  unsigned long long i = GetBkt(hash,C,bucket_number);
  C[i]=hash;
  return i;
}

// turn hashvalue into kmer
DnaString hash2kmer(const long long int & hash,const unsigned k){
  DnaString kmer="";
  long long int temp;
  for (size_t i = 0; i < k; i++) {
    temp=(hash>>((k-1-i)*2)) & (long long int)3;
    // std::cerr << "temp: " << temp <<"\n";
    if (temp==(long long int)0) {
      kmer+="A";
    }else if(temp==1){
      kmer+="C";
    }else if (temp==2){
      kmer+="G";
    }else{
      kmer+="T";
    }
  }
  return kmer;
}

//  Hashfunction for k-mer
std::pair <long long int, long long int> hashkMer(const DnaString & kmer, const unsigned k){
  long long int hash=0;
  long long int hash2=0;
  for (unsigned i=0;i<k;++i){
    hash= hash << 2 | (long long int)ordValue(kmer[i]);
    hash2= hash2 << 2 | (long long int)(3-ordValue(kmer[k-1-i]));
  }
  return std::make_pair(hash,hash2);
}

// Rolling hashfunction for k-mer
void rollinghashkMer(long long int & oldHash, long long int & oldHash2, const Dna5 & newnuc, const unsigned k, const long long int & maxhash){
  oldHash=((oldHash << 2) | (long long int)ordValue(newnuc)) & maxhash;
  oldHash2=(oldHash2 >> 2) | (long long int)(3-ordValue(newnuc)) << (k*2-2);
  return;
}
