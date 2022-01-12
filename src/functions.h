#ifndef FUNCTIONS_H
#define FUNCTIONS_H

# include "../seqan/seq_io.h"
# include "../seqan/sequence.h"
# include <iostream>
using namespace seqan;

std::string getBarcode(std::string id1, uint_fast8_t barcode_length);
CharString getID(std::string id);
std::vector<Dna5String> GetReads(std::pair<std::streampos,std::streampos> & BCI_positions, std::streampos endpos, const char* readfile1, const char* readfile2);
void LoadBarcodeIndex(std::string & Index_name, std::vector<std::string> & BCI_barcodes, std::vector<std::pair<std::streampos,std::streampos>> & BCI_positions);
std::vector<std::pair<Dna5String,Dna5String>> ReturnBarcodeReads(std::vector<std::string> & BCI_barcodes, std::vector<std::pair<std::streampos,std::streampos>> & BCI_positions, std::string & barcode, const char* readfile1, const char* readfile2);
void ReportWindow(std::vector<std::tuple<double,uint_fast8_t,uint32_t,uint32_t>> & best_windows, std::tuple<double,uint_fast8_t,uint32_t,uint32_t> & candidate);
int64_t ReturnSmaller(const int64_t hash1,const int64_t hash2,const int64_t random_seed);
bool IsSmaller(const int64_t hash1,const int64_t hash2,const int64_t random_seed);
int64_t InitMini(const DnaString & string, const uint_fast8_t k, std::pair <int64_t, int64_t> & hash, const int64_t & maxhash,const int64_t random_seed, int64_t & minimizer_position);
void AppendPos(std::vector<std::tuple <uint_fast8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, const int64_t & hash, const String<int32_t> & C,const String<uint32_t> & dir, const String<uint_fast8_t> & ref, const String<uint32_t> & pos, const uint_fast32_t bucket_number, uint_fast8_t & minimizer_active_bases, const int k_2);
DnaString hash2kmer(const int64_t & hash,const uint_fast8_t k);
std::pair <int64_t,int64_t> hashkMer(const DnaString & kmer, const uint_fast8_t k);
void rollinghashkMer(int64_t & oldHash, int64_t & oldHash2, const Dna5 & newnuc, const uint_fast8_t k, const int64_t & maxhash);
uint_fast32_t GetBkt(const int64_t & hash, const String<int32_t> & C, const uint_fast32_t bucket_number, const int k_2);
uint_fast32_t ReqBkt(const int64_t & hash, String<int32_t> & C, const uint_fast32_t bucket_number, const int k_2);
bool NInKmer(Dna5String kmer, int64_t & position);

class minimizer{
  public:
    uint64_t value;
    uint_fast8_t active_bases;
    int64_t position;
};

class minimizedSequence{
  public:
    uint_fast8_t k;
    uint_fast8_t m;
    Dna5String sequence;
    Dna5String kmer;
    int64_t last_position;
    int64_t position;
    std::pair<int64_t, int64_t> hash;
    int64_t random_seed;
    int64_t maxhash;
    minimizer new_minimizer;
    minimizer old_minimizer;
    bool at_end;

    minimizedSequence(Dna5String & sequenceC, uint_fast8_t & kC, uint_fast8_t & mC, int64_t & random_seedC, int64_t & maxhashC){ //constructor
      sequence=sequenceC;
      position=0;
      k=kC;
      m=mC;
      last_position=length(sequence)-m;
      random_seed=random_seedC;
      maxhash=maxhashC;
      kmer=infix(sequence,position,position+m);
      while(NInKmer(kmer,position)){
        kmer=infix(sequence,position,position+m);
      }
      hash=hashkMer(infix(sequence,position,position+k),k);
      new_minimizer.position=position;
      new_minimizer.value=InitMini(kmer, k, hash, maxhash, random_seed, new_minimizer.position);
      new_minimizer.active_bases=1;
      at_end=false;
      position++;
    }

    minimizer pop(){
      old_minimizer=new_minimizer;
      while(position <= last_position){
        if(position <= old_minimizer.position){ // if old minimizer in current window
          // iterate up the sequence, check if new minimizer
          if (sequence[position+m-1]!='N'){
            rollinghashkMer(hash.first,hash.second,sequence[position+m-1],k,maxhash);
          }else{
            position+=k+1;
            kmer=infix(sequence,position,position+k);
            while(NInKmer(kmer,position)){
              kmer=infix(sequence,position,position+k);
            }
            hash=hashkMer(kmer,k);
          }
          if(IsSmaller(old_minimizer.value,ReturnSmaller(hash.first,hash.second,random_seed),random_seed)){ // if minimizer stays active (new hash-values are bigger):
            old_minimizer.active_bases++;
          }else{ // create new_minimizer and return old_minimizer
            new_minimizer.value=ReturnSmaller(hash.first,hash.second,random_seed);
            new_minimizer.position=position+m-k;
            new_minimizer.active_bases=1;
            position++;
            return old_minimizer;
          }
        }else{
          // reinitialize new_minimizer, return old_minimizer
          kmer=infix(sequence,position,position+m);
          while(NInKmer(kmer,position)){
            kmer=infix(sequence,position,position+m);
          }
          new_minimizer.position=position;
          new_minimizer.value=InitMini(kmer, k, hash, maxhash, random_seed, new_minimizer.position);
          new_minimizer.active_bases=1;
          position++;
          return old_minimizer;
        }
        position++;
      }
      at_end=true;   // end of sequence reached: last minimizer returned
      return old_minimizer;
    }
};

#endif
