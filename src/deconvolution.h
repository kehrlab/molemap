#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H

# include <seqan/sequence.h>
# include <seqan/seq_io.h>
using namespace seqan;

int decon(int argc, char const **argv);

struct hash32_t{
  uint_fast8_t k;
  uint64_t fwd;
  uint64_t rev;
  uint64_t maxhash;

  hash32_t(const uint_fast8_t K){
    k=K;
    fwd=0;
    rev=0;
    maxhash=0;
    for (uint_fast8_t i=0;i<k;i++){
      maxhash= maxhash << 2 | 3;
    }
  }

  hash32_t(const DnaString & kmer){
    fwd=0;
    rev=0;
    k=length(kmer);
    maxhash=0;
    for (uint_fast8_t i=0;i<k;i++){
      maxhash= maxhash << 2 | 3;
    }
    for (uint_fast8_t i=0;i<k;++i){
      fwd= fwd << 2 | (int64_t)ordValue(kmer[i]);
      rev= rev << 2 | (int64_t)(3-ordValue(kmer[k-1-i]));
    }
  }

  void rehash(const DnaString & kmer){
    fwd=0;
    rev=0;
    for (uint_fast8_t i=0;i<length(kmer);++i){
      fwd= fwd << 2 | (int64_t)ordValue(kmer[i]);
      rev= rev << 2 | (int64_t)(3-ordValue(kmer[k-1-i]));
    }
    return;
  }

  void roll(const Dna5 & newnuc){
    fwd=((fwd << 2) | (int64_t)ordValue(newnuc)) & maxhash;
    rev=(rev >> 2) | (int64_t)(3-ordValue(newnuc)) << (k*2-2);
    return;
  }

  uint64_t bigger(){
    if(fwd > rev){
      return fwd;
    }
    return rev;
  }

  uint64_t hamming(uint64_t pos, uint64_t alt){ // alt is element of {1,2,3} pos is element of {0, ... , 31}
    uint64_t fwd_buff = fwd ^ (alt << (pos*2));
    uint64_t rev_buff = rev ^ (alt << (k*2-2-(pos*2)));

    if(fwd_buff > rev_buff){
      return fwd_buff;
    }
    return rev_buff;
  }

  void printKmer(){
    std::string kmer="";
    uint64_t bin;
    for(int i=k-1; i>=0; i--){
      bin = fwd>>(i*2) & 3;
      if(bin==0){kmer+="A";}
      if(bin==1){kmer+="C";}
      if(bin==2){kmer+="G";}
      if(bin==3){kmer+="T";}
    }
    std::cerr << kmer << "\t";
    kmer = "";
    for(int i=k-1; i>=0; i--){
      bin = rev>>(i*2) & 3;
      if(bin==0){kmer+="A";}
      if(bin==1){kmer+="C";}
      if(bin==2){kmer+="G";}
      if(bin==3){kmer+="T";}
    }
    std::cerr << kmer;
  }
};

struct coverageIndexBucket32_t{
  uint64_t id;
  uint32_t coverage;

  coverageIndexBucket32_t(){
    coverage = 0;
    id=0;
  }

  bool operator<(coverageIndexBucket32_t const& rhs){
      return coverage<rhs.coverage;
  }
};

struct coverageIndex32_t{
  std::vector<coverageIndexBucket32_t> table;
  uint64_t tableSize;

  coverageIndex32_t(){

  }

  void clear(){
    tableSize = 0;
    table.clear();
    return;
  }

  void resize(uint64_t tableSize_p){
    tableSize = tableSize_p;
    table.resize(tableSize);
    return;
  }

  bool getpos(uint64_t & hashId, uint64_t & pos){ // look up hashId in index. modifies pos and returns whether kmer is present
    pos = hashId % tableSize;
    for(int d=0; d<500; d++){
      // first try in loop
      if(table[pos].id == hashId){
        return true;
      } //return if correct bucket found
      if(table[pos].id == 0){ // if empty field hash is not in index
        return false;
      }
      pos = (pos^(hashId >> (d & 31))) % tableSize; //first shift in loop
      // second try in loop
      if(table[pos].id == hashId){
        return true;
      } //return if correct bucket found
      if(table[pos].id == 0){ // if empty field hash is not in index
        return false;
      }
      pos = (pos^(hashId << (d & 31))) % tableSize; // second shift in loop
    }
    std::cerr << "Warning: kmer index is overfilled.\n";
    return false;
  }

  bool insert(uint64_t & hashId){ // insert kmer into hashtable or increment coverage for kmer
    uint64_t pos;
    pos = hashId % tableSize;
    for(int d=0; d<500; d++){
      // first try in loop
      table[pos].id = (table[pos].id==0) ? hashId : table[pos].id; // occupies bucket if empty
      if(table[pos].id == hashId){
        table[pos].coverage++;
        return true;
      } //return if correct bucket found
      pos = (pos^(hashId >> (d & 31)) )% tableSize; // first shift in loop
      // second try in loop
      // #pragma omp atomic
      table[pos].id = (table[pos].id==0) ? hashId : table[pos].id; // occupies bucket if empty
      if(table[pos].id == hashId){
        // #pragma omp atomic
        table[pos].coverage++;
        return true;
      } //return if correct bucket found
      pos = (pos^(hashId << (d & 31))) % tableSize; // second shift in loop
    }
    std::cerr << "Warning: kmer index is overfilled.\n";
    return false;
  }

  bool lookup(uint64_t & hashId, uint32_t & coverage){
    uint64_t pos;
    if(this->getpos(hashId, pos)){
      coverage = table[pos].coverage;
      return true;
    }
    coverage = 0;
    return false;
  }

  uint32_t median(std::string & reference, uint32_t & k){
    // lookup reference kmer coverages and add to coverage list if coverage not 0
    std::vector<uint32_t> coverages; //vector of non zero coverages
    uint32_t coverage;
    coverages.reserve(reference.size());

    hash32_t hash(reference.substr(0,k));
    uint64_t hashId;

    for(std::string::iterator itr_ref = reference.begin()+k+1; itr_ref != reference.end(); itr_ref++){ // iterate over reference
      hashId=hash.bigger();
      if(this->lookup(hashId,coverage)){
        coverages.push_back(coverage);
      }
      hash.roll(*itr_ref);
    }

    // partial sort
    std::nth_element(coverages.begin(), coverages.begin()+coverages.size()/2, coverages.end());
    // get median
    uint32_t median = coverages[coverages.size()/2];
    // std::cerr << "\n\nMedian coverage: " << median << "\n\n";

    return median;
  }

};

#endif
