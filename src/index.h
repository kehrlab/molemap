#ifndef INDEX_H
#define INDEX_H

# include <seqan/seq_io.h>
# include <seqan/sequence.h>
using namespace seqan;

struct openAddressingKmerHashtable{
  String<uint32_t> dir;
  String<uint32_t> pos;
  String<uint8_t> ref;
  String<int32_t> C;

  uint32_t k;
  int k_2;
  uint32_t m;
  uint32_t bucket_number;

  std::vector<std::string> lookChrom;
  std::vector<std::string> chromLength;

  openAddressingKmerHashtable(uint32_t k_in, uint32_t m_in ,uint32_t bucket_number_in){
    k=k_in;
    m=m_in;
    bucket_number=bucket_number_in;
    seqan::resize(dir,bucket_number+1,0);
    seqan::resize(C,bucket_number,-1);
    this->set_k_2();
  }

  openAddressingKmerHashtable(){}

  void clear(){
    seqan::clear(dir);
    seqan::clear(pos);
    seqan::clear(ref);
    seqan::clear(C);
    return;
  }

  void resize(uint32_t bucket_number_in){
    bucket_number=bucket_number_in;
    seqan::resize(dir,bucket_number+1,0);
    seqan::resize(C,bucket_number,-1);
    return;
  }

  void setBucketNumber(){
    bucket_number=length(C);
    return;
  }

  void set_k_2(){
    if (k>16){
      k_2=(k-15)*2;
    }else{
      k_2=0;
    }
    return;
  }

  uint32_t GetBkt(const int64_t & hash){
    uint64_t i=(uint64_t)hash%(uint64_t)bucket_number;
    int64_t d=0;
    while(C[i]!=(hash>>k_2) and C[i]!=-1){
      i=(i^(hash>>((d*16)%31)));
      i=(i+2*d+1)%(int64_t)bucket_number;
      d++;
    }
    return i;
  }

  uint32_t ReqBkt(const int64_t & hash){
    uint64_t i=(uint64_t)hash%(uint64_t)bucket_number;
    for(int d=0; d<1000; d++){
      #pragma omp atomic
      C[i]+=((hash>>k_2)+1)*(C[i]==-1)/*true if empty bucket found --> occupies bucket*/;
      if(C[i]==(hash>>k_2)){
        return i;
      } //return if correct bucket found
      i=(i^(hash>>((d*16)%31)));
      i=(i+2*d+1)%(int64_t)bucket_number;
    }
    std::cerr << "ERROR: no free bucket found after 1000 tries. Please increase bucket_number.\n";
    return i;
  }
};


int index(int argc, char const **argv);

#endif
