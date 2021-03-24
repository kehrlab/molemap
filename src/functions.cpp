# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <iostream>
# include <fstream>
#include "functions.h"

using namespace seqan;

// Loads BarcodeIndex from file into string
void LoadBarcodeIndex(std::string & Index_name, std::vector<std::string> & BCI_barcodes, std::vector<std::pair<std::streampos,std::streampos>> & BCI_positions){
  //construct file names
  std::string IndBC=Index_name;
  IndBC.append("_bc.txt");
  std::string IndPos=Index_name;
  IndPos.append("_pos.txt");
  // open barcode file
  std::ifstream file_bc;
  file_bc.open(IndBC, std::ios::binary);
  //load Barcodes line by line
  std::string barcode;
  while (getline(file_bc,barcode)){
    BCI_barcodes.push_back(barcode);
  }
  file_bc.close();
  //determine length of position table
  BCI_positions.reserve(BCI_barcodes.size()+1);
  std::streampos pos1;
  std::streampos pos2;
  //fill Position table
  std::ifstream file_pos;
  std::string posstr1;
  std::string posstr2;
  while(std::getline(file_pos,posstr1)){
    std::getline(file_pos,posstr2);
    pos1=stoi(posstr1);
    pos2=stoi(posstr2);
    BCI_positions.push_back(std::make_pair(pos1,pos2));
  }
  return;
}

// retreives all reads of a given barcode using the BarcodeIndex
void ReturnBarcodeReads(std::vector<std::string> & BCI_barcodes, std::vector<std::pair<std::streampos,std::streampos>> & BCI_positions, std::string & barcode, const char* readfile1, const char* readfile2){
  // SeqFileIn file1(readfile1);
  // SeqFileIn file2(readfile2);
  std::ifstream file1;
  file1.open(readfile1);
  std::ifstream file2;
  file2.open(readfile2);
  std::streampos posfile1;
  std::streampos posfile2;
  // Dna5String read1;
  // Dna5String read2;
  std::string read1;
  std::string read2;
  CharString id;
  uint_fast32_t pos = std::distance(BCI_barcodes.begin(), std::lower_bound(BCI_barcodes.begin(), BCI_barcodes.end(),barcode));
  // file1.stream.file.seekg(std::get<0>(BCI_positions[pos]));
  // file2.stream.file.seekg(std::get<1>(BCI_positions[pos]));
  file1.seekg(std::get<0>(BCI_positions[pos]));
  std::cerr << "position: " << file1.tellg() << "\n";
  std::getline(file1,read1);
  std::cerr << read1 << "\n";
  file2.seekg(std::get<1>(BCI_positions[pos]));
  std::streampos endpos=std::get<0>(BCI_positions[pos+1]);
  std::cerr << "\nstart: " << std::get<0>(BCI_positions[pos]) << "\tend: " << endpos << "\n";
  while(file1.tellg()<endpos){
    // std::cerr << "\n" << __LINE__;
    file1.ignore(100000,'\n');
    file2.ignore(100000,'\n');
    std::getline(file2,read2);
    file1.ignore(100000,'\n');
    file2.ignore(100000,'\n');
    file1.ignore(100000,'\n');
    file2.ignore(100000,'\n');
    // readRecord(id,read1,file1);
    // readRecord(id,read2,file2);
    // std::cerr << "\n" << __LINE__<<\n;
    // std::cerr << read1 << "\t" << read2 << "\n";
  }
  return;
}

// checks if candidate should be inserted into best_windows and inserts it at the correct palce
void ReportWindow(std::vector<std::tuple<double,uint_least8_t,uint32_t,uint32_t>> & best_windows, std::tuple<double,uint_least8_t,uint32_t,uint32_t> & candidate){
  std::vector<std::tuple<double,uint_least8_t,uint32_t,uint32_t>>::iterator itrbw;
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
int64_t ReturnSmaller(const int64_t hash1,const int64_t hash2,const int64_t random_seed){
  if ((hash1^random_seed) < (hash2^random_seed)){
    return hash1^random_seed;
  } else {
    return hash2^random_seed;
  }
}

// initializes the minimizer
int64_t InitMini(const DnaString & string, const uint_fast8_t k, std::pair <int64_t, int64_t> & hash, const int64_t & maxhash,const int64_t random_seed, int64_t & minimizer_position){
  int64_t minimizer=ReturnSmaller(hash.first,hash.second,random_seed);
  // int64_t minimizer=std::min(hash.first,hash.second);
  // std::cerr << "minimizer: "<< minimizer << "min:" << std::min(hash.first,hash.second) << "\n";
  int64_t minimizer_pos=0;
  for (uint_fast32_t i=1;i<length(string)-k+1;i++){
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

//Insert k-mer positions into vector in sorted order
void AppendPos(std::vector<std::tuple <uint_least8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, const int64_t & hash, const String<int64_t> & C,const String<uint32_t> & dir,const String<std::pair <uint_least8_t,uint32_t>> & pos, const uint_fast32_t bucket_number,uint_fast8_t & minimizer_active_bases){
      // std::cerr <<"\nhash: " << hash << "\n";
      uint_fast32_t c=GetBkt(hash,C,bucket_number);
      // std::cerr << "c: " << c << "\n";
      uint_fast32_t abundance=dir[c+1]-dir[c];
      // std::cerr << "dir[c+1]: " << dir[c+1] << " dir[c]: " << dir[c] << "\n";
      // std::cerr << "abundance: " << abundance << "\n";
      kmer_list.reserve(kmer_list.size()+abundance);
      if (abundance<=10){
        for (uint_fast32_t i = dir[c];i!=dir[c+1];i++){
          kmer_list.push_back(std::make_tuple(pos[i].first,pos[i].second,abundance,minimizer_active_bases));
        }
      }
      return;
}

// return k-mer positions
std::vector<std::pair <uint_least8_t,uint32_t>> RetPos(const int64_t & hash, const String<int64_t> & C,const String<uint32_t> & dir,const String<std::pair <uint_least8_t,uint32_t>> & pos, const uint_fast32_t bucket_number){
      std::vector<std::pair <uint_least8_t,uint32_t>> positions;
      uint_fast32_t c=GetBkt(hash,C,bucket_number);
      for (uint_fast32_t i = dir[c];i!=dir[c+1];i++){
        positions.push_back(pos[i]);
      }
      return positions;
}

// Find correct Bucket
uint_fast32_t GetBkt(const int64_t & hash, const String<int64_t> & C, const uint_fast32_t bucket_number){
  // std::srand(hash);
  // uint64_t i=std::rand()%bucket_number;
  int64_t i=hash%(int64_t)bucket_number;
  int64_t d=0;
  // unsigned counter=0;
  while(C[i]!=hash and C[i]!=-1){
    // std::cerr <<counter <<"\n";
    // std::cerr << "i before: " << i << "\n";
    i=(i^(hash>>((d*16)%31)));
    i=(i+2*d+1)%(int64_t)bucket_number;
    // std::cerr << "i after: " << i << "\n";

    // counter+=1;
    d++;
    // if (counter > 100){   // error if bucket_number not high enough
    //   if (counter==101) {std::cerr<<"\nERROR: Bucket number to small.\n";
    //   std::cerr << " hash: " << hash << " seq: " << hash2kmer(hash,31) <<"\n";
    //   }
    //   // if (counter > 1000) {break;}
    // }
  }
  // std::cerr << d+1 << " ";
  // std::cerr << "tries: " << counter << "\n";
  return i;
}

// Request a Bucket
uint_fast32_t ReqBkt(const int64_t & hash, String<int64_t> & C, const uint_fast32_t bucket_number){
  uint_fast32_t i = GetBkt(hash,C,bucket_number);
  C[i]=hash;
  return i;
}

// turn hashvalue into kmer
DnaString hash2kmer(const int64_t & hash,const uint_fast8_t k){
  DnaString kmer="";
  int64_t temp;
  for (size_t i = 0; i < k; i++) {
    temp=(hash>>((k-1-i)*2)) & (int64_t)3;
    // std::cerr << "temp: " << temp <<"\n";
    if (temp==(int64_t)0) {
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
std::pair <int64_t, int64_t> hashkMer(const DnaString & kmer, const uint_fast8_t k){
  int64_t hash=0;
  int64_t hash2=0;
  for (uint_fast8_t i=0;i<k;++i){
    hash= hash << 2 | (int64_t)ordValue(kmer[i]);
    hash2= hash2 << 2 | (int64_t)(3-ordValue(kmer[k-1-i]));
  }
  return std::make_pair(hash,hash2);
}

// Rolling hashfunction for k-mer
void rollinghashkMer(int64_t & oldHash, int64_t & oldHash2, const Dna5 & newnuc, const uint_fast8_t k, const int64_t & maxhash){
  oldHash=((oldHash << 2) | (int64_t)ordValue(newnuc)) & maxhash;
  oldHash2=(oldHash2 >> 2) | (int64_t)(3-ordValue(newnuc)) << (k*2-2);
  return;
}
