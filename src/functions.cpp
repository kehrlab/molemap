# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <iostream>
# include <fstream>
# include "functions.h"

using namespace seqan;

bool compFunctionResult(result_t left, result_t right){
  if(left.chrom < right.chrom){return 1;}
  if(left.chrom > right.chrom){return 0;}
  if(left.start < right.start){return 1;}
  if(left.start > right.start){return 0;}
  if(left.end < right.end){return 1;}
  return 0;
}

void sortResults(std::vector<result_t> & results){
  std::sort(results.begin(), results.end(), compFunctionResult);
  return;
}


//retreive the barcode from 10x linked reads
std::string getBarcode(std::string id1, uint8_t barcode_length){
  std::size_t pos=id1.find(" ");
  std::string new_barcode;
  if(pos<1000){
    id1=id1.substr(id1.find(" "),10000);
    new_barcode=id1.substr(id1.find("BX:Z:")+5,barcode_length);
  }else{
    new_barcode="BAD_BARCODE_____";
  }
  return new_barcode;
}

CharString getPairedID(std::string id){
  id=id.substr(0, id.find_first_of(" \t\n"));
  id=id.substr(0, id.size()-2);
  return id;
}

CharString getID(std::string id){
  id=id.substr(0, id.find(" "));
  return id;
}

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
std::vector<std::pair<Dna5String,Dna5String>> ReturnBarcodeReads(std::vector<std::string> & BCI_barcodes, std::vector<std::pair<std::streampos,std::streampos>> & BCI_positions, std::string & barcode, const char* readfile1, const char* readfile2){
  SeqFileIn file1(readfile1);
  SeqFileIn file2(readfile2);
  Dna5String read1;
  Dna5String read2;
  CharString id;
  uint32_t pos = std::distance(BCI_barcodes.begin(), std::lower_bound(BCI_barcodes.begin(), BCI_barcodes.end(),barcode));
  file1.stream.file.seekg(std::get<0>(BCI_positions[pos]));
  file2.stream.file.seekg(std::get<1>(BCI_positions[pos]));
  std::streampos endpos=std::get<0>(BCI_positions[pos+1]);
  std::vector<std::pair<Dna5String,Dna5String>> reads;
  while(file1.stream.file.tellg()<endpos){
    readRecord(id,read1,file1);
    readRecord(id,read2,file2);
    reads.push_back(std::make_pair(read1,read2));
  }
  return reads;
}

// checks if candidate should be inserted into best_windows and inserts it at the correct palce
void ReportWindow(std::vector<std::tuple<double,uint8_t,uint32_t,uint32_t>> & best_windows, std::tuple<double,uint8_t,uint32_t,uint32_t> & candidate){
  std::vector<std::tuple<double,uint8_t,uint32_t,uint32_t>>::iterator itrbw;
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
    return hash1;
  } else {
    return hash2;
  }
}

bool IsSmaller(const int64_t hash1,const int64_t hash2,const int64_t random_seed){
  if ((hash1^random_seed) < (hash2^random_seed)){
    return true;
  } else {
    return false;
  }
}

// initializes the minimizer
int64_t InitMini(const DnaString & string, const uint32_t k, std::pair <int64_t, int64_t> & hash, const int64_t & maxhash,const int64_t random_seed, int64_t & minimizer_position){
  hash=hashkMer(infix(string,0,k),k);
  int64_t minimizer=ReturnSmaller(hash.first,hash.second,random_seed);
  int64_t minimizer_pos=0;
  for (uint32_t i=1;i<length(string)-k+1;i++){
      rollinghashkMer(hash.first,hash.second,string[i+k-1],k,maxhash);
      if (IsSmaller(hash.first,minimizer,random_seed)){
        minimizer=hash.first;
        minimizer_pos=i;
      }
      if (IsSmaller(hash.second,minimizer,random_seed)){
        minimizer=hash.second;
        minimizer_pos=i;
      }
  }
  minimizer_position+=minimizer_pos;
  return minimizer;
}

// Find correct Bucket
// uint32_t GetBkt(const int64_t & hash, const String<int32_t> & C, const uint32_t bucket_number, const int k_2){
//   uint64_t i=(uint64_t)hash%(uint64_t)bucket_number;
//   int64_t d=0;
//   while(C[i]!=(hash>>k_2) and C[i]!=-1){
//     i=(i^(hash>>((d*16)%31)));
//     i=(i+2*d+1)%(int64_t)bucket_number;
//     d++;
//   }
//   return i;
// }

// Request a Bucket multithread style
// uint32_t ReqBkt(const int64_t & hash, String<int32_t> & C, const uint32_t bucket_number, const int k_2){
//   uint64_t i=(uint64_t)hash%(uint64_t)bucket_number;
//   for(int d=0; d<1000; d++){
//     #pragma omp atomic
//     C[i]+=((hash>>k_2)+1)*(C[i]==-1)/*true if empty bucket found --> occupies bucket*/;
//     if(C[i]==(hash>>k_2)){
//       return i;
//     } //return if correct bucket found
//     i=(i^(hash>>((d*16)%31)));
//     i=(i+2*d+1)%(int64_t)bucket_number;
//   }
//   std::cerr << "ERROR: no free bucket found after 1000 tries. Please increase bucket_number.\n";
//   return i;
// }

// turn hashvalue into kmer
DnaString hash2kmer(const int64_t & hash,const uint8_t k){
  DnaString kmer="";
  int64_t temp;
  for (size_t i = 0; i < k; i++) {
    temp=(hash>>((k-1-i)*2)) & (int64_t)3;
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

bool NInKmer(Dna5String kmer, int64_t & position){
  for (uint8_t i=0;i<length(kmer);i++){
    if (kmer[i]=='N'){
      position+=i+1;
      return true;
    }
  }
  return false;
}

//  Hashfunction for k-mer
std::pair <int64_t, int64_t> hashkMer(const DnaString & kmer, const uint32_t k){
  int64_t hash=0;
  int64_t hash2=0;
  for (uint8_t i=0;i<k;++i){
    hash= hash << 2 | (int64_t)ordValue(kmer[i]);
    hash2= hash2 << 2 | (int64_t)(3-ordValue(kmer[k-1-i]));
  }
  return std::make_pair(hash,hash2);
}

// Rolling hashfunction for k-mer
void rollinghashkMer(int64_t & oldHash, int64_t & oldHash2, const Dna5 & newnuc, const uint8_t k, const int64_t & maxhash){
  oldHash=((oldHash << 2) | (int64_t)ordValue(newnuc)) & maxhash;
  oldHash2=(oldHash2 >> 2) | (int64_t)(3-ordValue(newnuc)) << (k*2-2);
  return;
}

int64_t getMaxHash(uint32_t & k){
  int64_t maxhash=0;
  for (uint8_t i=0;i<k;i++){
    maxhash= maxhash << 2 | 3;
  }
  return maxhash;
}

int64_t getRandSeed(uint32_t & k){
  std::srand(0);
  int64_t random_seed=0;
  for (uint8_t i=0;i<k;++i){
    random_seed= random_seed << 2 | (int64_t)(std::rand()%3);
  }
  return random_seed;
}
