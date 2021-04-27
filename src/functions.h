#ifndef FUNCTIONS_H
#define FUNCTIONS_H

# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <iostream>
using namespace seqan;

std::vector<Dna5String> GetReads(std::pair<std::streampos,std::streampos> & BCI_positions, std::streampos endpos, const char* readfile1, const char* readfile2);
void LoadBarcodeIndex(std::string & Index_name, std::vector<std::string> & BCI_barcodes, std::vector<std::pair<std::streampos,std::streampos>> & BCI_positions);
std::vector<std::pair<Dna5String,Dna5String>> ReturnBarcodeReads(std::vector<std::string> & BCI_barcodes, std::vector<std::pair<std::streampos,std::streampos>> & BCI_positions, std::string & barcode, const char* readfile1, const char* readfile2);
void ReportWindow(std::vector<std::tuple<double,uint_fast8_t,uint32_t,uint32_t>> & best_windows, std::tuple<double,uint_fast8_t,uint32_t,uint32_t> & candidate);
int64_t ReturnSmaller(const int64_t hash1,const int64_t hash2,const int64_t random_seed);
int64_t InitMini(const DnaString & string, const uint_fast8_t k, std::pair <int64_t, int64_t> & hash, const int64_t & maxhash,const int64_t random_seed, int64_t & minimizer_position);
void AppendPos(std::vector<std::tuple <uint_fast8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, const int64_t & hash, const String<int32_t> & C,const String<uint32_t> & dir, const String<uint_fast8_t> & ref, const String<uint32_t> & pos, const uint_fast32_t bucket_number, uint_fast8_t & minimizer_active_bases, const int k_2/*, pthread_mutex_t *lock*/);
DnaString hash2kmer(const int64_t & hash,const uint_fast8_t k);
std::pair <int64_t,int64_t> hashkMer(const DnaString & kmer, const uint_fast8_t k);
void rollinghashkMer(int64_t & oldHash, int64_t & oldHash2, const Dna5 & newnuc, const uint_fast8_t k, const int64_t & maxhash);
uint_fast32_t GetBkt(const int64_t & hash, const String<int32_t> & C, const uint_fast32_t bucket_number, const int k_2);
uint_fast32_t ReqBkt(const int64_t & hash, String<int32_t> & C, const uint_fast32_t bucket_number, const int k_2);

#endif
