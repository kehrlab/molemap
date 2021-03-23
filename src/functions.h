#ifndef FUNCTIONS_H
#define FUNCTIONS_H

# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <iostream>
using namespace seqan;

void LoadBarcodeIndex(std::string & Index_name, std::vector<std::string> & BCI_barcodes, std::vector<std::pair<std::streampos,std::streampos>> & BCI_positions);
void ReturnBarcodeReads(std::vector<std::string> & BCI_barcodes, std::vector<std::pair<std::streampos,std::streampos>> & BCI_positions, std::string & barcode, const char* readfile1, const char* readfile2);
void ReportWindow(std::vector<std::tuple<double,unsigned char,unsigned,unsigned>> & best_windows, std::tuple<double,unsigned char,unsigned,unsigned> & candidate);
int64_t ReturnSmaller(const int64_t hash1,const int64_t hash2,const int64_t random_seed);
int64_t InitMini(const DnaString & string, const unsigned k, std::pair <int64_t, int64_t> & hash, const int64_t & maxhash,const int64_t random_seed, int64_t & minimizer_position);
void AppendPos(std::vector<std::tuple <unsigned char,unsigned,unsigned,unsigned>> & kmer_list, const int64_t & hash, const String<int64_t> & C,const String<unsigned> & dir,const String<std::pair <unsigned char,unsigned>> & pos, const unsigned bucket_number, unsigned & minimizer_active_bases);
DnaString hash2kmer(const int64_t & hash,const unsigned k);
std::pair <int64_t,int64_t> hashkMer(const DnaString & kmer, const unsigned k);
void rollinghashkMer(int64_t & oldHash, int64_t & oldHash2, const Dna5 & newnuc, const unsigned k, const int64_t & maxhash);
unsigned GetBkt(const int64_t & hash, const String<int64_t> & C, const unsigned bucket_number);
unsigned ReqBkt(const int64_t & hash, String<int64_t> & C, const unsigned bucket_number);
std::vector<std::pair <unsigned char,unsigned>> RetPos(const int64_t & hash, const String<int64_t> & C,const String<unsigned> & dir,const String<std::pair <unsigned char,unsigned>> & pos, const unsigned bucket_number);

#endif
