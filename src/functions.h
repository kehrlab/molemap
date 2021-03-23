#ifndef FUNCTIONS_H
#define FUNCTIONS_H

# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <iostream>
using namespace seqan;

void LoadBarcodeIndex(std::string & Index_name, std::vector<std::string> & BCI_barcodes, std::vector<std::pair<std::streampos,std::streampos>> & BCI_positions);
void ReturnBarcodeReads(std::vector<std::string> & BCI_barcodes, std::vector<std::pair<std::streampos,std::streampos>> & BCI_positions, std::string & barcode, const char* readfile1, const char* readfile2);
void ReportWindow(std::vector<std::tuple<double,unsigned char,unsigned,unsigned>> & best_windows, std::tuple<double,unsigned char,unsigned,unsigned> & candidate);
long long int ReturnSmaller(const long long int hash1,const long long int hash2,const long long int random_seed);
long long int InitMini(const DnaString & string, const unsigned k, std::pair <long long int, long long int> & hash, const long long int & maxhash,const long long int random_seed, long long int & minimizer_position);
void AppendPos(std::vector<std::tuple <unsigned char,unsigned,unsigned,unsigned>> & kmer_list, const long long int & hash, const String<int long long> & C,const String<unsigned> & dir,const String<std::pair <unsigned char,unsigned>> & pos, const unsigned bucket_number, unsigned & minimizer_active_bases);
DnaString hash2kmer(const long long int & hash,const unsigned k);
std::pair <long long int,long long int> hashkMer(const DnaString & kmer, const unsigned k);
void rollinghashkMer(long long int & oldHash, long long int & oldHash2, const Dna5 & newnuc, const unsigned k, const long long int & maxhash);
unsigned GetBkt(const long long int & hash, const String<int long long> & C, const unsigned bucket_number);
unsigned ReqBkt(const long long int & hash, String<int long long> & C, const unsigned bucket_number);
std::vector<std::pair <unsigned char,unsigned>> RetPos(const long long int & hash, const String<int long long> & C,const String<unsigned> & dir,const String<std::pair <unsigned char,unsigned>> & pos, const unsigned bucket_number);

#endif
