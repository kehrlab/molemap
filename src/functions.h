#ifndef FUNCTIONS_H
#define FUNCTIONS_H

# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <iostream>
using namespace seqan;

void InsPos(std::vector<std::tuple <unsigned,unsigned,unsigned>> & kmer_list, const long long int & hash, const String<unsigned> & C,const String<unsigned> & dir,const String<std::pair <unsigned,unsigned>> & pos, const unsigned long long bucket_number);
std::pair <long long int,long long int> hashkMer(const DnaString & kmer, const unsigned k);
void rollinghashkMer(long long int & oldHash, long long int & oldHash2, const Dna & newnuc, const unsigned k, const long long int & maxhash);
unsigned long long GetBkt(const long long int & hash, const String<unsigned> & C, const unsigned long long bucket_number);
unsigned  long long ReqBkt(const long long int & hash, String<unsigned> & C, const unsigned long long bucket_number);
std::vector<std::pair <unsigned,unsigned>> RetPos(const long long int & hash, const String<unsigned> & C,const String<unsigned> & dir,const String<std::pair <unsigned,unsigned>> & pos, const unsigned long long bucket_number);

#endif
