#ifndef FUNCTIONS_H
#define FUNCTIONS_H

# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <iostream>
using namespace seqan;

std::pair <unsigned long long,unsigned long long> hashkMer(const DnaString & kmer, const unsigned k);
std::pair <unsigned long long,unsigned long long> rollinghashkMer(unsigned long long & oldHash, unsigned long long & oldHash2, const Dna & newnuc, const unsigned k);
unsigned long long GetBkt(const unsigned long long & hash, const String<unsigned> & C, const unsigned long long bucket_number);
unsigned  long long ReqBkt(const unsigned long long & hash, String<unsigned> & C, const unsigned long long bucket_number);
std::vector<std::pair <unsigned,unsigned>> RetPos(const unsigned long long & hash, const String<unsigned> & C,const String<unsigned> & dir,const String<std::pair <unsigned,unsigned>> & pos, const unsigned long long bucket_number);

#endif
