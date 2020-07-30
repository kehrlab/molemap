#ifndef FUNCTIONS_H
#define FUNCTIONS_H

# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <iostream>
using namespace seqan;

std::pair <unsigned,unsigned> hashkMer(const DnaString & kmer, const unsigned k);
std::pair <unsigned,unsigned> rollinghashkMer(unsigned & oldHash, unsigned & oldHash2, const Dna & newnuc, const unsigned k);
unsigned  GetBkt(const unsigned & hash, const String<unsigned> & C, const unsigned long long bucket_number);
unsigned  ReqBkt(const unsigned & hash, String<unsigned> & C, const unsigned long long bucket_number);
std::vector<std::pair <unsigned,unsigned>> RetPos(const DnaString & kmer, const String<unsigned> & C,const String<unsigned> & dir,const String<std::pair <unsigned,unsigned>> & pos, const unsigned long long bucket_number);

#endif
