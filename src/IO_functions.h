#ifndef IO_FUNCTIONS_H
#define IO_FUNCTIONS_H

# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include "parser.h"
# include "index.h"
using namespace seqan;

int loadReference(StringSet<Dna5String> & seqs, indexOptions & options);
void readKmerIndex(openAddressingKmerHashtable & Index, std::string & kmer_index_name);
int checkReadfile(std::string & readfileName);
// BamFileOut writeSamHeader(longmapOptions & options, openAddressingKmerHashtable & Index, std::fstream & output, int argc, char const ** argv);
void writeTmpHeader(longmapOptions & options, openAddressingKmerHashtable & Index, int argc, char const ** argv);
void writeSamHeader(longmapOptions & options);
void writeBamHeader(longmapOptions & options);
// void writeSamOutput(std::vector<BamAlignmentRecord> & records, std::vector<std::string> & lookChrom,  longmapOptions & options);
void writeSamOutput(std::vector<BamAlignmentRecord> & records, longmapOptions & options);
void writeBamOutput(std::vector<BamAlignmentRecord> & records, longmapOptions & options);
void writeSamCout(std::vector<BamAlignmentRecord> & records, longmapOptions & options);
void writeBamCout(std::vector<BamAlignmentRecord> & records, longmapOptions & options);
void writeHistogram(std::vector<uint32_t> & histogram, std::string & output_file);

#endif
