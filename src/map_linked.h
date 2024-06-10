#ifndef MAPLINKED_H
#define MAPLINKED_H

# include <seqan/seq_io.h>
# include <seqan/sequence.h>
using namespace seqan;

bool SearchID(SeqFileIn & file, CharString id, std::streampos startpos, std::streampos endpos, CharString (*getIdFunction)(std::string));
int maplinked(int argc, char const **argv);

#endif
