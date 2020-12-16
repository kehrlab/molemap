# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/index.h>
# include <map>
# include <iostream>
# include <time.h>
# include <fstream>
# include "./src/functions.h"
using namespace seqan;

/*
g++ countKmers.cpp -o countK
*/

/*
k can be 31 at max because of the hash function
*/

int main(int argc, char *argv[]){

  // auto begin = std::chrono::high_resolution_clock::now();

  // for (int a=0;a<100;a++){

  if(argc!=4){
    std::cerr << "Usage: ./countK InputFILE k Index_name \n\n";
    exit(-1);
  }

// reading the FastQ file

  StringSet<CharString> ids;
  StringSet<Dna5String> seqs;

  try {
    SeqFileIn file(argv[1]);

    readRecords(ids, seqs, file);

    close(file);
  }
  catch (ParseError const & e){
    std::cerr << "ERROR: input record is badly formatted. " << e.what() << std::endl;
  }
  catch (IOError const & e){
    std::cerr << "ERROR: input file can not be opened. " << e.what() << std::endl;
  }

  std::cerr << "Genome read in. \n";

  // defining key parameters

  const unsigned k=8;//std::stoi(argv[2]); // length of k-mer



  typedef Index< StringSet<Dna5String>, IndexQGram<UngappedShape<k>, OpenAddressing > > TIndex;
  TIndex index(seqs);

  std::cerr << "Index build. \n";

  //write index to file

  const char * IndName=argv[3];
  save(index,IndName);

  std::cerr << "Index writen to file.\n";

  // Kontrollausgabe

  hash(indexShape(index), "TTAAAAGC");
  for (unsigned i = 0; i < length(getOccurrences(index, indexShape(index))); ++i){
    std::cout << getSeqNo(getOccurrences(index, indexShape(index))[i]) << std::endl;
    std::cout << getSeqOffset(getOccurrences(index, indexShape(index))[i]) << std::endl;
  }



}
