# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <map>
# include <iostream>

using namespace seqan;

/*
g++ countKmers_rolling_hash.cpp -o countK
*/
unsigned hashEightMer(const DnaString & kmer);
unsigned rollinghashEightMer(const unsigned & oldHash, const Dna & oldnuc, const Dna & newnuc);

int main(int argc, char *argv[]){

  if(argc!=2){
    std::cerr << "Usage: ./countK InputFILE \n\n";
    exit(-1);
  }

// reading the FastQ file

  StringSet<CharString> ids;
  StringSet<DnaString> seqs;

  try {
    SeqFileIn file(argv[1]);

    readRecords(ids, seqs, file);

    close(file);
  }
  catch (ParseError const & e){
    std::cerr << "ERROR: input record is badly formatted. " << e.what() << std::endl;
  }
  catch (IOError const & e){
    std::cerr << "ERROR: file can not be opened. " << e.what() << std::endl;
  }

// concatination of all sequences

DnaString seq;
seq=concat(seqs);

// map erzeugen. durch sequenzen iterieren. und k mere zählen;
  // DnaString seq=seqs[0];

  std::map<unsigned,std::vector<unsigned>> kmerPositions;
  std::map<unsigned,std::vector<unsigned>>::iterator kit;


  unsigned key=hashEightMer(infix(seq,0,8));

  for (unsigned i = 0;i<length(seq)-7;++i){
      kit=kmerPositions.find(key);
      if (kit==kmerPositions.end()){
        kmerPositions.insert({key,{i}});
      }
      else {
        kit->second.push_back(i);
      }

      key=rollinghashEightMer(key,seq[i],seq[i+8]);

  }

// Kontrollausgebe

  DnaString test_kmer="GGGAAAAA";

  key=hashEightMer(test_kmer);
  std::cout << "hash value: " << key <<"\n";
  std::vector<unsigned>::iterator itrv;
  std::cout << "positions:";
  kit=kmerPositions.find(key);
  if (kit!=kmerPositions.end())
    for (itrv=kit->second.begin();itrv!=kit->second.end();++itrv){
      std::cout << " " << *itrv;
    }
  else{
    std::cout << "DNA string not in Dataset";
    }
}



//  Hashfunction for 8-mer
unsigned hashEightMer(const DnaString & kmer){
  unsigned hash=0;
  for (int i=0;i<8;++i){
    hash+=ordValue(kmer[i])*pow(4,i);
  }
  return hash;
}
// Rolling hashfunction for 8-mer
unsigned rollinghashEightMer(const unsigned & oldHash, const Dna & oldnuc, const Dna & newnuc){
  unsigned hash;

  hash=(oldHash-ordValue(oldnuc))/4+ordValue(newnuc)*16384; //16384=pow(4,7)

  return hash;
}
