# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <map>
# include <iostream>

using namespace seqan;

/*
g++ countKmers_index.cpp -o countK
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

  // defining k

  unsigned k=8;

  // concatination of all sequences

  DnaString seq;
  seq=concat(seqs);

  // building index storage

  std::vector<unsigned> dir(pow(4,k)+1,0);       // pow(4,k) depending on k-mer size
  std::vector<unsigned> pos(length(seq)-k+1,0);
  // std::vector<unsigned>::iterator itrv;
  std::vector<unsigned>::reverse_iterator itrvr;

  // counting k-mers

  unsigned c=hashEightMer(infix(seq,0,k));

  for (unsigned i = 0;i<length(seq)-k+1;++i){
    dir[c+1]+=1;
    c=rollinghashEightMer(c,seq[i],seq[i+k]);
  }
  unsigned sum=length(seq)-k+1;

  // cumulative sum

  for (itrvr=dir.rbegin();itrvr!=dir.rend();++itrvr){
    sum-=*itrvr;
    *itrvr=sum;
  }

  // filling pos

  c=hashEightMer(infix(seq,0,k));

  for (unsigned i = 0;i<length(seq)-k+1;++i){
    pos[dir[c+1]]=i;
    dir[c+1]++;
    c=rollinghashEightMer(c,seq[i],seq[i+k]);
  }

  // Kontrollausgabe

  DnaString testDNA="TTAGTAAA";
  std::cout << hashEightMer(testDNA) << "\n";
  std::cout << dir[hashEightMer(testDNA)+1] << "\n";
  for (int i=dir[hashEightMer(testDNA)];i!=dir[hashEightMer(testDNA)+1];i++){
    std::cout << pos[i] <<" ";
  }
  std::cout << "\n";

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
