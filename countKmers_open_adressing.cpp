# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <map>
# include <iostream>
using namespace seqan;

/*
g++ countKmers_open_adressing.cpp -o countK
*/
unsigned hashEightMer(const DnaString & kmer);
unsigned rollinghashEightMer(unsigned & oldHash, const Dna & newnuc);
unsigned  GetBkt(const unsigned & hash, const std::vector<unsigned> & C);
unsigned  ReqBkt(const unsigned & hash, std::vector<unsigned> & C);
std::vector<unsigned> RetPos(const DnaString & kmer, const std::vector<unsigned> & C,const std::vector<unsigned> & dir,const std::vector<unsigned> & pos);

int main(int argc, char *argv[]){


  for (int a=0;a<100;a++){

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
    std::cerr << "ERROR: input file can not be opened. " << e.what() << std::endl;
  }

  // defining key parameters

  unsigned k=8; // length of k-mer
  unsigned bucket_number=30000;

  // concatination of all sequences

  DnaString seq=concat(seqs);

  // building index storage

  std::vector<unsigned> dir(bucket_number,0);       // pow(4,k) depending on k-mer size
  std::vector<unsigned> pos(length(seq),0);         // length(seq)-k+1 runns into error
  std::vector<unsigned> C(pow(4,k)+1,-1);
  std::vector<unsigned>::iterator itrv;
  std::vector<unsigned>::reverse_iterator itrvr;

  // counting k-mers

  unsigned hash=hashEightMer(infix(seq,0,k));
  unsigned c;

  for (unsigned i = 0;i<length(seq)-k+1;++i){
    c=ReqBkt(hash,C);
    dir[c+1]+=1;
    hash=rollinghashEightMer(hash,seq[i+k]);
  }
  unsigned sum=length(seq)-k+1;

  // cumulative sum

  for (itrvr=dir.rbegin();itrvr!=dir.rend();++itrvr){
    sum-=*itrvr;
    *itrvr=sum;
  }

  // filling pos

  hash=hashEightMer(infix(seq,0,k));
  for (unsigned i = 0;i<length(seq)-k+1;++i){
    c=GetBkt(hash,C);
    pos[dir[c+1]]=i;
    dir[c+1]++;
    hash=rollinghashEightMer(hash,seq[i+k]);
  }

  // Kontrollausgabe


  // DnaString testDNA="TTAGTAAA";
  // std::vector<unsigned> positions=RetPos(testDNA,C,dir,pos);
  // for (itrv=positions.begin();itrv!=positions.end();itrv++){
    // std::cout << *itrv <<" ";
  // }
  // std::cout << "\n";
}
}

// return k-mer positions
std::vector<unsigned> RetPos(const DnaString & kmer, const std::vector<unsigned> & C,const std::vector<unsigned> & dir,const std::vector<unsigned> & pos){
      std::vector<unsigned> positions;
      unsigned hash=hashEightMer(kmer);
      int c=GetBkt(hash,C);
      for (unsigned i = dir[c];i!=dir[c+1];i++){
        positions.push_back(pos[i]);
      }
      return positions;
}

// Find correct Bucket
unsigned  GetBkt(const unsigned & hash, const std::vector<unsigned> & C){
  unsigned bucket_number=30000;
  std::srand(hash);
  unsigned i=std::rand()%bucket_number;
  unsigned d=0;
  while(C[i]!=hash and C[i]!=-1){
    i=(i+2*d+1)%bucket_number;
    d++;

  }
  return i;
}

// Request a Bucket
unsigned  ReqBkt(const unsigned & hash, std::vector<unsigned> & C){
  unsigned i = GetBkt(hash,C);
  C[i]=hash;
  return i;
}

//  Hashfunction for 8-mer
unsigned hashEightMer(const DnaString & kmer){
  unsigned hash=0;
  for (int i=0;i<8;++i){
    hash= hash << 2 | ordValue(kmer[i]);
  }
  return hash;
}

// Rolling hashfunction for 8-mer
unsigned rollinghashEightMer(unsigned & oldHash, const Dna & newnuc){
  oldHash=((oldHash << 2) | ordValue(newnuc)) % ((unsigned long)2 << (15));
  return oldHash;
}
