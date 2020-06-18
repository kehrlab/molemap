# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <map>
# include <iostream>
# include <time.h>
# include <fstream>
# include <string>
using namespace seqan;

/*
g++ abund_from_index.cpp -o read_ind
*/


unsigned hashkMer(const DnaString & kmer, const unsigned & k);
unsigned  GetBkt(const unsigned & hash, const std::vector<unsigned> & C, const unsigned bucket_number);
std::vector<unsigned> RetPos(const DnaString & kmer, const std::vector<unsigned> & C,const std::vector<unsigned> & dir,const std::vector<unsigned> & pos, const unsigned bucket_number);





int main(int argc, char *argv[]){

  // auto begin = std::chrono::high_resolution_clock::now();

  // for (int a=0;a<100;a++){

  if(argc!=2){
    std::cerr << "Usage: ./read_ind InputFILE \n\n";
    exit(-1);
  }
  // defining key parameters

  unsigned k; // length of k-mer
  unsigned bucket_number; // should depend on k and the length of the indexed sequence
  unsigned pos_length;

// reading the index file

  std::ifstream index;
  index.open(argv[1]);
  index>>bucket_number>>k>>pos_length;

  std::vector<unsigned> dir(bucket_number,0);       // pow(4,k) depending on k-mer size
  std::vector<unsigned> pos(pos_length,0);         // length(seq)-k+1 runns into error
  std::vector<unsigned> C(pow(4,k)+1,0);
  std::vector<unsigned>::iterator itrv;
  std::vector<unsigned>::iterator itrv2;

  //reading C
  for (itrv=C.begin();itrv!=C.end();itrv++){
      index >> *itrv;
  }
  // reading dir
  for (itrv=dir.begin();itrv!=dir.end();itrv++){
      index >> *itrv;
  }
  // reading pos
  for (itrv=pos.begin();itrv!=pos.end();itrv++){
      index >> *itrv;
  }
  index.close();

  std::cout<<"k: "<<k << "\n";
  //Determining abundance of k-mers

  std::vector<unsigned> abundance(bucket_number-1,0);
  for (itrv=dir.begin(),itrv2=abundance.begin();itrv!=dir.end()-1;itrv++,itrv2++){
    *itrv2=*(itrv+1)-*itrv;
  }

  std::sort(abundance.rbegin(),abundance.rend());

  for (itrv=abundance.begin();itrv!=abundance.end();itrv++){
    std::cout<<*itrv<<" ";
  }

  // writing abundances to file


    std::ofstream abund;
    abund.open("abundances.txt");

    for (itrv=abundance.begin();itrv!=abundance.end();itrv++){
      abund<<*itrv<<" ";
    }

    abund.close();

// Test
  // std::cout << "\n\ninput " << k << "-mer: ";
  // std::string mer;
  // std::cin >> mer;
  // DnaString kmer=mer;
  // std::cout << "\n";


  // std::vector<unsigned> positions=RetPos(kmer,C,dir,pos,bucket_number);
  // for (itrv=positions.begin();itrv!=positions.end();itrv++){
  //   std::cout << *itrv <<" ";
  // }

  std::cout << "\n";

}

// return k-mer positions
std::vector<unsigned> RetPos(const DnaString & kmer, const std::vector<unsigned> & C,const std::vector<unsigned> & dir,const std::vector<unsigned> & pos, const unsigned bucket_number){
      std::vector<unsigned> positions;
      unsigned hash=hashkMer(kmer,length(kmer));
      int c=GetBkt(hash,C,bucket_number);
      for (unsigned i = dir[c];i!=dir[c+1];i++){
        positions.push_back(pos[i]);

      }
      return positions;
}

// Find correct Bucket
unsigned  GetBkt(const unsigned & hash, const std::vector<unsigned> & C, const unsigned bucket_number){
  std::srand(hash);
  unsigned i=std::rand()%bucket_number;
  unsigned d=0;
  unsigned counter=0;
  while(C[i]!=hash and C[i]!=-1){
    counter+=1;
    i=(i+2*d+1)%bucket_number;
    d++;
    if (counter > bucket_number){   // error if bucket_number not high enough
      std::cerr<<"\nERROR: Bucket number to small.\n";
      break;}
  }
  return i;
}

//  Hashfunction for 8-mer
unsigned hashkMer(const DnaString & kmer, const unsigned & k){
  unsigned hash=0;
  for (int i=0;i<k;++i){
    hash= hash << 2 | ordValue(kmer[i]);
  }
  return hash;
}
