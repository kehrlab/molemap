# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <map>
# include <iostream>
# include <time.h>
# include <fstream>
using namespace seqan;

/*
g++ ssaha_run_length.cpp -o ssaha
*/

std::pair <unsigned,unsigned> hashkMer(const DnaString & kmer, const unsigned & k);
std::vector<unsigned> RetPos(const DnaString & kmer, const std::vector<unsigned> & C,const std::vector<unsigned> & dir,const std::vector<unsigned> & pos, const unsigned bucket_number);
unsigned  GetBkt(const unsigned & hash, const std::vector<unsigned> & C, const unsigned bucket_number);




int main(int argc, char *argv[]){

if(argc!=3){
  std::cerr << "Usage: ./ssaha IndexFILE kmerFile \n\n";
  exit(-1);
}

unsigned k; // length of k-mers in index
unsigned bucket_number;
unsigned pos_length;


// reading the k-mers

StringSet<CharString> ids;

StringSet<Dna5String> k_mers;

try {
  SeqFileIn file(argv[2]);

  readRecords(ids, k_mers, file);

  close(file);
}
catch (ParseError const & e){
  std::cerr << "ERROR: input record is badly formatted. " << e.what() << std::endl;
}
catch (IOError const & e){
  std::cerr << "ERROR: input file can not be opened. " << e.what() << std::endl;
}

// reading the Index

std::ifstream index;
index.open(argv[1]);
index>>bucket_number>>k>>pos_length;

std::vector<unsigned> dir(bucket_number,0);       // pow(4,k) depending on k-mer size
std::vector<unsigned> pos(pos_length,0);         // length(seq)-k+1 runns into error
std::vector<unsigned> C(pow(4,k)+1,0);
std::vector<unsigned>::iterator itrv;

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

std::vector<unsigned> all_run_lengths;
std::vector<unsigned> best_run_lengths;

typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
for (TStringSetIterator it = begin(k_mers); it!=end(k_mers); ++it){


  // Searching for the kmer

    // building the Master_list

  std::vector<std::tuple<unsigned,unsigned>> Master_list;
  std::vector<std::tuple<unsigned,unsigned>>::iterator itrM; // better to declare this in the for loops?


  for (int t=0;t<(length(*it)-k);t++){
    std::vector<unsigned> positions=RetPos(infix(*it,t,t+k), C, dir, pos, bucket_number);

    for (itrv=positions.begin();itrv!=positions.end();itrv++){
      Master_list.push_back(std::make_tuple(*itrv-t,*itrv));
    }
  }

    // sort the Master_list

    sort(Master_list.begin(),Master_list.end());

    // search for runs in the Master_list

    unsigned run=0;

    std::vector<unsigned> run_lengths;
    for (itrM=Master_list.begin();itrM!=(Master_list.end()-1);itrM++){
        std::cerr << "itrM:" << std::get<0>(*itrM) << "\n";
        std::cerr << "(itrM+1):" << std::get<0>(*(itrM+1))<< "\n";
        if (std::get<0>(*itrM)==std::get<0>(*(itrM+1))){
          run++;
        }
        else {
          run_lengths.push_back(run);
          all_run_lengths.push_back(run);
          if (run>=2){
            std::cout << "run lenght: " << run << "\nposition: " << std::get<0>(*itrM) << " - "<< std::get<1>(*itrM)+k << "\n";
          }
          run=0;}

    }


    best_run_lengths.push_back(*max_element(run_lengths.begin(),run_lengths.end()));


}

// saving run information to file

std::ofstream all_runs;
all_runs.open("all_run_lengths.txt");
for (itrv=all_run_lengths.begin();itrv!=all_run_lengths.end();++itrv){
  all_runs << *itrv;
}

std::ofstream best_runs;
best_runs.open("best_run_lengths.txt");
for (itrv=best_run_lengths.begin();itrv!=best_run_lengths.end();++itrv){
  best_runs << *itrv;
}

// Kontrollausgabe

//
// std::cout << "k: " << k << "\n";
// std::cout << "kmer: " << kmer << "\n\n";
//
// for (itrM=Master_list.begin();itrM!=Master_list.end();itrM++){
//   std::cout << std::get<0>(*itrM) << ", "  << std::get<1>(*itrM) << "\n";
// }
// std::cout << "\n";


}


// Functions

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

std::pair <unsigned,unsigned> hashkMer(const DnaString & kmer, const unsigned & k){
  unsigned hash=0;
  unsigned hash2=0;
  for (int i=0;i<k;++i){
    hash= hash << 2 | ordValue(kmer[i]);
    hash2= hash2 << 2 | (3-ordValue(kmer[k-1-i]));
  }
  return std::make_pair(hash,hash2);
}

std::vector<unsigned> RetPos(const DnaString & kmer, const std::vector<unsigned> & C,const std::vector<unsigned> & dir,const std::vector<unsigned> & pos, const unsigned bucket_number){
      std::vector<unsigned> positions;
      std::pair <unsigned,unsigned> hash=hashkMer(kmer,length(kmer));
      int c=GetBkt(std::min(hash.first,hash.second),C,bucket_number);
      for (unsigned i = dir[c];i!=dir[c+1];i++){
        positions.push_back(pos[i]);
      }
      return positions;
}
