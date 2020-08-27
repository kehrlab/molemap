# include <iostream>
# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <fstream>
# include "./src/functions.h"
using namespace seqan;

/*
g++ BarcodeMapper.cpp -o bcmap
*/

int main(int argc, char *argv[]){

if(argc!=3){
  std::cerr << "Usage: ./bcmap readFile k \n\n";
  exit(-1);
}

/*
defining Parameters
*/

int window_size=200;   // size of the genomic windows to wich the reads are matched
int window_count=1000;   // amount of saved candidate windows

/*
loading in the reads
*/

StringSet<CharString> ids;
StringSet<Dna5String> reads;

// reading from fastq/fasta files:

try {
  SeqFileIn file(argv[1]);

  readRecords(ids, reads, file);

  close(file);
}
catch (ParseError const & e){
  std::cerr << "ERROR: input record is badly formatted. " << e.what() << std::endl;
}
catch (IOError const & e){
  std::cerr << "ERROR: input file can not be opened. " << e.what() << std::endl;
}

/*
reading the Index
*/

String<unsigned> dir;
String<std::pair <unsigned,unsigned>> pos;
String<unsigned> C;


String<std::pair <unsigned,unsigned>, External<ExternalConfigLarge<>> > extpos;
if (!open(extpos, "index_pos.txt", OPEN_RDONLY)){
  throw std::runtime_error("Could not open index counts file." );
}
assign(pos, extpos, Exact());
close(extpos);

String<unsigned, External<> > extdir;
if (!open(extdir, "index_dir.txt", OPEN_RDONLY)){
  throw std::runtime_error("Could not open index counts file." );
}
assign(dir, extdir, Exact());
close(extdir);

String<unsigned, External<> > extC;
if (!open(extC, "index_C.txt", OPEN_RDONLY)){
  throw std::runtime_error("Could not open index counts file." );
}
assign(C, extC, Exact());
close(extC);

unsigned k=std::stoi(argv[2]); // length of k-mers in index
unsigned bucket_number=length(C);


/*
Searching for all kmers of reads with the same Barcode
*/

// building the kmer_list for a specific Barcode (maybe exclude very frequent k-mers?)
std::vector<std::tuple<unsigned,unsigned,unsigned>> kmer_list;   // (i,j,a)   i=reference (Chromosome), j=position of matching k-mer in reference, a=abundance of k-mer in reference
std::vector<std::tuple<unsigned,unsigned,unsigned>>::iterator itrk;
std::vector<std::pair<unsigned,unsigned>>::iterator itrp;

std::cerr << "Index and reads loaded.\n";

typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
for (TStringSetIterator it = begin(reads); it!=end(reads); ++it){ // Iterating over the reads
  std::cerr << 0;
  // find k-mers and append positions to kmer_list

  if(int(length(*it)-k)>0){
    std::cerr <<1;
    for (int t=0;t<(length(*it)-k);t++){
      std::cerr <<2;
      std::vector<std::pair <unsigned,unsigned>> positions=RetPos(infix(*it,t,t+k), C, dir, pos, bucket_number);
      std::cerr <<3;
      for (itrp=positions.begin();itrp!=positions.end();itrp++){
        std::cerr <<4;
        kmer_list.push_back(std::make_tuple((*itrp).first,(*itrp).second,positions.size()));
        std::cerr << "5\n";
      }
    }
  }
  else {continue;}
}

std::cerr << "k-mers listed.\n"
// std::cerr<<"kmer_list is build. \n";
//sorting k-mers by position in reference

sort(kmer_list.begin(),kmer_list.end());

// iterate over k-mers (sliding window)

#define REF(X) std::get<0>(*(X))
#define POS(X) std::get<1>(*(X))
#define ABU(X) std::get<2>(*(X))

unsigned slider=1;
double window_quality=0;
if (ABU(kmer_list.begin())!=1){
  window_quality+=1/log(ABU(kmer_list.begin()));
} else {window_quality+=2;}


std::vector<std::tuple<double,unsigned,unsigned>> best_windows(window_count,std::make_tuple(0,0,0)); //(maping_quality, reference, position in referende)
std::vector<std::tuple<double,unsigned,unsigned>>::iterator itrbw;

// std::cerr<<"iteration prepared. \n";

// for(itrk=kmer_list.begin();itrk!=kmer_list.end(); itrk++){
//   std::cout<<"\nref: " << std::get<0>(*itrk) << "\tpos: " << std::get<1>(*itrk)<< "\tabu: "<< std::get<2>(*itrk);
// }


for(itrk=kmer_list.begin()+1;itrk!=kmer_list.end();itrk++){ // iterating over kmer_list
    // trimm the begining of the window
    if (ABU(itrk-1)!=1){
      window_quality-=1/log(ABU(itrk-1));
    } else {window_quality-=2;}
    // expanding window to maximum length
    while(REF(itrk)==REF(itrk+slider) && POS(itrk+slider)-POS(itrk)<=window_size){ // while k-mers inside sliding window
        if (ABU(itrk+slider)!=1){                                 // updating window quality
          window_quality+=1/log(ABU(itrk+slider));
        } else {window_quality+=2;}
        slider++;
    }
    slider--;

    // checking if current window qualifies
    int inserted=0;
    if( window_quality > std::get<0>(best_windows.front()) ) { // if current window better than worst window:
      for (itrbw=best_windows.begin();itrbw!=best_windows.end();itrbw++){                             // iterate over best_windows
        if (std::get<1>(*itrbw)==REF(itrk) && abs((int)POS(itrk)-(int)std::get<2>(*itrbw))<=window_size){ // if overlapping window: keep better window and break loop.
          if (window_quality > std::get<0>(*itrbw)){
            best_windows.erase(itrbw);
          }
          else{inserted=1;}
          break;
        }
      }
      if(inserted==0){                                                                                     // if no overlapping window was found:
        for (itrbw=best_windows.begin()+1;itrbw!=best_windows.end();itrbw++){                             // iterate over best_windows
          if(window_quality < std::get<0>(*itrbw)){                                                       // if (as soon as) quality is worse than quality in best_windows
              best_windows.insert(itrbw,std::make_tuple(window_quality, REF(itrk), POS(itrk)));      // insert new window there
              best_windows.erase(best_windows.begin());                                                // and delete worst window
              inserted=1;
              break;
          }
        }
        if(inserted==0){
          best_windows.push_back(std::make_tuple(window_quality, REF(itrk), POS(itrk)));      // if no better window in best_windows insert new window at the end
          best_windows.erase(best_windows.begin());
        }
      }
    }
}


// std::cerr<<"best_windows found. \n";

// Konttrollausgabe


for(itrbw=best_windows.begin();itrbw!=best_windows.end(); itrbw++){
  std::cout<<"\nquality: " << std::get<0>(*itrbw) << "\tref: " << std::get<1>(*itrbw)<< "\tpos: "<< std::get<2>(*itrbw);
}

}
