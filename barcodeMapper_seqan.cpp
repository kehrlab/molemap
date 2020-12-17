# include <iostream>
# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/index.h>
# include <fstream>
# include "./src/functions.h"
# include <time.h>
using namespace seqan;

/*
g++ BarcodeMapper.cpp -o bcmap
*/


int main(int argc, char *argv[]){

if(argc!=4){
  std::cerr << "Usage: ./bcmap readFile k Index_name\n\n";
  exit(-1);
}

/*
defining Parameters
*/

int window_size=5000;   // size of the genomic windows to wich the reads are matched
int window_count=100;   // amount of saved candidate windows

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

std::cerr << "reads loaded.\n";

/*
reading the Index
*/

typedef Index< StringSet<Dna5String>, IndexQGram<UngappedShape<8>, OpenAddressing > > TIndex;
TIndex Index;
const char * IndName=argv[3];
open(Index, IndName);
// Shape< StringSet<Dna5String>, IndexQGram<UngappedShape<8>, OpenAddressing > > myShape;
/*
Searching for all kmers of reads with the same Barcode
*/

// building the kmer_list for a specific Barcode (maybe exclude very frequent k-mers?)
std::vector<std::tuple<unsigned,unsigned,unsigned>> kmer_list;   // (i,j,a)   i=reference (Chromosome), j=position of matching k-mer in reference, a=abundance of k-mer in reference
std::vector<std::tuple<unsigned,unsigned,unsigned>>::const_iterator itrk;
std::vector<std::pair<unsigned,unsigned>>::const_iterator itrp;

std::cerr << "Index and reads loaded.\n";

auto tbegin = std::chrono::high_resolution_clock::now();

typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
std::cerr << __LINE__ << " ";
for (TStringSetIterator it = begin(reads); it!=end(reads); ++it){ // Iterating over the reads
  std::cerr << __LINE__ << " ";
  hashInit(indexShape(Index), begin(*it));
  std::cerr << __LINE__ << " ";
  for (unsigned t = 0; t < length(*it) - length(indexShape(Index)) + 1; ++t){
    std::cerr << __LINE__ << " ";
    hashNext(indexShape(Index), begin(*it) + t);
    std::cerr << __LINE__ << " ";
    unsigned abundance=length(getOccurrences(Index, indexShape(Index)));
    std::cerr << __LINE__ << " ";
    for (unsigned i = 0; i < abundance; ++i){
      std::cerr << __LINE__ << " ";
      kmer_list.push_back(std::make_tuple((unsigned)getSeqNo(getOccurrences(Index, indexShape(Index))[i]),(unsigned)getSeqOffset(getOccurrences(Index, indexShape(Index))[i]),abundance));
      std::cerr << __LINE__ << " ";
    }
  }
}
std::cerr << "k-mers listed.  \n";


auto tend = std::chrono::high_resolution_clock::now();
std::cout << "\ntime: "<<(float)std::chrono::duration_cast<std::chrono::milliseconds>(tend-tbegin).count()/1000 << " s\n";// << "ns" << std::endl;
tbegin = std::chrono::high_resolution_clock::now();
//sorting k-mers by position in reference

sort(kmer_list.begin(),kmer_list.end());

tend = std::chrono::high_resolution_clock::now();
std::cout << "\nsorting time: "<<(float)std::chrono::duration_cast<std::chrono::milliseconds>(tend-tbegin).count()/1000 << " s\n";// << "ns" << std::endl;
tbegin = std::chrono::high_resolution_clock::now();

// iterate over k-mers (sliding window)

float lookQual[100]= {0,1024,6.24989, 0.624853, 0.195309, 0.0926038, 0.0541504, 0.0358415, 0.0257197, 0.0195267, 0.0154498, 0.0126139, 0.0105548, 0.00900754, 0.00781189, 0.0068662, 0.00610341, 0.00547777, 0.00495714, 0.00451843, 0.00414462, 0.003823, 0.00354385, 0.00329967, 0.00308456, 0.00289387, 0.00272383, 0.00257141, 0.00243412, 0.0023099, 0.00219705, 0.00209414, 0.00199997, 0.0019135, 0.00183386, 0.00176031, 0.0016922, 0.00162897, 0.00157012, 0.00151524, 0.00146395, 0.00141593, 0.00137087, 0.00132852, 0.00128865, 0.00125106, 0.00121556, 0.00118199, 0.00115019, 0.00112005, 0.00109142, 0.00106421, 0.00103832, 0.00101365, 0.000990122, 0.00096766, 0.000946195, 0.000925665, 0.00090601, 0.000887177, 0.000869117, 0.000851784, 0.000835136, 0.000819134, 0.000803742, 0.000788926, 0.000774656, 0.000760902, 0.000747638, 0.000734837, 0.000722477, 0.000710537, 0.000698994, 0.00068783, 0.000677027, 0.000666568, 0.000656437, 0.000646619, 0.0006371, 0.000627866, 0.000618906, 0.000610208, 0.00060176, 0.000593551, 0.000585573, 0.000577815, 0.000570269, 0.000562926, 0.000555778, 0.000548817, 0.000542037, 0.000535431, 0.000528992, 0.000522713, 0.000516589, 0.000510615, 0.000504785, 0.000499093, 0.000493536, 0.000488108};

#define REF(X) std::get<0>(*(X))
#define POS(X) std::get<1>(*(X))
#define ABU(X) std::get<2>(*(X))

// for(itrk=kmer_list.begin();itrk!=kmer_list.end();itrk++){
//   std::cerr << "(" << REF(itrk) <<"," << POS(itrk) <<","<<ABU(itrk)<< ")" << " ";
// }

unsigned slider=1;
double window_quality=0;
if(ABU(kmer_list.begin())>99){
  window_quality+=0.00032;
}else{
  window_quality+=lookQual[ABU(kmer_list.begin())]; // lookQual = 1/(log(abund)^5)
}


std::vector<std::tuple<double,unsigned,int>> best_windows(window_count,std::make_tuple(0,0,-(window_size+10))); //(maping_quality, reference, position in referende)
std::vector<std::tuple<double,unsigned,int>>::iterator itrbw;
// std::cerr<<"iteration prepared. \n";

// for(itrk=kmer_list.begin();itrk!=kmer_list.end(); itrk++){
//   std::cout<<"\nref: " << std::get<0>(*itrk) << "\tpos: " << std::get<1>(*itrk)<< "\tabu: "<< std::get<2>(*itrk);
// }


for(itrk=kmer_list.begin()+1;itrk!=kmer_list.end();itrk++){ // iterating over kmer_list
    // trimm the begining of the window
    if(ABU(itrk-1)>99){
      window_quality-=0.00032;
    }else{
      window_quality-=lookQual[ABU(itrk-1)];
    }
    // expanding window to maximum length
    while(REF(itrk)==REF(itrk+slider) && POS(itrk+slider)-POS(itrk)<=window_size){ // while k-mers inside sliding window
        if(ABU(itrk+slider)>99){
          window_quality+=0.00032;
        }else {
          window_quality+=lookQual[ABU(itrk+slider)];
        }
        slider++;
    }
    slider--;

    // checking if current window qualifies
    int inserted=0;

    if( window_quality > std::get<0>(best_windows.front()) && POS(itrk)!=POS(itrk-1)) { // if current window better than worst window:
       for (itrbw=best_windows.begin();itrbw!=best_windows.end();itrbw++){                             // iterate over best_windows
         if (std::get<1>(*itrbw)==REF(itrk) && abs((int)POS(itrk)-(int)std::get<2>(*itrbw))<=window_size){ // if overlapping window: keep better window and break loop.
           if (window_quality > std::get<0>(*itrbw)){
             best_windows.erase(itrbw);
             for (itrbw=best_windows.begin()+1;itrbw!=best_windows.end();itrbw++){                             // iterate over best_windows
                if(window_quality < std::get<0>(*itrbw)){                                                       // if (as soon as) quality is worse than quality in best_windows
                    best_windows.insert(itrbw,std::make_tuple(window_quality, REF(itrk), POS(itrk)));      // insert new window there
                    inserted=1;
                    break;
                }
              }
              if(inserted==0){
                best_windows.push_back(std::make_tuple(window_quality, REF(itrk), POS(itrk)));      // if no better window in best_windows insert new window at the end
              }


           }
           inserted=1;
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

// trimm unused parts of best_windows

while(std::get<0>(*best_windows.begin())==0){
  best_windows.erase(best_windows.begin());
}
std::cerr<<"best_windows found. ";
tend = std::chrono::high_resolution_clock::now();
std::cout <<"\ntime: "<< (float)std::chrono::duration_cast<std::chrono::milliseconds>(tend-tbegin).count()/1000 << " s\n";// << "ns" << std::endl;

// Konttrollausgabe


for(itrbw=best_windows.begin();itrbw!=best_windows.end(); itrbw++){
  std::cout<<"\nquality: " << std::get<0>(*itrbw) << "\tref: " << std::get<1>(*itrbw)<< "\tpos: "<< std::get<2>(*itrbw);
}

}
