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

int window_size=100000;   // size of the genomic windows to wich the reads are matched
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

/*
reading the Index
*/

String<unsigned> dir;
String<std::pair <unsigned,unsigned>> pos;
String<unsigned> C;


String<std::pair <unsigned,unsigned>, External<ExternalConfigLarge<>> > extpos;
if (!open(extpos, "index_pos_21.txt", OPEN_RDONLY)){
  throw std::runtime_error("Could not open index counts file." );
}
assign(pos, extpos, Exact());
close(extpos);

String<unsigned, External<> > extdir;
if (!open(extdir, "index_dir_21.txt", OPEN_RDONLY)){
  throw std::runtime_error("Could not open index counts file." );
}
assign(dir, extdir, Exact());
close(extdir);

String<unsigned, External<> > extC;
if (!open(extC, "index_C_21.txt", OPEN_RDONLY)){
  throw std::runtime_error("Could not open index counts file." );
}
assign(C, extC, Exact());
close(extC);

unsigned k=std::stoi(argv[2]); // length of k-mers in index
unsigned long long bucket_number=length(C);


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
  // find k-mers and append positions to kmer_list

  if(int(length(*it)-k)>0){
    for (int t=0;t<(length(*it)-k);t++){
      std::vector<std::pair <unsigned,unsigned>> positions=RetPos(infix(*it,t,t+k), C, dir, pos, bucket_number);
      for (itrp=positions.begin();itrp!=positions.end();itrp++){
        kmer_list.push_back(std::make_tuple((*itrp).first,(*itrp).second,positions.size()));
      }
    }
  }
  else {continue;}
}

std::cerr << "k-mers listed.\n";
// std::cerr<<"kmer_list is build. \n";
//sorting k-mers by position in reference

sort(kmer_list.begin(),kmer_list.end());

// iterate over k-mers (sliding window)

std::cerr << "\nbefore log building.";
float lookLog[100]= {0,1,0.693147,1.09861,1.38629,1.60944,1.79176,1.94591,2.07944,2.19722,2.30259,2.3979,2.48491,2.56495,2.63906,2.70805,2.77259,2.83321,2.89037,2.94444,2.99573,3.04452,3.09104,3.13549,3.17805,3.21888,3.2581,3.29584,3.3322,3.3673,3.4012,3.43399,3.46574,3.49651,3.52636,3.55535,3.58352,3.61092,3.63759,3.66356,3.68888,3.71357,3.73767,3.7612,3.78419,3.80666,3.82864,3.85015,3.8712,3.89182,3.91202,3.93183,3.95124,3.97029,3.98898,4.00733,4.02535,4.04305,4.06044,4.07754,4.09434,4.11087,4.12713,4.14313,4.15888,4.17439,4.18965,4.20469,4.21951,4.23411,4.2485,4.26268,4.27667,4.29046,4.30407,4.31749,4.33073,4.34381,4.35671,4.36945,4.38203,4.39445,4.40672,4.41884,4.43082,4.44265,4.45435,4.46591,4.47734,4.48864,4.49981,4.51086,4.52179,4.5326,4.54329,4.55388,4.56435,4.57471,4.58497,4.59512};
std::cerr << "\nafter log building.";
std::cerr << "\n array size: " << sizeof(lookLog)/sizeof(lookLog[0]);
#define REF(X) std::get<0>(*(X))
#define POS(X) std::get<1>(*(X))
#define ABU(X) std::get<2>(*(X))

unsigned slider=1;
double window_quality=0;
if (ABU(kmer_list.begin())==1){                                 // updating window quality
  window_quality+=2;}
else if(ABU(kmer_list.begin())>99){
  window_quality+=0.2;}
else{
  window_quality+=1/lookLog[ABU(kmer_list.begin())];}

std::cerr << "\nafter first log usage.\n";

std::vector<std::tuple<double,unsigned,unsigned>> best_windows(window_count,std::make_tuple(0,0,0)); //(maping_quality, reference, position in referende)
std::vector<std::tuple<double,unsigned,unsigned>>::iterator itrbw;

// std::cerr<<"iteration prepared. \n";

// for(itrk=kmer_list.begin();itrk!=kmer_list.end(); itrk++){
//   std::cout<<"\nref: " << std::get<0>(*itrk) << "\tpos: " << std::get<1>(*itrk)<< "\tabu: "<< std::get<2>(*itrk);
// }


for(itrk=kmer_list.begin()+1;itrk!=kmer_list.end();itrk++){ // iterating over kmer_list
    // trimm the begining of the window
    if (ABU(itrk-1)==1){                                 // updating window quality
      window_quality-=2;
    }else if(ABU(itrk-1)>99){
      window_quality-=0.2;
    }else{
      window_quality-=1/lookLog[ABU(itrk-1)];
    }
    // expanding window to maximum length
    while(REF(itrk)==REF(itrk+slider) && POS(itrk+slider)-POS(itrk)<=window_size){ // while k-mers inside sliding window
        if (ABU(itrk+slider)==1){                                 // updating window quality
          window_quality+=2;}
        else if(ABU(itrk+slider)>99){window_quality+=0.2;}
        else {window_quality+=1/lookLog[ABU(itrk+slider)];}
        slider++;
    }
    slider--;

    // checking if current window qualifies
    int inserted=0;
    if( window_quality > std::get<0>(best_windows.front()) ) { // if current window better than worst window:
       for (itrbw=best_windows.begin();itrbw!=best_windows.end();itrbw++){                             // iterate over best_windows
         if (std::get<1>(*itrbw)==REF(itrk) && abs((int)POS(itrk)-(int)std::get<2>(*itrbw))<=window_size){ // if overlapping window: keep better window and break loop.
           if (window_quality > std::get<0>(*itrbw)){
             *itrbw=std::make_tuple(window_quality, REF(itrk), POS(itrk));
             sort(best_windows.begin(),best_windows.end());
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


// std::cerr<<"best_windows found. \n";

// Konttrollausgabe


for(itrbw=best_windows.begin();itrbw!=best_windows.end(); itrbw++){
  std::cout<<"\nquality: " << std::get<0>(*itrbw) << "\tref: " << std::get<1>(*itrbw)<< "\tpos: "<< std::get<2>(*itrbw);
}

}
