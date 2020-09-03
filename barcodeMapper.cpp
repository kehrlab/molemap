# include <iostream>
# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <fstream>
# include "./src/functions.h"
using namespace seqan;

/*
g++ BarcodeMapper.cpp -o bcmap
*/

float lookLog(unsigned num);

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

#define REF(X) std::get<0>(*(X))
#define POS(X) std::get<1>(*(X))
#define ABU(X) std::get<2>(*(X))

unsigned slider=1;
double window_quality=0;
if (ABU(kmer_list.begin())!=1){
  window_quality+=1/lookLog(ABU(kmer_list.begin()));
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
      window_quality-=1/lookLog(ABU(itrk-1));
    } else {window_quality-=2;}
    // expanding window to maximum length
    while(REF(itrk)==REF(itrk+slider) && POS(itrk+slider)-POS(itrk)<=window_size){ // while k-mers inside sliding window
        if (ABU(itrk+slider)!=1){                                 // updating window quality
          window_quality+=1/lookLog(ABU(itrk+slider));
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


float lookLog(unsigned num){
  switch(num){
    case 2: return 0.693147;
    case 3: return 1.09861;
    case 4: return 1.38629;
    case 5: return 1.60944;
    case 6: return 1.79176;
    case 7: return 1.94591;
    case 8: return 2.07944;
    case 9: return 2.19722;
    case 10: return 2.30259;
    case 11: return 2.3979;
    case 12: return 2.48491;
    case 13: return 2.56495;
    case 14: return 2.63906;
    case 15: return 2.70805;
    case 16: return 2.77259;
    case 17: return 2.83321;
    case 18: return 2.89037;
    case 19: return 2.94444;
    case 20: return 2.99573;
    case 21: return 3.04452;
    case 22: return 3.09104;
    case 23: return 3.13549;
    case 24: return 3.17805;
    case 25: return 3.21888;
    case 26: return 3.2581;
    case 27: return 3.29584;
    case 28: return 3.3322;
    case 29: return 3.3673;
    case 30: return 3.4012;
    case 31: return 3.43399;
    case 32: return 3.46574;
    case 33: return 3.49651;
    case 34: return 3.52636;
    case 35: return 3.55535;
    case 36: return 3.58352;
    case 37: return 3.61092;
    case 38: return 3.63759;
    case 39: return 3.66356;
    case 40: return 3.68888;
    case 41: return 3.71357;
    case 42: return 3.73767;
    case 43: return 3.7612;
    case 44: return 3.78419;
    case 45: return 3.80666;
    case 46: return 3.82864;
    case 47: return 3.85015;
    case 48: return 3.8712;
    case 49: return 3.89182;
    case 50: return 3.91202;
    case 51: return 3.93183;
    case 52: return 3.95124;
    case 53: return 3.97029;
    case 54: return 3.98898;
    case 55: return 4.00733;
    case 56: return 4.02535;
    case 57: return 4.04305;
    case 58: return 4.06044;
    case 59: return 4.07754;
    case 60: return 4.09434;
    case 61: return 4.11087;
    case 62: return 4.12713;
    case 63: return 4.14313;
    case 64: return 4.15888;
    case 65: return 4.17439;
    case 66: return 4.18965;
    case 67: return 4.20469;
    case 68: return 4.21951;
    case 69: return 4.23411;
    case 70: return 4.2485;
    case 71: return 4.26268;
    case 72: return 4.27667;
    case 73: return 4.29046;
    case 74: return 4.30407;
    case 75: return 4.31749;
    case 76: return 4.33073;
    case 77: return 4.34381;
    case 78: return 4.35671;
    case 79: return 4.36945;
    case 80: return 4.38203;
    case 81: return 4.39445;
    case 82: return 4.40672;
    case 83: return 4.41884;
    case 84: return 4.43082;
    case 85: return 4.44265;
    case 86: return 4.45435;
    case 87: return 4.46591;
    case 88: return 4.47734;
    case 89: return 4.48864;
    case 90: return 4.49981;
    case 91: return 4.51086;
    case 92: return 4.52179;
    case 93: return 4.5326;
    case 94: return 4.54329;
    case 95: return 4.55388;
    case 96: return 4.56435;
    case 97: return 4.57471;
    case 98: return 4.58497;
    case 99: return 4.59512;
    default: return 5;
  }
  std::cerr << "no case worked!";
  return 5;
}
