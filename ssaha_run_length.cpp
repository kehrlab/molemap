# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <map>
# include <iostream>
# include <time.h>
# include <fstream>
# include "./src/functions.h"
using namespace seqan;

/*
g++ ssaha_run_length.cpp -o ssaha
*/



int main(int argc, char *argv[]){

if(argc!=3){
  std::cerr << "Usage: ./ssaha readFile k \n\n";
  exit(-1);
}


/*
loading in the reads
*/

StringSet<CharString> ids;
StringSet<Dna5String> k_mers;

// // reading from fastq/fasta files:
//
// try {
//   SeqFileIn file(argv[1]);
//
//   readRecords(ids, k_mers, file);
//
//   close(file);
// }
// catch (ParseError const & e){
//   std::cerr << "ERROR: input record is badly formatted. " << e.what() << std::endl;
// }
// catch (IOError const & e){
//   std::cerr << "ERROR: input file can not be opened. " << e.what() << std::endl;
// }

// reading in from BAM/SAM files

BamFileIn file;
if (!open(file, toCString(argv[1])))
{
    std::cerr << "ERROR: Could not open " << toCString(argv[1]) << std::endl;
    return 1;
}
BamHeader header;
readHeader(header,file);

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

std::vector<unsigned> all_run_lengths;
std::vector<unsigned> best_run_lengths;


// // reading from fastq/fasta files:
//
// typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
// for (TStringSetIterator it = begin(k_mers); it!=end(k_mers); ++it){ // iTERATING OVER THE reads
//   std::cerr << *it << "\n";

// reading in from BAM/SAM files

while (!atEnd(file))
{
  BamAlignmentRecord record;
  readRecord(record, file);
  Dna5String seq=record.seq;
  Dna5String * it=&seq;

/*
 Searching for the kmer
*/

    // building the Master_list
  std::vector<std::tuple<unsigned,int,unsigned>> Master_list;
  std::vector<std::tuple<unsigned,int,unsigned>>::iterator itrM;
  std::vector<std::pair<unsigned,unsigned>>::iterator itrp;

  if(int(length(*it)-k)>0){
    for (int t=0;t<(length(*it)-k);t++){
      std::vector<std::pair <unsigned,unsigned>> positions=RetPos(infix(*it,t,t+k), C, dir, pos, bucket_number);
      for (itrp=positions.begin();itrp!=positions.end();itrp++){
        Master_list.push_back(std::make_tuple((*itrp).first,(*itrp).second-t,(*itrp).second));
      }
    }
  }
  else {continue;}

    // sort the Master_list

    sort(Master_list.begin(),Master_list.end());

    // search for runs in the Master_list

    unsigned run=0;
    std::vector<unsigned> run_lengths;
    if (Master_list.empty()!=1){
      for (itrM=Master_list.begin();itrM!=(Master_list.end()-1);itrM++){
          if (std::get<0>(*itrM)==std::get<0>(*(itrM+1))){
            if (std::get<1>(*itrM)==std::get<1>(*(itrM+1))){
              run++;
            }
          }
          else {
            run_lengths.push_back(run);
            all_run_lengths.push_back(run);
            if (run>=2){
              std::cout << "run lenght: " << run << "\nposition: (" << std::get<0>(*itrM) << ") " << std::get<2>(*itrM)-run << " - "<< std::get<2>(*itrM)+k << "\n";
            }
            run=0;}

      }
    }

    if (run_lengths.empty()!=1){
      best_run_lengths.push_back(*max_element(run_lengths.begin(),run_lengths.end()));
    }

}

// saving run information to file

std::vector<unsigned>::iterator itrv;
std::ofstream all_runs;
all_runs.open("all_run_lengths.txt");
for (itrv=all_run_lengths.begin();itrv!=all_run_lengths.end();++itrv){
  all_runs << *itrv << " ";
}

std::ofstream best_runs;
best_runs.open("best_run_lengths.txt");
for (itrv=best_run_lengths.begin();itrv!=best_run_lengths.end();++itrv){
  best_runs << *itrv << " ";
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
