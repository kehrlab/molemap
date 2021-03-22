# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <iostream>
# include <fstream>
# include "./src/functions.h"
# include <time.h>
using namespace seqan;

/*
g++ BarcodeMapper.cpp -o bcmap
*/
void MapKmerList(std::vector<std::tuple<char,unsigned,unsigned,unsigned>> & kmer_list, unsigned & max_window_size, unsigned & max_gap_size, unsigned & window_count, const char* file);

struct bcmapOptions{
  std::string readfile1;
  std::string readfile2;
  std::string index_name;
  std::string bci_name;
  unsigned k;
  unsigned mini_window_size;

  bcmapOptions() :
  k(31), mini_window_size(35)
  {}
  };

seqan::ArgumentParser::ParseResult parseCommandLine(bcmapOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bcmap");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "readfile1"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "readfile2"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "index_name"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "bci_"));

    // Define Options
    addOption(parser, seqan::ArgParseOption(
        "k", "kmer_length", "length of k-mers in kmer indexed",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "k", "31");
    addOption(parser, seqan::ArgParseOption(
        "m", "mini_window_size", "length of minimizing window",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "m", "35");


    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;}

    // Extract option values.
    getOptionValue(options.k, parser, "k");
    getOptionValue(options.mini_window_size, parser, "m");

    getArgumentValue(options.readfile1, parser, 0);
    getArgumentValue(options.readfile2, parser, 1);
    getArgumentValue(options.index_name, parser, 2);
    getArgumentValue(options.bci_name, parser, 3);

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv){

if(argc < 5){
  std::cerr << "Usage: ./bcmap readFile1 readFile2 Index_name BCI_name [-flags] \n\n";
  exit(-1);
}
// parsing command line arguments
bcmapOptions options;
seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

// If parsing was not successful then exit with code 1 if there were errors.
// Otherwise, exit with code 0 (e.g. help was printed).
if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;

std::cout <<'\n'
          << "k          \t" << options.k << '\n'
          << "m          \t" << options.mini_window_size << '\n'
          << "readfile1  \t" << options.readfile1 << '\n'
          << "readfile2  \t" << options.readfile2 << '\n'
          << "index_name \t" << options.index_name << '\n'
          << "bci_name   \t" << options.bci_name << '\n\n';

unsigned k = options.k;
unsigned mini_window_size = options.mini_window_size;

/*
defining Parameters
*/

unsigned max_window_size=200000;  //5000;   // maximum size of the genomic windows to wich the reads are matched
unsigned max_gap_size=20000;     // maximum gap size between two adjacent k_mer hits
unsigned window_count=100;   // amount of saved candidate windows

const char* resultfile="bc_windows.txt";



/*
reading the Index
*/

std::cerr << "Reading in the k-mer index ";
auto tbegin = std::chrono::high_resolution_clock::now();

String<unsigned> dir;
String<std::pair <char,unsigned>> pos;
String<int long long> C;
//
std::string IndPos=options.index_name;
IndPos.append("_pos.txt");
std::string IndDir=options.index_name;
IndDir.append("_dir.txt");
std::string IndC=options.index_name;
IndC.append("_C.txt");

String<std::pair <char,unsigned>, External<ExternalConfigLarge<>> > extpos;
if (!open(extpos, IndPos.c_str(), OPEN_RDONLY)){
  throw std::runtime_error("Could not open index counts file." );
}
assign(pos, extpos, Exact());
close(extpos);
std::cerr <<".";

String<unsigned, External<> > extdir;
if (!open(extdir, IndDir.c_str(), OPEN_RDONLY)){
  throw std::runtime_error("Could not open index counts file." );
}
assign(dir, extdir, Exact());
close(extdir);
std::cerr <<".";


String<int long long, External<> > extC;
if (!open(extC, IndC.c_str(), OPEN_RDONLY)){
  throw std::runtime_error("Could not open index counts file." );
}
assign(C, extC, Exact());
close(extC);


long long int maxhash;
for (int i=0;i<k;i++){
  maxhash= maxhash << 2 | 3;
}

std::srand(0);
long long int random_seed=0;
for (unsigned i=0;i<k;++i){
  random_seed= random_seed << 2 | (long long int)(std::rand()%3);
}

unsigned bucket_number=length(C);

std::cerr <<". done!";
std::cerr << " in: " << (float)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-tbegin).count()/1000 << "s\n";
/*
loading in the reads
*/


// loading reads from fastq/fasta files:


try {         // opening read-files
  SeqFileIn file1(toCString(options.readfile1));
  SeqFileIn file2(toCString(options.readfile2));
  close(file1);
  close(file2);
}
catch (ParseError const & e){
  std::cerr << "ERROR: input record is badly formatted. " << e.what() << std::endl;
}
catch (IOError const & e){
  std::cerr << "ERROR: input file can not be opened. " << e.what() << std::endl;
}

std::cerr << "read file checked.\n";

/*
Searching for all kmers of reads with the same Barcode
*/

// building the kmer_list for a specific Barcode
std::vector<std::tuple<char,unsigned,unsigned,unsigned>> kmer_list;   // (i,j,a,m_a)   i=reference (Chromosome), j=position of matching k-mer in reference, a=abundance of k-mer in reference, m_a=minimizer_active_bases
std::vector<std::tuple<char,unsigned,unsigned,unsigned>>::const_iterator itrk;
// auto tbegin = std::chrono::high_resolution_clock::now();

std::string barcode;
std::string new_barcode;
std::string meta;
StringSet<Dna5String> reads;
typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
resize(reads, 2, Exact());
Dna5String read1;
Dna5String read2;
CharString id1;
CharString id2;

// opening read files
SeqFileIn file1(toCString(options.readfile1));
SeqFileIn file2(toCString(options.readfile2));

// preparing barcode Index
std::vector<std::string> BCI_barcodes;
std::vector<std::pair<std::streampos,std::streampos>> BCI_positions;
std::streampos BCI_pos1;
std::streampos BCI_pos2;

std::cerr << "Processing read file:\n";
tbegin = std::chrono::high_resolution_clock::now();

while (atEnd(file1)!=1) { // proceeding through files
  BCI_pos1=file1.stream.file.tellg();
  readRecord(id1, read1, file1);
  assignValue(reads,0,read1);
  meta=toCString(id1);
  new_barcode=meta.substr(meta.find("RX:Z:")+5,16);
  // std::cerr << "\n" << new_barcode << "\n";
  if (barcode!=new_barcode && !kmer_list.empty()) { //If Barcode changes: map kmer_list and reinitialize kmer_list
    //append Barcode Index
    BCI_pos2=file2.stream.file.tellg();
    BCI_barcodes.push_back(new_barcode);
    BCI_positions.push_back(std::make_pair(BCI_pos1,BCI_pos2));
    // map barcode and clear k_mer list
    sort(kmer_list.begin(),kmer_list.end());
    MapKmerList(kmer_list,max_window_size,max_gap_size,window_count,resultfile);
    kmer_list.clear();
    std::cerr << "barcode processed in: " << (float)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-tbegin).count()/1000 << "s\n";
    tbegin = std::chrono::high_resolution_clock::now();
  }

  readRecord(id2, read2, file2);
  assignValue(reads,1,read2);
  barcode=new_barcode;

  for (TStringSetIterator it = begin(reads); it!=end(reads); ++it){                                            // Iterating over the reads
    std::pair <long long int, long long int> hash = hashkMer(infix(*it,0,k),k);                                // calculation of the hash value for the first k-mer
    long long int minimizer_position=0;
    long long int minimizer = InitMini(infix(*it,0,mini_window_size), k, hash, maxhash, random_seed, minimizer_position);          // calculating the minimizer of the first window
    unsigned minimizer_active_bases=1;
    if (length(*it)>mini_window_size){
      for (unsigned t=0;t<(length(*it)-1-mini_window_size);t++){                                                   // iterating over all kmers
        if (t!=minimizer_position){                 // if old minimizer in current window
          rollinghashkMer(hash.first,hash.second,(*it)[t+mini_window_size],k,maxhash); // inline?!
          if (minimizer > ReturnSmaller(hash.first,hash.second,random_seed)){ // if new value replaces current minimizer
            AppendPos(kmer_list, minimizer, C, dir, pos, bucket_number,minimizer_active_bases);
            minimizer=ReturnSmaller(hash.first,hash.second,random_seed);
            minimizer_position=t+1+mini_window_size-k;
            minimizer_active_bases=0;
          }
          minimizer_active_bases++;
        }else{
          AppendPos(kmer_list, minimizer, C, dir, pos, bucket_number, minimizer_active_bases);                                                                                                  // if old minimizer no longer in window
          minimizer_position=t+1;
          hash=hashkMer(infix(*it,t+1,t+1+k),k);
          minimizer=InitMini(infix(*it,t+1,t+1+mini_window_size), k, hash, maxhash, random_seed, minimizer_position); // find minimizer in current window by reinitialization
          unsigned minimizer_active_bases=1;
        }
      }
      AppendPos(kmer_list, minimizer, C, dir, pos, bucket_number, minimizer_active_bases);   // append last minimizer                                                                                               // if old minimizer no longer in window
    }
  }
}
if (!kmer_list.empty()) {
  sort(kmer_list.begin(),kmer_list.end());
  MapKmerList(kmer_list,max_window_size,max_gap_size,window_count,resultfile);
  std::cerr << "barcode processed in: " << (float)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-tbegin).count()/1000 << "s\n";
  tbegin = std::chrono::high_resolution_clock::now();
}

close(file1);
close(file2);

std::cerr << "Writing BarcodeIndex to file ...";

// write Barcode Index to file
std::string IndBC=options.bci_name;
IndBC.append("_bc.txt");
IndPos=options.bci_name;
IndPos.append("_pos.txt");

std::ofstream file_bc;
file_bc.open(IndBC, std::ios::binary);
for (std::vector<std::string>::const_iterator it=BCI_barcodes.begin(); it!=BCI_barcodes.end(); it++){
  file_bc << *it << "\n";
}
file_bc.close();

std::ofstream file_pos;
file_pos.open(IndPos, std::ios::binary);
for (std::vector<std::pair<std::streampos,std::streampos>>::const_iterator it=BCI_positions.begin(); it!=BCI_positions.end(); it++){
  file_bc << std::get<0>(*it) << "\n" << std::get<1>(*it) << "\n";
}
file_pos.close();

std::cerr << "done!\n";

//Kontrollausgabe
// std::cerr << "\nKontrollausgabe:\n";
// std::cerr << "\n" << __LINE__;
//
// BCI_barcodes.clear();
// BCI_positions.clear();
// std::string testbarcode = "AAACACCGTAGATTAG";
//
// // open(file1,argv[1]);
// // open(file2,argv[2]);
// // std::string readfile1 = argv[1];
// // std::string readfile2 = argv[2];
// std::string Index_name = argv[4];
// std::cerr << "\n" << __LINE__;
//
// LoadBarcodeIndex(Index_name,BCI_barcodes,BCI_positions);
// std::cerr << "\n" << __LINE__;
//
// ReturnBarcodeReads(BCI_barcodes,BCI_positions,testbarcode,argv[1],argv[2]);
// std::cerr << "\n" << __LINE__;
//
// close(file1);
// close(file2);
return 0;
} //main


// maps k-mer list to reference genome and returns best fitting genomic windows
void MapKmerList(std::vector<std::tuple<char,unsigned,unsigned,unsigned>> & kmer_list, unsigned & max_window_size, unsigned & max_gap_size, unsigned & window_count, const char* file){

    std::vector<std::tuple<char,unsigned,unsigned,unsigned>>::const_iterator itrk;

    float lookQual[100]= {0,1024,6.24989, 0.624853, 0.195309, 0.0926038, 0.0541504, 0.0358415, 0.0257197, 0.0195267, 0.0154498, 0.0126139, 0.0105548, 0.00900754, 0.00781189, 0.0068662, 0.00610341, 0.00547777, 0.00495714, 0.00451843, 0.00414462, 0.003823, 0.00354385, 0.00329967, 0.00308456, 0.00289387, 0.00272383, 0.00257141, 0.00243412, 0.0023099, 0.00219705, 0.00209414, 0.00199997, 0.0019135, 0.00183386, 0.00176031, 0.0016922, 0.00162897, 0.00157012, 0.00151524, 0.00146395, 0.00141593, 0.00137087, 0.00132852, 0.00128865, 0.00125106, 0.00121556, 0.00118199, 0.00115019, 0.00112005, 0.00109142, 0.00106421, 0.00103832, 0.00101365, 0.000990122, 0.00096766, 0.000946195, 0.000925665, 0.00090601, 0.000887177, 0.000869117, 0.000851784, 0.000835136, 0.000819134, 0.000803742, 0.000788926, 0.000774656, 0.000760902, 0.000747638, 0.000734837, 0.000722477, 0.000710537, 0.000698994, 0.00068783, 0.000677027, 0.000666568, 0.000656437, 0.000646619, 0.0006371, 0.000627866, 0.000618906, 0.000610208, 0.00060176, 0.000593551, 0.000585573, 0.000577815, 0.000570269, 0.000562926, 0.000555778, 0.000548817, 0.000542037, 0.000535431, 0.000528992, 0.000522713, 0.000516589, 0.000510615, 0.000504785, 0.000499093, 0.000493536, 0.000488108};

    #define REF(X) std::get<0>(*(X))
    #define POS(X) std::get<1>(*(X))
    #define ABU(X) std::get<2>(*(X))
    #define ACT(X) std::get<3>(*(X))

    // std::cerr<<__LINE__<<"\n";
    std::vector<std::tuple<double,char,unsigned,unsigned>> best_windows(window_count,std::make_tuple(0,0,0,0)); //(maping_quality, reference, start position in referende, end position)
    std::vector<std::tuple<double,char,unsigned,unsigned>>::iterator itrbw;
    // std::cerr<<"iteration prepared. \n";

    char reference=REF(kmer_list.begin());
    std::vector<std::tuple<char,unsigned,unsigned,unsigned>>::const_iterator itrstart=kmer_list.begin();
    unsigned start_position=POS(kmer_list.begin());
    unsigned end_position=POS(kmer_list.begin());
    double window_quality=0;
    std::tuple<double,char,unsigned,unsigned> candidate=std::make_tuple(0,0,0,4294967295); //(maping_quality, reference, start position in referende, end position)

    if(ABU(kmer_list.begin())>99){        // calculating the quality of the first k-mer hit
      window_quality+=0.00032*ACT(kmer_list.begin());
    }else{
      window_quality+=lookQual[ABU(kmer_list.begin())]*ACT(kmer_list.begin()); // lookQual = (1/(log(abund)^5))*minimizer_active_bases
    }

    for(itrk=kmer_list.begin()+1;itrk!=kmer_list.end();itrk++){ //iterating over kmer listed

      if (/*end position*/std::get<3>(candidate) < start_position) { // if current window no longer overlaps the qualifiing window
        ReportWindow(best_windows,candidate);
      }

      if (reference==REF(itrk) && (POS(itrk)-start_position) < max_window_size && (POS(itrk)-end_position) < max_gap_size) { //checking break criteria
        //append window by kmer_hit
        if(ABU(itrk)>99){
          window_quality+=0.00032*ACT(itrk);
        }else{
          window_quality+=lookQual[ABU(itrk)]*ACT(itrk);
        }
        end_position=POS(itrk);

      }else if (REF(itrk)!=reference || (POS(itrk)-end_position) > max_gap_size){  // if k_mer hit from next reference or gapsize to large: report current window or candiadate window and initialize new window

        if(window_quality > std::get<0>(candidate)){ // report current window or candidate
          candidate=std::make_tuple(window_quality,reference,start_position,end_position);
          ReportWindow(best_windows,candidate);
        }else{
          ReportWindow(best_windows,candidate);
        }

        if(ABU(itrk)>99){ // initialize new window
          window_quality=0.00032*ACT(itrk);
        }else{
          window_quality=lookQual[ABU(itrk)]*ACT(itrk);
        }
        itrstart=itrk;
        reference=REF(itrk);
        start_position=POS(itrk);
        end_position=POS(itrk);
      }else{ // maximum window size criterion hurt: shrink window from start and save better one as candidate
        if (window_quality > std::get<0>(candidate)) { // check if current window better than candidate: if yes: replace candidate
          candidate=std::make_tuple(window_quality,reference,start_position,end_position);
        }
        if(ABU(itrk)>99){ //Append window by new kmer_hit
          window_quality+=0.00032*ACT(itrk);
        }else{
          window_quality+=lookQual[ABU(itrk)]*ACT(itrk);
        }
        end_position=POS(itrk);

        while (POS(itrk)-POS(itrstart)>max_window_size){ //shrinking window untill max_window_size criterion met
          if(ABU(itrstart)>99){
            window_quality-=0.00032*ACT(itrstart);
          }else{
            window_quality-=lookQual[ABU(itrstart)*ACT(itrstart)];
          }
          itrstart++;
        }
        start_position=POS(itrstart);
      }
    }

    candidate=std::make_tuple(window_quality,REF(itrk),POS(itrstart),POS(itrk));
    ReportWindow(best_windows,candidate); //reporting last window


    /*--------------------------------------------------------------------------------------------------*/


    //filter low quality windows
    double qualityThreshold=20000;
    if (std::get<0>(*(best_windows.end()-1))!=0) {
      while(std::get<0>(*best_windows.begin())<qualityThreshold && !best_windows.empty()){
        best_windows.erase(best_windows.begin());
      }
    }else{return;}

    // std::cerr << "len before: " << best_windows.size()<< "\t";
    // filter short windows
    // unsigned lengthThreshold=1000;
    // std::vector<int> toshort;
    // for (int i = best_windows.size()-1; i>=0; i--){
    //   if ((std::get<3>(best_windows[i])-std::get<2>(best_windows[i]))<lengthThreshold){
    //     toshort.push_back(i);
    //   }
    // }
    //
    // for (int i=(toshort.size()-1);i>=0;i--) {
    //   best_windows.erase(best_windows.begin()+toshort[i]);
    // }

    // std::cerr << "len during: " << best_windows.size()<< "\t";

    // std::cerr << "len after: " << best_windows.size()<< "\t";


    // std::cerr<<"best_windows found. ";
    // tend = std::chrono::high_resolution_clock::now();
    // std::cout <<"\ntime: "<< (float)std::chrono::duration_cast<std::chrono::milliseconds>(tend-tbegin).count()/1000 << " s\n";// << "ns" << std::endl;

    /*--------------------------------------------------------------------------------------------------*/
    // Output
    std::fstream results;
    results.open(file,std::ios::out::app);

    for(itrbw=best_windows.begin();itrbw!=best_windows.end(); itrbw++){

      std::string qual=std::to_string(std::get<0>(*itrbw));
      for (int i=qual.length();i<=18;i++) {qual+=" ";}
      std::string ref=std::to_string(std::get<1>(*itrbw));
      for (int i=ref.length();i<=8;i++) {ref+=" ";}
      std::string start=std::to_string(std::get<2>(*itrbw));
      for (int i=start.length();i<=13;i++) {start+=" ";}
      std::string end=std::to_string(std::get<3>(*itrbw));
      for (int i=end.length();i<=13;i++) {end+=" ";}
      std::string len=std::to_string(std::get<3>(*itrbw)-std::get<2>(*itrbw));
      for (int i=len.length();i<=13;i++) {len+=" ";}
      // results.write("\nquality: ",sizeof("\nquality: "));
      // results.write(&qual[0],sizeof(&qual[0]));
      // results.write("\tref: ",sizeof("\tref: "));
      // results.write(&ref[0],sizeof(&ref[0]));
      // results.write("\tstart: ",sizeof("\tstart: "));
      // results.write(&start[0],sizeof(&start[0]));
      // results.write("\tend: ",sizeof("\tend: "));
      // results.write(&end[0],sizeof(&end[0]));
      // results.write("\tlength: ",sizeof("\tlength: "));
      // results.write(&len[0],sizeof(&len[0]));
      results<<"\nquality: " << qual << "\tref: " << ref << "\tstart: "<< start << "\tend: " << end << "\tlength: " << len;
    }
    results.close();
    results<<"\n";
    return;
  } //MapKmerList
