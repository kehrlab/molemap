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
void MapKmerList(std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, uint_fast32_t & max_window_size, uint_fast32_t & max_gap_size, uint_fast8_t & window_count, const char* file, DnaString barcode, unsigned qualityThreshold, unsigned lengthThreshold);

struct bcmapOptions{
  std::string readfile1;
  std::string readfile2;
  std::string index_name;
  std::string bci_name;
  unsigned k;
  unsigned mini_window_size;
  std::string output_file;
  unsigned q;
  unsigned l;
  bcmapOptions() :
  k(31), mini_window_size(35), output_file("barcode_windows.bed"),l(1000) , q(20000)
  {}
  };

seqan::ArgumentParser::ParseResult parseCommandLine(bcmapOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bcmap");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "Path to readfile1.fastq"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "Path to readfile2.fastq"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Index_name[IN]"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Barcode_index_name[OUT]"));

    // Define Options
    addOption(parser, seqan::ArgParseOption(
        "k", "kmer_length", "Length of kmers in index.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "k", "31");
    setMinValue(parser, "k", "8");
    setMaxValue(parser, "k", "31");
    addOption(parser, seqan::ArgParseOption(
        "m", "mini_window_size", "Length of minimizing window.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "m", "35");
    addOption(parser, seqan::ArgParseOption(
        "o", "output", "Path to the output file.",
        seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "o", "barcode_windows.bed");
    addOption(parser, seqan::ArgParseOption(
        "q", "quality", "Quality threshold for genomic windows.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "q", "20000");
    addOption(parser, seqan::ArgParseOption(
        "l", "length", "Length threshold for genomic windows.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "l", "1000");

    setShortDescription(parser, "Map barcodes to reference.");
    setVersion(parser, "0.1");
    setDate(parser, "March 24 2021");
    addDescription(parser,
               "Barcodes will be mapped to reference genome."
               "Returns genomic windows from which barcoded reads most likely originate."
               "Each window is rated by a quality score."
               "Requires readfiles to be sorted by barcode (use bcctools)."
               "Requires reference to be indexed by 'countK'.");
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;}

    // Extract option values.
    getOptionValue(options.k, parser, "k");
    getOptionValue(options.mini_window_size, parser, "m");
    getOptionValue(options.output_file, parser, "o");

    getArgumentValue(options.readfile1, parser, 0);
    getArgumentValue(options.readfile2, parser, 1);
    getArgumentValue(options.index_name, parser, 2);
    getArgumentValue(options.bci_name, parser, 3);

    return seqan::ArgumentParser::PARSE_OK;
}

omp_lock_t out_lock;
// omp_init_lock(&out_lock);

int main(int argc, char const ** argv){

// parsing command line arguments
bcmapOptions options;
seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
if (res != seqan::ArgumentParser::PARSE_OK)
    return res == seqan::ArgumentParser::PARSE_ERROR;
std::cout <<'\n'
          << "readfile1        \t" << options.readfile1 << '\n'
          << "readfile2        \t" << options.readfile2 << '\n'
          << "index_name       \t" << options.index_name << '\n'
          << "barcodeindex_name\t" << options.bci_name << '\n'
          << "k                \t" << options.k << '\n'
          << "minimizer window \t" << options.mini_window_size << '\n'
          << "output file      \t" << options.output_file << '\n'
          << "quality threshold\t" << options.q << '\n'
          << "length threshold \t" << options.l << "\n\n";

uint_fast8_t k = options.k;
int k_2 = k+1;

uint_fast8_t mini_window_size = options.mini_window_size;

/*
defining Parameters
*/

uint_fast32_t max_window_size=200000;  //5000;   // maximum size of the genomic windows to wich the reads are matched
uint_fast32_t max_gap_size=20000;     // maximum gap size between two adjacent k_mer hits
uint_fast8_t window_count=100;   // amount of saved candidate windows


/*
reading the Index
*/

std::cerr << "Reading in the k-mer index";
auto tbegin = std::chrono::high_resolution_clock::now();

String<uint32_t> dir;
String<uint32_t> pos;
String<uint_fast8_t> ref;
String<int32_t> C;
//
std::string IndPos=options.index_name;
IndPos.append("_pos.txt");
std::string IndRef=options.index_name;
IndRef.append("_ref.txt");
std::string IndDir=options.index_name;
IndDir.append("_dir.txt");
std::string IndC=options.index_name;
IndC.append("_C.txt");

//reading Index files in parallel
#pragma omp parallel for
for(int i=0;i<4;i++){
  if (i==0){
    String<uint32_t, External<ExternalConfigLarge<>> > extpos;
    if (!open(extpos, IndPos.c_str(), OPEN_RDONLY)){
      throw std::runtime_error("Could not open index position file." );
    }
    assign(pos, extpos, Exact());
    close(extpos);
    std::cerr <<".";
  }
  if (i==1){
    String<uint_fast8_t, External<ExternalConfigLarge<>> > extref;
    if (!open(extref, IndRef.c_str(), OPEN_RDONLY)){
      throw std::runtime_error("Could not open index position file." );
    }
    assign(ref, extref, Exact());
    close(extref);
    std::cerr <<".";
  }
  if (i==2){
    String<uint32_t, External<> > extdir;
    if (!open(extdir, IndDir.c_str(), OPEN_RDONLY)){
      throw std::runtime_error("Could not open index directory file." );
    }
    assign(dir, extdir, Exact());
    close(extdir);
    std::cerr <<".";
  }
  if (i==3){
    String<int32_t, External<> > extC;
    if (!open(extC, IndC.c_str(), OPEN_RDONLY)){
      throw std::runtime_error("Could not open index counts file." );
    }
    assign(C, extC, Exact());
    close(extC);
    std::cerr << ".";
  }
} //for omp

// std::cerr << " in: " << (float)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-tbegin).count()/1000 << "s\n";

int64_t maxhash;
for (uint_fast8_t i=0;i<k;i++){
  maxhash= maxhash << 2 | 3;
}

std::srand(0);
int64_t random_seed=0;
for (uint_fast8_t i=0;i<k;++i){
  random_seed= random_seed << 2 | (int64_t)(std::rand()%3);
}

uint_fast32_t bucket_number=length(C);
// std::cerr << "bucket_number: " << bucket_number << "\n";


std::cerr <<"..done.\n";
// std::cerr << " in: " << (float)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-tbegin).count()/1000 << "s\n";

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

/*
Searching for all kmers of reads with the same Barcode
*/

// building the kmer_list for a specific Barcode
// std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> kmer_list;   // (i,j,a,m_a)   i=reference (Chromosome), j=position of matching k-mer in reference, a=abundance of k-mer in reference, m_a=minimizer_active_bases
std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrk;
// auto tbegin = std::chrono::high_resolution_clock::now();

std::string barcode;
std::string new_barcode;
std::string meta;
// StringSet<Dna5String> reads;
// resize(reads, 2, Exact());
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

// multithread handling
std::vector<std::vector<std::vector<Dna5String>>> readSet;
std::vector<std::vector<Dna5String>>::iterator itrreadSetG;
std::vector<Dna5String>::iterator it;
Dna5String read_overflow;
Dna5String barcode_overflow;
readSet.resize(3,{});
std::vector<std::vector<DnaString>> barcodeSet;
std::vector<DnaString>::iterator itrbarcG;
barcodeSet.resize(3,{});
std::vector<std::vector<std::streampos>> BCI_posSet;
BCI_posSet.resize(3,{});
uint32_t thread=0;    // "thread" for reading in reads from file1
uint32_t thread2=0;   // "thread" for reading in reads from file2
uint32_t thread3=0;   // "thread" for processing reads
uint32_t max_readCount=10000;
uint32_t readCount;
omp_lock_t lock;
omp_init_lock(&lock);
omp_init_lock(&out_lock);


std::cerr << "Processing read file...";

tbegin = std::chrono::high_resolution_clock::now();

// std::cerr << __LINE__ << "\n";

while (!atEnd(file1)){ //read first batch of reads from file1
  BCI_pos1=file1.stream.file.tellg();
  readRecord(id1, read1, file1);
  meta=toCString(id1);
  new_barcode=meta.substr(meta.find("RX:Z:")+5,16);
  if (barcode!=new_barcode){
    BCI_barcodes.push_back(new_barcode);
    BCI_posSet[thread].push_back(BCI_pos1);
    barcode=new_barcode;
    if (readCount>max_readCount){
      thread=(thread+1)%3; // iterate thread
      barcodeSet[thread].push_back(barcode);  //write barcode to Set of next batch
      readSet[thread].push_back({read1});
      readCount=0;
      break;
    }else{ //write read to readset of new barcode
      barcodeSet[thread].push_back(barcode);
      readSet[thread].push_back({read1});
    }
  }else{ //appned read to readset of current barcode
    readSet[thread].back().push_back(read1);
  }
  readCount++;
}

// std::cerr << __LINE__ << "\n";

#pragma omp parallel for  //read 2nd batch of reads from file1 and read first batch of reads from file2
for(int i=0;i<2;i++){

  if (i==0){  //read 2nd batch of reads from file1
    // std::cerr << __LINE__ << "\n";
    while (!atEnd(file1)){ //read first batch of reads from file1
      BCI_pos1=file1.stream.file.tellg();
      readRecord(id1, read1, file1);
      meta=toCString(id1);
      new_barcode=meta.substr(meta.find("RX:Z:")+5,16);
      if (barcode!=new_barcode){
        BCI_barcodes.push_back(new_barcode);
        BCI_posSet[thread].push_back(BCI_pos1);
        barcode=new_barcode;
        if (readCount>max_readCount){
          thread=(thread+1)%3; // iterate thread
          barcodeSet[thread].push_back(barcode);  //write barcode to Set of next batch
          readSet[thread].push_back({read1});
          readCount=0;
          break;
        }else{ //write read to readset of new barcode
          barcodeSet[thread].push_back(barcode);
          readSet[thread].push_back({read1});
        }
      }else{ //appned read to readset of current barcode
        readSet[thread].back().push_back(read1);
      }
      readCount++;
    }
    // std::cerr << __LINE__ << "\n";
  }

  if (i==1){  // read first batch of reads from file2
    // std::cerr << __LINE__ << "\n";
    for (uint32_t barc=0; barc<barcodeSet[thread2].size(); barc++){
      uint32_t r_count=readSet[thread][barc].size();
      BCI_pos2=file2.stream.file.tellg();
      BCI_positions.push_back(std::make_pair(BCI_posSet[thread2][barc],BCI_pos2));
      for (uint32_t read = 0; read < r_count; read++) {
        readRecord(id2, read2, file2);
        readSet[thread2][barc].push_back(read2);
      }
    }
    BCI_posSet[thread2].clear();
    thread2=(thread2+1)%3;
    // std::cerr << __LINE__ << "\n";

  }
}

// std::cerr << __LINE__ << "\n";

while (!atEnd(file1)){ // reading and processing next batch of reads until file endpos
  #pragma omp parallel sections
  {
  // for(int i=0;i<3;i++){

    // if (i==0){   // read next batch of reads from file1
    #pragma omp section
    {
      auto tbegin = std::chrono::high_resolution_clock::now();
      // std::cerr << __LINE__ << "\n";
      omp_set_lock(&lock);
      while (!atEnd(file1)){
        BCI_pos1=file1.stream.file.tellg();
        readRecord(id1, read1, file1);
        meta=toCString(id1);
        new_barcode=meta.substr(meta.find("RX:Z:")+5,16);
        if (barcode!=new_barcode){
          BCI_barcodes.push_back(new_barcode);
          BCI_posSet[thread].push_back(BCI_pos1);
          barcode=new_barcode;
          if (readCount>max_readCount){
            thread=(thread+1)%3; // iterate thread
            barcode_overflow=barcode;//write barcode to Set of next batch
            read_overflow=read1;
            readCount=0;
            break;
          }else{ //write read to readset of new barcode
            barcodeSet[thread].push_back(barcode);
            readSet[thread].push_back({read1});
          }
        }else{ //append read to readset of current barcode
          readSet[thread].back().push_back(read1);
        }
        readCount++;
      }
      omp_unset_lock(&lock);
      std::cerr << " reading file1 in: " << (float)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-tbegin).count()/1000 << "s\n";
      // std::cerr << __LINE__ << "\n";
    }

    // if (i==1){   // read next batch of reads from file2
    #pragma omp section
    {
      auto tbegin2 = std::chrono::high_resolution_clock::now();
      // std::cerr << __LINE__ << "\n";
      for (uint32_t barc=0; barc<barcodeSet[thread2].size(); barc++){
        uint32_t r_count=readSet[thread2][barc].size();
        BCI_pos2=file2.stream.file.tellg();
        BCI_positions.push_back(std::make_pair(BCI_posSet[thread2][barc],BCI_pos2));
        for (uint32_t read = 0; read < r_count; read++) {
          readRecord(id2, read2, file2);
          readSet[thread2][barc].push_back(read2);
        }
      }
      BCI_posSet[thread2].clear();
      thread2=(thread2+1)%3;
      std::cerr << " reading file2 in: " << (float)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-tbegin2).count()/1000 << "s\n";
      // std::cerr << __LINE__ << "\n";
    }
  }

    // if (i==2){   // process reads and write results to file
    //   #pragma omp section
    // {
      auto tbegin3 = std::chrono::high_resolution_clock::now();
      itrbarcG=barcodeSet[thread3].begin();
      itrreadSetG = readSet[thread3].begin();
      // #pragma omp parallel for
      #pragma omp parallel for
      for (int i=0; i!=barcodeSet[thread3].size(); i++) {// for all barcodes in set
        std::vector<std::vector<Dna5String>>::iterator itrreadSet = itrreadSetG+i;
        std::vector<DnaString>::iterator itrbarc = itrbarcG+i;
        // std::cerr << __LINE__ << "\n";
        std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> kmer_list;   // (i,j,a,m_a)   i=reference (Chromosome), j=position of matching k-mer in reference, a=abundance of k-mer in reference, m_a=minimizer_active_bases
        for (it = (*itrreadSet).begin(); it!=(*itrreadSet).end(); ++it){                                            // Iterating over the reads
          // std::cerr << __LINE__ << "\n";
          std::pair <int64_t, int64_t> hash = hashkMer(infix(*it,0,k),k);                                // calculation of the hash value for the first k-mer
          int64_t minimizer_position=0;
          int64_t minimizer = InitMini(infix(*it,0,mini_window_size), k, hash, maxhash, random_seed, minimizer_position);          // calculating the minimizer of the first window
          uint_fast8_t minimizer_active_bases=1;
          // std::cerr << __LINE__ << "\n";
          if (length(*it)>mini_window_size){
            // std::cerr << __LINE__ << "\n";
            for (uint_fast32_t t=0;t<(length(*it)-1-mini_window_size);t++){
              // std::cerr << __LINE__ << "\n";
              if (t!=minimizer_position){                 // if old minimizer in current window
                rollinghashkMer(hash.first,hash.second,(*it)[t+mini_window_size],k,maxhash); // inline?!
                if (minimizer > ReturnSmaller(hash.first,hash.second,random_seed)){ // if new value replaces current minimizer
                  AppendPos(kmer_list, minimizer, C, dir, ref, pos, bucket_number,minimizer_active_bases,k_2);
                  minimizer=ReturnSmaller(hash.first,hash.second,random_seed);
                  minimizer_position=t+1+mini_window_size-k;
                  minimizer_active_bases=0;
                }
                minimizer_active_bases++;
              }else{
                AppendPos(kmer_list, minimizer, C, dir, ref, pos, bucket_number, minimizer_active_bases,k_2);
                minimizer_position=t+1;
                hash=hashkMer(infix(*it,t+1,t+1+k),k);
                minimizer=InitMini(infix(*it,t+1,t+1+mini_window_size), k, hash, maxhash, random_seed, minimizer_position); // find minimizer in current window by reinitialization
                minimizer_active_bases=1;
              }
              // std::cerr << __LINE__ << "\n";
            }
            AppendPos(kmer_list, minimizer, C, dir, ref, pos, bucket_number, minimizer_active_bases,k_2);   // append last minimizer                                                                                               // if old minimizer no longer in window
          }
          // std::cerr << __LINE__ << "\n";
        } //for (itrreads = *(itrreadSetG).begin();
        // std::cerr << __LINE__ << "\n";
        if (!kmer_list.empty()) {
          // std::cerr << __LINE__ << "\n";
          sort(kmer_list.begin(),kmer_list.end());
          // std::cerr << __LINE__ << "\n";
          // std::cerr << "size readSet: " << readSet[thread3].size() << "\n";
          // std::cerr << "size barcodeSet: " << barcodeSet[thread3].size() << "\n";
          // std::cerr << "pos: " << (int)(itrreadSetG-readSet[thread3].begin()) << "\n";
          // std::cerr << "itrbarcGG: "<< *itrbarcG << "\n";
          MapKmerList(kmer_list,max_window_size,max_gap_size,window_count,toCString(options.output_file),*itrbarc, options.q, options.l);
          // std::cerr << __LINE__ << "\n";
        }
        // std::cerr << __LINE__ << "\n";
        // itrbarcG++;
      } //for (itrreadSetG = readSet[thread3].begin();
      // std::cerr << __LINE__ << "\n";
      readSet[thread3].clear();
      barcodeSet[thread3].clear();
      omp_set_lock(&lock);
      barcodeSet[thread3].push_back(barcode_overflow);
      readSet[thread3].push_back({read_overflow});
      omp_unset_lock(&lock);
      thread3=(thread3+1)%3;
      std::cerr << " processing reads in: " << (float)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-tbegin3).count()/1000 << "s\n";
      // std::cerr << __LINE__ << "\n";
    // } //if (i==2)
  // }
  // std::cerr << __LINE__ << "!!!!!!!!!!!!!\n";
// }   //pragma omp parallel

}   //while (!atEnd(file1))

// std::cerr << __LINE__ << "\n";

#pragma omp parallel for
for(int i=0;i<2;i++){

  if (i==0){  // read last batch of reads from file2
    for (uint32_t barc=0; barc<barcodeSet[thread2].size(); barc++){
      uint32_t r_count=readSet[thread][barc].size();
      BCI_pos2=file2.stream.file.tellg();
      BCI_positions.push_back(std::make_pair(BCI_posSet[thread2][barc],BCI_pos2));
      for (uint32_t read = 0; read < r_count; read++) {
        readRecord(id2, read2, file2);
        readSet[thread2][barc].push_back(read2);
      }
    }
    BCI_posSet[thread2].clear();
    thread2=(thread2+1)%3;
  }

  if (i==1){   // process reads and write results to file
    // itrbarcG=barcodeSet[thread3].begin();
    for (itrreadSetG = readSet[thread3].begin(), itrbarcG=barcodeSet[thread3].begin(); itrreadSetG != readSet[thread3].end(); itrreadSetG++ ,itrbarcG++) {// for all barcodes in set
      std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> kmer_list;   // (i,j,a,m_a)   i=reference (Chromosome), j=position of matching k-mer in reference, a=abundance of k-mer in reference, m_a=minimizer_active_bases
      for (it = (*itrreadSetG).begin(); it!=(*itrreadSetG).end(); ++it){                                            // Iterating over the reads
        std::pair <int64_t, int64_t> hash = hashkMer(infix(*it,0,k),k);                                // calculation of the hash value for the first k-mer
        int64_t minimizer_position=0;
        int64_t minimizer = InitMini(infix(*it,0,mini_window_size), k, hash, maxhash, random_seed, minimizer_position);          // calculating the minimizer of the first window
        uint_fast8_t minimizer_active_bases=1;
        // std::cerr << __LINE__ << "\n";
        if (length(*it)>mini_window_size){
          for (uint_fast32_t t=0;t<(length(*it)-1-mini_window_size);t++){
            // std::cerr << __LINE__ << "\n";                                                  // iterating over all kmers
            if (t!=minimizer_position){                 // if old minimizer in current window
              rollinghashkMer(hash.first,hash.second,(*it)[t+mini_window_size],k,maxhash); // inline?!
              // std::cerr << __LINE__ << "\n";
              if (minimizer > ReturnSmaller(hash.first,hash.second,random_seed)){ // if new value replaces current minimizer
                AppendPos(kmer_list, minimizer, C, dir, ref, pos, bucket_number,minimizer_active_bases,k_2);
                minimizer=ReturnSmaller(hash.first,hash.second,random_seed);
                minimizer_position=t+1+mini_window_size-k;
                minimizer_active_bases=0;
              }
              minimizer_active_bases++;
            }else{
              // std::cerr << __LINE__ << "\n";
              AppendPos(kmer_list, minimizer, C, dir, ref, pos, bucket_number, minimizer_active_bases,k_2);
              // std::cerr << __LINE__ << "\n";                                                                                                  // if old minimizer no longer in window
              minimizer_position=t+1;
              // std::cerr << __LINE__ << "\n";
              hash=hashkMer(infix(*it,t+1,t+1+k),k);
              // std::cerr << __LINE__ << "\n";
              minimizer=InitMini(infix(*it,t+1,t+1+mini_window_size), k, hash, maxhash, random_seed, minimizer_position); // find minimizer in current window by reinitialization
              // std::cerr << __LINE__ << "\n";
              minimizer_active_bases=1;
              // std::cerr << __LINE__ << "\n";
            }
          }
          AppendPos(kmer_list, minimizer, C, dir, ref, pos, bucket_number, minimizer_active_bases,k_2);   // append last minimizer                                                                                               // if old minimizer no longer in window
        }
      } //for (it = *(itrreadSetG).begin();
      if (!kmer_list.empty()) {
        sort(kmer_list.begin(),kmer_list.end());
        MapKmerList(kmer_list,max_window_size,max_gap_size,window_count,toCString(options.output_file),*itrbarcG, options.q, options.l);
      }
      // itrbarcG++;
    } //for (itet = readSet[thread3].begin();
    readSet[thread3].clear();
    barcodeSet[thread3].clear();
    thread3=(thread3+1)%3;

  } //if (i==1)
}

// std::cerr << __LINE__ << "\n";

// process last batch of reads
// itrbarcG=barcodeSet[thread3].begin();
for (itrreadSetG = readSet[thread3].begin(), itrbarcG=barcodeSet[thread3].begin(); itrreadSetG != readSet[thread3].end(); itrreadSetG++ ,itrbarcG++) {// for all barcodes in set
  // std::cerr << __LINE__ << "\n";
  std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> kmer_list;   // (i,j,a,m_a)   i=reference (Chromosome), j=position of matching k-mer in reference, a=abundance of k-mer in reference, m_a=minimizer_active_bases
  for (it = (*itrreadSetG).begin(); it!=(*itrreadSetG).end(); ++it){                                            // Iterating over the reads
    // std::cerr << __LINE__ << "\n";
    std::pair <int64_t, int64_t> hash = hashkMer(infix(*it,0,k),k);                                // calculation of the hash value for the first k-mer
    int64_t minimizer_position=0;
    int64_t minimizer = InitMini(infix(*it,0,mini_window_size), k, hash, maxhash, random_seed, minimizer_position);          // calculating the minimizer of the first window
    uint_fast8_t minimizer_active_bases=1;
    // std::cerr << __LINE__ << "\n";
    if (length(*it)>mini_window_size){
      for (uint_fast32_t t=0;t<(length(*it)-1-mini_window_size);t++){
        // std::cerr << __LINE__ << "\n";                                                  // iterating over all kmers
        if (t!=minimizer_position){                 // if old minimizer in current window
          rollinghashkMer(hash.first,hash.second,(*it)[t+mini_window_size],k,maxhash); // inline?!
          // std::cerr << __LINE__ << "\n";
          if (minimizer > ReturnSmaller(hash.first,hash.second,random_seed)){ // if new value replaces current minimizer
            AppendPos(kmer_list, minimizer, C, dir, ref, pos, bucket_number,minimizer_active_bases,k_2);
            minimizer=ReturnSmaller(hash.first,hash.second,random_seed);
            minimizer_position=t+1+mini_window_size-k;
            minimizer_active_bases=0;
          }
          minimizer_active_bases++;
        }else{
          // std::cerr << __LINE__ << "\n";
          AppendPos(kmer_list, minimizer, C, dir, ref, pos, bucket_number, minimizer_active_bases,k_2);
          // std::cerr << __LINE__ << "\n";                                                                                                  // if old minimizer no longer in window
          minimizer_position=t+1;
          // std::cerr << __LINE__ << "\n";
          hash=hashkMer(infix(*it,t+1,t+1+k),k);
          // std::cerr << __LINE__ << "\n";
          minimizer=InitMini(infix(*it,t+1,t+1+mini_window_size), k, hash, maxhash, random_seed, minimizer_position); // find minimizer in current window by reinitialization
          // std::cerr << __LINE__ << "\n";
          minimizer_active_bases=1;
          // std::cerr << __LINE__ << "\n";
        }
      }
      AppendPos(kmer_list, minimizer, C, dir, ref, pos, bucket_number, minimizer_active_bases,k_2);   // append last minimizer                                                                                               // if old minimizer no longer in window
    }
  } //for (itrreads = *(itrreadSetG).begin();
  // std::cerr << __LINE__ << "\n";
  if (!kmer_list.empty()) {
    // std::cerr << __LINE__ << "\n";
    sort(kmer_list.begin(),kmer_list.end());
    // std::cerr << "Barcode: " << *itrbarcG << "\n";
    MapKmerList(kmer_list,max_window_size,max_gap_size,window_count,toCString(options.output_file),*itrbarcG, options.q, options.l);
  }
  // std::cerr << __LINE__ << "\n";
  // itrbarcG++;
} //for (itrreadSetG = readSet[thread3].begin();

readSet[thread3].clear();
barcodeSet[thread3].clear();
thread3=(thread3+1)%3;

// std::cerr << __LINE__ << "\n";

close(file1);
close(file2);
// std::cerr << "\nbarcode processed in: " << (float)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-tbegin).count()/1000 << "s";
// tbegin = std::chrono::high_resolution_clock::now();
std::cerr << ".........done.\n";
std::cerr << "Writing BarcodeIndex to file...";

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

std::cerr << ".done.\n";
std::cerr << "Barcodes mapped sucessfully!\n";

// Kontrollausgabe
// std::cerr << "\nKontrollausgabe:\n";
//
// BCI_barcodes.clear();
// BCI_positions.clear();
// std::string testbarcode = "AAACACCGTAGATTAG";
//
// LoadBarcodeIndex(options.bci_name,BCI_barcodes,BCI_positions);
// std::vector<std::pair<Dna5String,Dna5String>> barcodedreads;
// barcodedreads=ReturnBarcodeReads(BCI_barcodes,BCI_positions,testbarcode,toCString(options.readfile1),toCString(options.readfile2));
// for (int i=0;i<barcodedreads.size();i++){
//   std::cerr << std::get<0>(barcodedreads[i]) << "\n" << std::get<1>(barcodedreads[i]) << "\n\n";
// }
// close(file1);
// close(file2);

return 0;
} //main


// maps k-mer list to reference genome and returns best fitting genomic windows
void MapKmerList(std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, uint_fast32_t & max_window_size, uint_fast32_t & max_gap_size, uint_fast8_t & window_count, const char* file, DnaString barcode, unsigned qualityThreshold, unsigned lengthThreshold){
    // std::cerr << __LINE__ << "\n";

    std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrk;

    float lookQual[100]= {0,1024,6.24989, 0.624853, 0.195309, 0.0926038, 0.0541504, 0.0358415, 0.0257197, 0.0195267, 0.0154498, 0.0126139, 0.0105548, 0.00900754, 0.00781189, 0.0068662, 0.00610341, 0.00547777, 0.00495714, 0.00451843, 0.00414462, 0.003823, 0.00354385, 0.00329967, 0.00308456, 0.00289387, 0.00272383, 0.00257141, 0.00243412, 0.0023099, 0.00219705, 0.00209414, 0.00199997, 0.0019135, 0.00183386, 0.00176031, 0.0016922, 0.00162897, 0.00157012, 0.00151524, 0.00146395, 0.00141593, 0.00137087, 0.00132852, 0.00128865, 0.00125106, 0.00121556, 0.00118199, 0.00115019, 0.00112005, 0.00109142, 0.00106421, 0.00103832, 0.00101365, 0.000990122, 0.00096766, 0.000946195, 0.000925665, 0.00090601, 0.000887177, 0.000869117, 0.000851784, 0.000835136, 0.000819134, 0.000803742, 0.000788926, 0.000774656, 0.000760902, 0.000747638, 0.000734837, 0.000722477, 0.000710537, 0.000698994, 0.00068783, 0.000677027, 0.000666568, 0.000656437, 0.000646619, 0.0006371, 0.000627866, 0.000618906, 0.000610208, 0.00060176, 0.000593551, 0.000585573, 0.000577815, 0.000570269, 0.000562926, 0.000555778, 0.000548817, 0.000542037, 0.000535431, 0.000528992, 0.000522713, 0.000516589, 0.000510615, 0.000504785, 0.000499093, 0.000493536, 0.000488108};

    std::vector<std::string> lookChrom={"chr1", "chr10", "chr11", "chr11_KI270721v1_random", "chr12", "chr13", "chr14", "chr14_GL000009v2_random", "chr14_GL000194v1_random", "chr14_GL000225v1_random", "chr14_KI270722v1_random", "chr14_KI270723v1_random", "chr14_KI270724v1_random", "chr14_KI270725v1_random", "chr14_KI270726v1_random", "chr15", "chr15_KI270727v1_random", "chr16", "chr16_KI270728v1_random", "chr17", "chr17_GL000205v2_random", "chr17_KI270729v1_random", "chr17_KI270730v1_random", "chr18", "chr19", "chr1_KI270706v1_random", "chr1_KI270707v1_random", "chr1_KI270708v1_random", "chr1_KI270709v1_random", "chr1_KI270710v1_random", "chr1_KI270711v1_random", "chr1_KI270712v1_random", "chr1_KI270713v1_random", "chr1_KI270714v1_random", "chr2", "chr20", "chr21", "chr22", "chr22_KI270731v1_random", "chr22_KI270732v1_random", "chr22_KI270733v1_random", "chr22_KI270734v1_random", "chr22_KI270735v1_random", "chr22_KI270736v1_random", "chr22_KI270737v1_random", "chr22_KI270738v1_random", "chr22_KI270739v1_random", "chr2_KI270715v1_random", "chr2_KI270716v1_random", "chr3", "chr3_GL000221v1_random", "chr4", "chr4_GL000008v2_random", "chr5", "chr5_GL000208v1_random", "chr6", "chr7", "chr8", "chr9", "chr9_KI270717v1_random", "chr9_KI270718v1_random", "chr9_KI270719v1_random", "chr9_KI270720v1_random", "chrM", "chrUn_GL000195v1", "chrUn_GL000213v1", "chrUn_GL000214v1", "chrUn_GL000216v2", "chrUn_GL000218v1", "chrUn_GL000219v1", "chrUn_GL000220v1", "chrUn_GL000224v1", "chrUn_GL000226v1", "chrUn_KI270302v1", "chrUn_KI270303v1", "chrUn_KI270304v1", "chrUn_KI270305v1", "chrUn_KI270310v1", "chrUn_KI270311v1", "chrUn_KI270312v1", "chrUn_KI270315v1", "chrUn_KI270316v1", "chrUn_KI270317v1", "chrUn_KI270320v1", "chrUn_KI270322v1", "chrUn_KI270329v1", "chrUn_KI270330v1", "chrUn_KI270333v1", "chrUn_KI270334v1", "chrUn_KI270335v1", "chrUn_KI270336v1", "chrUn_KI270337v1", "chrUn_KI270338v1", "chrUn_KI270340v1", "chrUn_KI270362v1", "chrUn_KI270363v1", "chrUn_KI270364v1", "chrUn_KI270366v1", "chrUn_KI270371v1", "chrUn_KI270372v1", "chrUn_KI270373v1", "chrUn_KI270374v1", "chrUn_KI270375v1", "chrUn_KI270376v1", "chrUn_KI270378v1", "chrUn_KI270379v1", "chrUn_KI270381v1", "chrUn_KI270382v1", "chrUn_KI270383v1", "chrUn_KI270384v1", "chrUn_KI270385v1", "chrUn_KI270386v1", "chrUn_KI270387v1", "chrUn_KI270388v1", "chrUn_KI270389v1", "chrUn_KI270390v1", "chrUn_KI270391v1", "chrUn_KI270392v1", "chrUn_KI270393v1", "chrUn_KI270394v1", "chrUn_KI270395v1", "chrUn_KI270396v1", "chrUn_KI270411v1", "chrUn_KI270412v1", "chrUn_KI270414v1", "chrUn_KI270417v1", "chrUn_KI270418v1", "chrUn_KI270419v1", "chrUn_KI270420v1", "chrUn_KI270422v1", "chrUn_KI270423v1", "chrUn_KI270424v1", "chrUn_KI270425v1", "chrUn_KI270429v1", "chrUn_KI270435v1", "chrUn_KI270438v1", "chrUn_KI270442v1", "chrUn_KI270448v1", "chrUn_KI270465v1", "chrUn_KI270466v1", "chrUn_KI270467v1", "chrUn_KI270468v1", "chrUn_KI270507v1", "chrUn_KI270508v1", "chrUn_KI270509v1", "chrUn_KI270510v1", "chrUn_KI270511v1", "chrUn_KI270512v1", "chrUn_KI270515v1", "chrUn_KI270516v1", "chrUn_KI270517v1", "chrUn_KI270518v1", "chrUn_KI270519v1", "chrUn_KI270521v1", "chrUn_KI270522v1", "chrUn_KI270528v1", "chrUn_KI270529v1", "chrUn_KI270530v1", "chrUn_KI270538v1", "chrUn_KI270539v1", "chrUn_KI270544v1", "chrUn_KI270548v1", "chrUn_KI270579v1", "chrUn_KI270580v1", "chrUn_KI270581v1", "chrUn_KI270582v1", "chrUn_KI270583v1", "chrUn_KI270584v1", "chrUn_KI270587v1", "chrUn_KI270588v1", "chrUn_KI270589v1", "chrUn_KI270590v1", "chrUn_KI270591v1", "chrUn_KI270593v1", "chrUn_KI270741v1", "chrUn_KI270742v1", "chrUn_KI270743v1", "chrUn_KI270744v1", "chrUn_KI270745v1", "chrUn_KI270746v1", "chrUn_KI270747v1", "chrUn_KI270748v1", "chrUn_KI270749v1", "chrUn_KI270750v1", "chrUn_KI270751v1", "chrUn_KI270752v1", "chrUn_KI270753v1", "chrUn_KI270754v1", "chrUn_KI270755v1", "chrUn_KI270756v1", "chrUn_KI270757v1", "chrX", "chrY", "chrY_KI270740v1_random"};

    #define REF(X) std::get<0>(*(X))
    #define POS(X) std::get<1>(*(X))
    #define ABU(X) std::get<2>(*(X))
    #define ACT(X) std::get<3>(*(X))

    std::vector<std::tuple<double,uint_fast8_t,uint32_t,uint32_t>> best_windows(window_count,std::make_tuple(0,0,0,0)); //(maping_quality, reference, start position in referende, end position)
    std::vector<std::tuple<double,uint_fast8_t,uint32_t,uint32_t>>::iterator itrbw;
    // std::cerr<<"iteration prepared. \n";

    uint_fast8_t reference=REF(kmer_list.begin());
    std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrstart=kmer_list.begin();
    uint_fast32_t start_position=POS(kmer_list.begin());
    uint_fast32_t end_position=POS(kmer_list.begin());
    double window_quality=0;
    std::tuple<double,uint_fast8_t,uint32_t,uint32_t> candidate=std::make_tuple(0,0,0,4294967295); //(maping_quality, reference, start position in referende, end position)
    // std::cerr << __LINE__ << "\n";

    if(ABU(kmer_list.begin())>99){        // calculating the quality of the first k-mer hit
      window_quality+=0.00032*ACT(kmer_list.begin());
    }else{
      window_quality+=lookQual[ABU(kmer_list.begin())]*ACT(kmer_list.begin()); // lookQual = (1/(log(abund)^5))*minimizer_active_bases
    }
    // std::cerr << __LINE__ << "\n";

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
    // std::cerr << __LINE__ << "\n";

    //filter low quality windows
    if (std::get<0>(*(best_windows.end()-1))!=0) {
      while(std::get<0>(*best_windows.begin())<qualityThreshold && !best_windows.empty()){
        best_windows.erase(best_windows.begin());
      }
    }else{return;}

    // filter short windows
    std::vector<int> toshort;
    for (int i = 0; i!=best_windows.size(); i++){
      if ((std::get<3>(best_windows[i])-std::get<2>(best_windows[i]))<lengthThreshold){
        toshort.push_back(i);
      }
    }
    for (int i=(toshort.size()-1);i>=0;--i) {
      best_windows.erase(best_windows.begin()+toshort[i]);
    }
    // std::cerr << __LINE__ << "\n";

    // Output

    omp_set_lock(&out_lock);

    std::fstream results;
    results.open(file,std::ios::out | std::ios::app);

    for(itrbw=best_windows.begin();itrbw!=best_windows.end(); itrbw++){

      std::string qual=std::to_string((int)std::get<0>(*itrbw));
      std::string ref=lookChrom[std::get<1>(*itrbw)];
      std::string start=std::to_string(std::get<2>(*itrbw));
      std::string end=std::to_string(std::get<3>(*itrbw));
      std::string len=std::to_string(std::get<3>(*itrbw)-std::get<2>(*itrbw));
      results<< ref << "\t"<< start << "\t" << end <<"\t" << barcode <<"\t" << qual <<"\t" << len << "\n";
      // results<< "ref: " << ref << "\tstart: "<< start << "\tend: " << end <<"\tbarcode: " << barcode <<"\tquality: " << qual <<"\tlength: " << len << "\n";
    }
    results.close();
    omp_unset_lock(&out_lock);
    // std::cerr << __LINE__ << "\n";
    return;
  } //MapKmerList
