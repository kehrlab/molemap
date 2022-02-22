# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <iostream>
# include <fstream>
# include "functions.h"
# include "map.h"
# include <time.h>
using namespace seqan;

/*
g++ BarcodeMapper.cpp -o bcmap
*/
void MapKmerList(std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, uint_fast32_t & max_window_size, uint_fast32_t & max_gap_size, uint_fast8_t & window_count, const char* file, std::string barcode, unsigned scoreThreshold, unsigned lengthThreshold, std::string & results, std::vector<std::string> & lookChrom, std::vector<uint32_t> & histogram);
uint_fast8_t getBarcodeLength(std::string & readfile1, std::streampos & readfile1_size);
bool SearchID(SeqFileIn & file, CharString id, std::streampos startpos, std::streampos endpos);
std::string skipToNextBarcode(SeqFileIn & file, CharString & id1, uint_fast8_t barcode_length);
void skipToNextBarcode2(SeqFileIn & file1, SeqFileIn & file2, std::string & barcode, uint_fast8_t barcode_length);
void trimmWindow(std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrstart, std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrk, std::tuple<double,uint_fast8_t,uint32_t,uint32_t> & candidate, std::vector<float> & lookQual);

struct bcmapOptions{
  std::string readfile1;
  std::string readfile2;
  std::string kmer_index_name;
  std::string read_index_name;
  unsigned k;
  unsigned mini_window_size;
  unsigned max_window_size;
  unsigned max_gap_size;
  std::string output_file;
  unsigned s;
  unsigned l;
  unsigned threads;
  bcmapOptions() :
  kmer_index_name("Index"), read_index_name("ReadIndex"), k(31), mini_window_size(61), max_window_size(300000), max_gap_size(20000),output_file("barcode_windows.bed"),l(10000) , s(0), threads(16)
  {}
  };

seqan::ArgumentParser::ParseResult parseCommandLine(bcmapOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bcmap map");

    // Define arguments.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "readfile1.fastq"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "readfile2.fastq"));
    // addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Index_name[IN]"));
    // addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Barcode_index_name[OUT]"));

    // Define Options
    addOption(parser, seqan::ArgParseOption(
        "i", "kmer_index_name", "Name of the folder in which the kmer index is stored.",
        seqan::ArgParseArgument::STRING, "kmer_index_name[IN]"));
    setDefaultValue(parser, "i", "Index");
    addOption(parser, seqan::ArgParseOption(
        "r", "Read_index_name", "Name of the ReadIndex.",
        seqan::ArgParseArgument::STRING, "Index_name[IN]"));
    setDefaultValue(parser, "r", "ReadIndex");
    addOption(parser, seqan::ArgParseOption(
        "k", "kmer_length", "Length of kmers in index.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "k", "31");
    setMinValue(parser, "k", "8");
    setMaxValue(parser, "k", "31");
    addOption(parser, seqan::ArgParseOption(
        "m", "mini_window_size", "Length of minimizing window.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "m", "61");
    addOption(parser, seqan::ArgParseOption(
        "w", "max_window_size", "Maximum length of genomic windows.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "w", "300000");
    addOption(parser, seqan::ArgParseOption(
        "g", "max_gap_size", "Maximum gap between minimizer hits of same genomic window.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "g", "20000");
    addOption(parser, seqan::ArgParseOption(
        "o", "output", "Path to the output file.",
        seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "o", "barcode_windows.bed");
    addOption(parser, seqan::ArgParseOption(
        "s", "score_threshold", "Minimum score threshold for genomic windows.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "s", "0");
    addOption(parser, seqan::ArgParseOption(
        "l", "length", "Length threshold for genomic windows.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "l", "10000");
    setMinValue(parser, "l", "5000");
    addOption(parser, seqan::ArgParseOption(
        "t", "threads", "Number of threads available.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "t", "16");

    seqan::addUsageLine(parser,"readfile.1.fq readfile.2.fq [OPTIONS]");
    setShortDescription(parser, "Map barcodes to reference.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addDescription(parser,
               "Barcodes will be mapped to reference genome. "
               "Returns genomic windows from which barcoded reads most likely originate. "
               "Each window is rated by a quality score. "
               "Requires readfiles to be sorted by barcode (use bcctools). "
               "Requires reference to be indexed using the 'map' command. ");
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;}

    // Extract argument and option values.
    getArgumentValue(options.readfile1, parser, 0);
    getArgumentValue(options.readfile2, parser, 1);
    // getArgumentValue(options.index_name, parser, 2);
    // getArgumentValue(options.bci_name, parser, 3);

    getOptionValue(options.kmer_index_name, parser, "i");
    getOptionValue(options.read_index_name, parser, "r");
    getOptionValue(options.k, parser, "k");
    getOptionValue(options.mini_window_size, parser, "m");
    getOptionValue(options.max_window_size, parser, "w");
    getOptionValue(options.max_gap_size, parser, "g");
    getOptionValue(options.output_file, parser, "o");
    getOptionValue(options.s, parser, "s");
    getOptionValue(options.l, parser, "l");
    getOptionValue(options.threads, parser, "t");

    return seqan::ArgumentParser::PARSE_OK;
}

int map(int argc, char const ** argv){

  // parsing command line arguments
  bcmapOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;
  std::cout <<'\n'
            << "readfile1        \t" << options.readfile1 << '\n'
            << "readfile2        \t" << options.readfile2 << '\n'
            << "kmer_index_name  \t" << options.kmer_index_name << '\n'
            << "output file      \t" << options.output_file << '\n'
            << "read_index_name  \t" << options.read_index_name << '\n'
            << "threads          \t" << options.threads << '\n'
            << "score threshold  \t" << options.s << '\n'
            << "k                \t" << options.k << '\n'
            << "minimizer window \t" << options.mini_window_size << '\n'
            << "max window size  \t" << options.max_window_size << '\n'
            << "max gap size     \t" << options.max_gap_size << '\n'
            << "length threshold \t" << options.l << '\n';

  uint_fast8_t k = options.k;

  int k_2;
  if (k>16){
    k_2=(k-15)*2;
  }else{
    k_2=0;
  }
  // int k_2 = k+1;

  uint_fast8_t mini_window_size = options.mini_window_size;

  // defining Parameters

  uint_fast32_t max_window_size=options.max_window_size; //200000;  //5000;   // maximum size of the genomic windows to wich the reads are matched
  uint_fast32_t max_gap_size=options.max_gap_size; //20000;     // maximum gap size between two adjacent k_mer hits
  uint_fast8_t window_count=50;   // amount of saved candidate windows

  //checking file size and barcode size
  SeqFileIn file1(toCString(options.readfile1));
  SeqFileIn file2(toCString(options.readfile2));
  file1.stream.file.seekg(0, std::ios::end);
  file2.stream.file.seekg(0, std::ios::end);
  std::streampos readfile1_size=file1.stream.file.tellg();
  std::streampos readfile2_size=file2.stream.file.tellg();
  close(file1);
  close(file2);
  uint_fast8_t barcode_length=getBarcodeLength(options.readfile1, readfile1_size);
  if(barcode_length == 0){
    return 0;
  }


  // reading the Index
  std::cerr << "Reading in the k-mer index";

  String<uint32_t> dir;
  String<uint32_t> pos;
  String<uint_fast8_t> ref;
  String<int32_t> C;

  std::string IndPos=options.kmer_index_name;
  IndPos.append("/pos.txt");
  std::string IndRef=options.kmer_index_name;
  IndRef.append("/ref.txt");
  std::string IndDir=options.kmer_index_name;
  IndDir.append("/dir.txt");
  std::string IndC=options.kmer_index_name;
  IndC.append("/C.txt");
  std::string IndFai=options.kmer_index_name;
  IndFai.append("/fai.txt");
  std::vector<std::string> lookChrom;


  //reading Index files in parallel
  auto tbegin = std::chrono::high_resolution_clock::now();

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

      std::ifstream input;
      input.open(toCString(IndFai), std::ios::in);
      std::string line;
      while(getline(input,line)){
        lookChrom.push_back(line);
      }
      input.close();

      std::cerr <<".";
    }
    if (i==2){
      String<uint32_t, External<> > extdir;
      if (!open(extdir, IndDir.c_str(), OPEN_RDONLY)){
        throw std::runtime_error("Could not open index directory file." );
      }
      assign(dir, extdir, Exact());
      close(extdir);
      std::cerr << ".";
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

  std::cerr <<"..done.\n";

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

  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  omp_lock_t lock;
  omp_init_lock(&lock);

  std::cerr << "Processing read file...";

  std::vector<uint32_t> histogram(200,0);

  #pragma omp parallel for ordered
  for (int t=0; t<options.threads; t++){
    //declare variables
    std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> kmer_list;   // (i,j,a,m_a)   i=reference (Chromosome), j=position of matching k-mer in reference, a=abundance of k-mer in reference, m_a=minimizer_active_bases
    StringSet<Dna5String> reads;
    resize(reads, 2, Exact());
    Dna5String read1;
    Dna5String read2;
    CharString id1;
    CharString id2;
    std::string barcode;
    std::string new_barcode;
    std::string results;
    std::streampos pos_temp;
    std::streampos BCI_1s;
    std::streampos BCI_2s;
    std::vector<std::tuple<std::string,std::streampos,std::streampos>> BCI_local; // (barcode, BCI_1s, BCI_2s, count)     #Deprecated:(barcode, BCI_1s, BCI_1e, BCI_2s, BCI_2e)
    std::vector<uint32_t> histogram_local(200,0);
    //open readfiles
    SeqFileIn file1(toCString(options.readfile1));
    SeqFileIn file2(toCString(options.readfile2));

    // calculate end and starposition in file for this thread
    std::streampos startpos=readfile1_size/options.threads*t;
    std::streampos endpos=readfile1_size/options.threads*(t+1);

    //move file 1 to start position
    if (t!=0){
      std::string line;
      file1.stream.file.seekg(startpos);
      while(line!="+"){
        getline(file1.stream.file,line);
      }
      file1.stream.file.ignore(10000,'\n');
      barcode=skipToNextBarcode(file1, id1, barcode_length);
      } else {
      file1.stream.file.seekg(0);
      try{
        readRecord(id1, read1, file1);
      }
      catch (Exception const & e){
        std::cerr << "ERROR: " << e.what() << std::endl;
      }
      // readRecord(id1, read1, file1);
      barcode=getBarcode(toCString(id1),barcode_length);
      file1.stream.file.seekg(0);
    }

    //align file2 with file1
    if (t!=0){
      std::streampos startpos2=readfile2_size*t/options.threads;
      std::streampos endpos1=startpos2;
      std::streampos startpos1=startpos2;
      std::streampos endpos2=startpos2;
      if(startpos2>100000){
        startpos1-=100000;
      }else{
        startpos1=0;
      }
      if(readfile2_size-startpos2>100000){
        endpos2+=100000;
      }else{
        endpos2=readfile2_size;
      }
      while(true){
        if(SearchID(file2, getID(toCString(id1)),startpos1, endpos1)){
          break;
        }
        if(SearchID(file2, getID(toCString(id1)),startpos2, endpos2)){
          break;
        }
        endpos1=startpos1;
        if(startpos1>100000){startpos1-=100000;}else{startpos1=0;}
        startpos2=endpos2;
        if(readfile2_size-endpos2>100000){endpos2+=100000;}else{endpos2=readfile2_size;}
      }
    }

    // skip to first valid barcode
    while (barcode[0]=='*' && !atEnd(file1)) {
      skipToNextBarcode2(file1,file2,barcode,barcode_length);
    }

    BCI_1s=file1.stream.file.tellg();
    BCI_2s=file2.stream.file.tellg();

    // if skiped beyond the boundaries of threads scope: end thread
    if (file1.stream.file.tellg()>endpos){
      continue;
    }

    //proceed through readfile untill endpos
    while (!atEnd(file1)) { // proceeding through files
      pos_temp=file1.stream.file.tellg();
      readRecord(id1, read1, file1);
      new_barcode=getBarcode(toCString(id1),barcode_length);
      if (barcode!=new_barcode){ //If Barcode changes: map kmer_list and reinitialize kmer_list
        //append Barcode Index
        // BCI_1e=pos_temp;
        // BCI_2e=file2.stream.file.tellg();
        BCI_local.push_back(std::make_tuple(barcode, BCI_1s, BCI_2s));
        BCI_1s=pos_temp;

        // map barcode and clear k_mer list
        if (!kmer_list.empty()) {
          std::sort(kmer_list.begin(),kmer_list.end());
          MapKmerList(kmer_list,max_window_size,max_gap_size,window_count,toCString(options.output_file),barcode, options.s, options.l, results, lookChrom, histogram_local);
          // std::cerr << "thread: " << t << "line:  "<<  __LINE__ << "\n";

          kmer_list.clear();
          if (results.size()>100000) {
            if (omp_test_lock(&lock)){
              std::fstream output;
              output.open(options.output_file,std::ios::out | std::ios::app);
              output << results;
              results="";
              output.close();
              omp_unset_lock(&lock);
            }
          }
        }

        if (file1.stream.file.tellg()>endpos){ // break if range for this thread exceeded
          break;
        }
        while (new_barcode[0]=='*' && !atEnd(file1)) { // proceed until valid barcode
          readRecord(id2, read2, file2);
          skipToNextBarcode2(file1,file2,new_barcode,barcode_length);
          BCI_1s=file1.stream.file.tellg();
          readRecord(id1, read1, file1);
        }
        BCI_2s=file2.stream.file.tellg();
      }

      //process next read
      readRecord(id2, read2, file2);
      assignValue(reads,0,read1);
      assignValue(reads,1,read2);
      barcode=new_barcode;
      for (TStringSetIterator it = begin(reads); it!=end(reads); ++it){  // Iterating over the two reads of the read pair
        // initialize minimizer class
        minimizer mini;
        minimizedSequence miniSeq(*it,k,mini_window_size,random_seed,maxhash);
        // return and save all minimizers to kmer_list
        while(!miniSeq.at_end){
          mini=miniSeq.pop();
          AppendPos(kmer_list, mini.value^random_seed, C, dir, ref, pos, bucket_number, mini.active_bases, k_2);
        }
      }

    }
    if (!kmer_list.empty()) { // only ever happens for the last thread on the last barcode
      // BCI_local.push_back(std::make_tuple(barcode, BCI_1s, readfile1_size, BCI_2s, readfile2_size));
      std::sort(kmer_list.begin(),kmer_list.end());
      MapKmerList(kmer_list,max_window_size,max_gap_size,window_count,toCString(options.output_file),barcode, options.s, options.l, results, lookChrom, histogram_local);
    }

    close(file1);
    close(file2);

    // writing Barcode index to file and summing up local histograms
    #pragma omp ordered
    {
      std::ofstream file_bci;
      file_bci.open(options.read_index_name , std::ios::app/*, std::ios::binary*/);
      for (int i=0; i<BCI_local.size(); i++){
        file_bci  << std::get<0>(BCI_local[i]) << "\t"
                  << std::get<1>(BCI_local[i]) << "\t"
                  << std::get<2>(BCI_local[i]) << "\n";
      }
      file_bci.close();

      for (int i=0; i<histogram.size(); i++){
        histogram[i]+=histogram_local[i];
      }
    }



    omp_set_lock(&lock);
    std::fstream output;
    output.open(options.output_file,std::ios::out | std::ios::app);
    output << results;
    results="";
    output.close();
    omp_unset_lock(&lock);

  }
  // write last entry of Barcode index
  std::ofstream file_bci;
  file_bci.open(options.read_index_name , std::ios::app/*, std::ios::binary*/);
  file_bci << "ZZZZZZZZZZZZZZZZ\t" << readfile1_size << "\t" << readfile2_size;
  file_bci.close();

  // write histigram to file
  std::ofstream file_histogram;
  file_histogram.open(options.output_file+".hist");
  for (int i=0; i<histogram.size(); i++){
    file_histogram << std::to_string(histogram[i]) << "\n";
  }
  file_histogram.close();

  std::cerr << ".........done.\n";

  std::cerr << "Barcodes mapped sucessfully!\n\n";

  return 0;
} //map(argc,argv)

// maps k-mer list to reference genome and returns best fitting genomic windows
void MapKmerList(std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, uint_fast32_t & max_window_size, uint_fast32_t & max_gap_size, uint_fast8_t & window_count, const char* file, std::string barcode, unsigned scoreThreshold, unsigned lengthThreshold, std::string & results, std::vector<std::string> & lookChrom, std::vector<uint32_t> & histogram){
  unsigned qualityThreshold=50;
  std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrk;

  // std::vector<float> lookQual={0,20,3.0, 0.75, 0.38, 0.24, 0.17, 0.14, 0.11, 0.09, 0.08};
  std::vector<float> lookQual={0,5,1.44, 0.91, 0.72, 0.62, 0.56, 0.51, 0.48, 0.46, 0.43, 0.42, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.35, 0.34, 0.33}; // log(1/abundance)

  #define REF(X) std::get<0>(*(X))
  #define POS(X) std::get<1>(*(X))
  #define ABU(X) std::get<2>(*(X))
  #define ACT(X) std::get<3>(*(X))

  std::vector<std::tuple<double,uint_fast8_t,uint32_t,uint32_t>> best_windows(window_count,std::make_tuple(0,0,0,0)); //(maping_quality, reference, start position in referende, end position)
  std::vector<std::tuple<double,uint_fast8_t,uint32_t,uint32_t>>::iterator itrbw;
  uint_fast8_t reference=REF(kmer_list.begin());
  uint_fast32_t start_position=POS(kmer_list.begin());
  uint_fast32_t end_position=POS(kmer_list.begin());
  std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrstart=kmer_list.begin();
  double window_quality=0;
  std::tuple<double,uint_fast8_t,uint32_t,uint32_t> candidate=std::make_tuple(0,0,0,4294967295); //(maping_quality, reference, start position in referende, end position)

  window_quality+=lookQual[ABU(kmer_list.begin())]*ACT(kmer_list.begin()); // lookQual = (1/(log(abund)^3))*minimizer_active_bases
  for(itrk=kmer_list.begin()+1;itrk!=kmer_list.end();itrk++){ //iterating over kmer listed

    if (/*end position*/std::get<3>(candidate) < start_position) { // if current window no longer overlaps the qualifiing window
      if (std::get<0>(candidate)>qualityThreshold){
        trimmWindow(kmer_list, itrstart, itrk, candidate, lookQual);
        ReportWindow(best_windows,candidate);
      }
    }

    if (reference==REF(itrk) && (POS(itrk)-start_position) < max_window_size && (POS(itrk)-end_position) < max_gap_size) { //checking break criteria
      //append window by kmer_hit
      window_quality+=lookQual[ABU(itrk)]*ACT(itrk);
      end_position=POS(itrk);

    }else if (REF(itrk)!=reference || (POS(itrk)-end_position) >= max_gap_size){  // if k_mer hit from next reference or gapsize to large: report current window or candiadate window and initialize new window

      if(window_quality > std::get<0>(candidate)){ // report current window or candidate
        candidate=std::make_tuple(window_quality,reference,start_position,end_position);
        if (std::get<0>(candidate)>qualityThreshold){
          trimmWindow(kmer_list, itrstart, itrk, candidate, lookQual);
          ReportWindow(best_windows,candidate);
        }
      }else{
        if (std::get<0>(candidate)>qualityThreshold){
          trimmWindow(kmer_list, itrstart, itrk, candidate, lookQual);
          ReportWindow(best_windows,candidate);
        }
      }
      window_quality=lookQual[ABU(itrk)]*ACT(itrk);
      itrstart=itrk;
      reference=REF(itrk);
      start_position=POS(itrk);
      end_position=POS(itrk);

    }else if ((POS(itrk)-start_position) >= max_window_size){ // maximum window size criterion hurt: shrink window from start and save better one as candidate
      if (window_quality > std::get<0>(candidate)) { // check if current window better than candidate: if yes: replace candidate
        candidate=std::make_tuple(window_quality,reference,start_position,end_position);
      }
      window_quality+=lookQual[ABU(itrk)]*ACT(itrk);
      end_position=POS(itrk);
      while (POS(itrk)-POS(itrstart)>max_window_size){ //shrinking window untill max_window_size criterion met
        window_quality-=lookQual[ABU(itrstart)]*ACT(itrstart);
        itrstart++;
      }
      start_position=POS(itrstart);

    }else{
      std::cerr << "\nREF(itrk): " << REF(itrk) << " reference: " << reference << "\nPOS(itrk): " << POS(itrk) << " start_position: " << start_position << " end_position: " << end_position << "\n";
      std::cerr << "POS(itrstart): " << POS(itrstart) << " REF(itrstart): " << REF(itrstart) << "\n";
    }
  }
  if(window_quality > std::get<0>(candidate)){ // report last window or last candidate
    candidate=std::make_tuple(window_quality,reference,start_position,end_position);
    if (std::get<0>(candidate)>qualityThreshold){
      trimmWindow(kmer_list, itrstart, itrk, candidate, lookQual);
      ReportWindow(best_windows,candidate);
    }
  }else{
    if (std::get<0>(candidate)>qualityThreshold){
      trimmWindow(kmer_list, itrstart, itrk, candidate, lookQual);
      ReportWindow(best_windows,candidate);
    }
  }
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

  for(itrbw=best_windows.begin();itrbw!=best_windows.end(); itrbw++){
    float qual=roundf((float)(std::get<0>(*itrbw)/(std::get<3>(*itrbw)-std::get<2>(*itrbw)))*100); // quality/length_of_mapping--> QualityPerBase
    if(qual>200){
      histogram[200]++;
    }else{
      histogram[qual]++;
    }
    if(qual>=scoreThreshold){
      std::string ref=lookChrom[std::get<1>(*itrbw)];
      std::string start=std::to_string(std::get<2>(*itrbw));
      std::string end=std::to_string(std::get<3>(*itrbw));
      results+=(ref + "\t"+ start + "\t" + end + "\t" + barcode + "\t" + std::to_string((int)qual) + "\n");
    }
  }
  return;
} //MapKmerList

uint_fast8_t getBarcodeLength(std::string & readfile1, std::streampos & readfile1_size){
  SeqFileIn file1(toCString(readfile1));
  Dna5String read1;
  CharString id;
  readRecord(id,read1,file1);
  std::string id1;
  std::string line="X";
  uint_fast8_t barcode_length;
  std::string barcode;
  for (uint64_t i=round((long long unsigned)readfile1_size/20); i<(uint64_t)readfile1_size;i+=round((uint64_t)readfile1_size/10)){ // check multiple positions throughout the readfile
    file1.stream.file.seekg(i);
    while(line!="+"){
      getline(file1.stream.file,line);
    }
    file1.stream.file.ignore(10000,'\n');
    readRecord(id, read1, file1);
    id1=toCString(id);
    barcode=id1.substr(id1.find("BX:Z:")+5,id1.find(' ',id1.find("BX:Z:")+5)-id1.find("BX:Z:")-5);              // determine BX:Z: entry
    if (barcode[0]=='A' || barcode[0]=='C' || barcode[0]=='T' || barcode[0]=='G'){
      close(file1);
      std::cerr << "barcode length   \t" << barcode.size() << "\n\n";
      return barcode.size();
    }
  }
  std::cerr << "\nERROR: Incorrect BX:Z: tag in readfile. Barcodes can not be determined. Please check your input file formating.\n";
  return 0;
}

//skips file1 to start of next barcode and returns the barcode
std::string skipToNextBarcode(SeqFileIn & file, CharString & id1, uint_fast8_t barcode_length){
  CharString id;
  Dna5String read;
  readRecord(id,read,file);
  std::string barcode=getBarcode(toCString(id),barcode_length);
  std::string new_barcode=barcode;
  std::streampos pos;
  while(barcode==new_barcode){
    pos=file.stream.file.tellg();
    readRecord(id,read,file);
    new_barcode=getBarcode(toCString(id),barcode_length);
  }
  file.stream.file.seekg(pos);
  id1=id;
  return new_barcode;
}

//skips both files till next barcode
void skipToNextBarcode2(SeqFileIn & file1, SeqFileIn & file2, std::string & barcode, uint_fast8_t barcode_length){
  CharString id;
  Dna5String read;
  std::streampos pos;
  pos=file1.stream.file.tellg();
  readRecord(id,read,file1);
  std::string new_barcode=getBarcode(toCString(id),barcode_length);
  while(new_barcode[0]=='*' && !atEnd(file1)){
    pos=file1.stream.file.tellg();
    readRecord(id,read,file1);
    new_barcode=getBarcode(toCString(id),barcode_length);
    readRecord(id,read,file2);
  }
  file1.stream.file.seekg(pos);
  barcode=new_barcode;
  return;
}

// searches for id in readfile and returns read and sets fileposition accordingly
bool SearchID(SeqFileIn & file, CharString id, std::streampos startpos, std::streampos endpos){
  CharString new_id;
  Dna5String read;
  std::string line;
  std::streampos pos;
  file.stream.file.seekg(startpos);
  while(line!="+"){
    getline(file.stream.file,line);
  }
  file.stream.file.ignore(10000,'\n');
  while(new_id!=id){
    if(pos>endpos){
      return false;
    }
    pos=file.stream.file.tellg();
    readRecord(new_id,read,file);
    new_id=getID(toCString(new_id));
  }
  file.stream.file.seekg(pos);
  return true;
}

// trimm edges of windows if low score per base | usage: trimmWindow(kmer_list, itrstart, itrk, candidate)
void trimmWindow(std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrstart, std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrk, std::tuple<double,uint_fast8_t,uint32_t,uint32_t> & candidate, std::vector<float> & lookQual){
  double edge_quality=0;
  double min_edge_qual=100; //200
  uint32_t edge_len=2000; //2000
  uint32_t counter=0;
  double window_quality=std::get<0>(candidate);
  itrk--;

  //trimm start
  while(edge_quality<min_edge_qual){
    if(itrstart+counter>=itrk){
      candidate=std::make_tuple(0,0,0,4294967295);
      return;
    }
    edge_quality+=lookQual[ABU(itrstart+counter)]*ACT(itrstart+counter);
    while(POS(itrstart+counter)-POS(itrstart)>edge_len){ //shrinking window until edge_len criterion met
      edge_quality-=lookQual[ABU(itrstart)]*ACT(itrstart);
      window_quality-=lookQual[ABU(itrstart)]*ACT(itrstart);
      itrstart++;
      counter--;
    }
    counter++;
  }
  edge_quality=0;
  counter=0;

  //trimm end
  while(edge_quality<min_edge_qual){
    if(itrk-counter<=itrstart){ //no area in candidate has the required minimum score/base
      candidate=std::make_tuple(0,0,0,4294967295);
      return;
    }
    edge_quality+=lookQual[ABU(itrk-counter)]*ACT(itrk-counter);
    while(POS(itrk)-POS(itrk-counter)>edge_len){ //shrinking window until edge_len criterion met
      edge_quality-=lookQual[ABU(itrk)]*ACT(itrk);
      window_quality-=lookQual[ABU(itrk)]*ACT(itrk);
      itrk--;
      counter--;
    }
    counter++;
  }
  candidate=std::make_tuple(window_quality,REF(itrstart),POS(itrstart),POS(itrk));
  return;
}
