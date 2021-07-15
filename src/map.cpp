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
void MapKmerList(std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, uint_fast32_t & max_window_size, uint_fast32_t & max_gap_size, uint_fast8_t & window_count, const char* file, std::string barcode, unsigned qualityThreshold, unsigned lengthThreshold, std::string & results, std::vector<std::string> & lookChrom);
std::string skipToNextBarcode(SeqFileIn & file, CharString & id1);
void skipToNextBarcode2(SeqFileIn & file1, SeqFileIn & file2, std::string & barcode);

void SearchID(SeqFileIn & file, CharString id, std::streampos startpos, std::streampos endpos);

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
  unsigned threads;
  bcmapOptions() :
  k(31), mini_window_size(61), output_file("barcode_windows.bed"),l(1000) , q(20000), threads(16)
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
        "i", "index_name", "Name of the folder in which the index is stored.",
        seqan::ArgParseArgument::STRING, "Index_name[IN]"));
    setDefaultValue(parser, "i", "Index");
    addOption(parser, seqan::ArgParseOption(
        "b", "Barcode_index_name", "Name of the BarcodeIndex.",
        seqan::ArgParseArgument::STRING, "Index_name[IN]"));
    setDefaultValue(parser, "b", "BarcodeIndex");
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

    getOptionValue(options.index_name, parser, "i");
    getOptionValue(options.bci_name, parser, "b");
    getOptionValue(options.k, parser, "k");
    getOptionValue(options.mini_window_size, parser, "m");
    getOptionValue(options.output_file, parser, "o");
    getOptionValue(options.q, parser, "q");
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
            << "index_name       \t" << options.index_name << '\n'
            << "barcodeindex_name\t" << options.bci_name << '\n'
            << "threads          \t" << options.threads << '\n'
            << "k                \t" << options.k << '\n'
            << "minimizer window \t" << options.mini_window_size << '\n'
            << "output file      \t" << options.output_file << '\n'
            << "quality threshold\t" << options.q << '\n'
            << "length threshold \t" << options.l << "\n\n";

  uint_fast8_t k = options.k;
  int k_2 = k+1;

  uint_fast8_t mini_window_size = options.mini_window_size;

  // defining Parameters

  uint_fast32_t max_window_size=300000;  //5000;   // maximum size of the genomic windows to wich the reads are matched
  uint_fast32_t max_gap_size=20000;     // maximum gap size between two adjacent k_mer hits
  uint_fast8_t window_count=100;   // amount of saved candidate windows

  // reading the Index
  std::cerr << "Reading in the k-mer index";

  String<uint32_t> dir;
  String<uint32_t> pos;
  String<uint_fast8_t> ref;
  String<int32_t> C;

  std::string IndPos=options.index_name;
  IndPos.append("/pos.txt");
  std::string IndRef=options.index_name;
  IndRef.append("/ref.txt");
  std::string IndDir=options.index_name;
  IndDir.append("/dir.txt");
  std::string IndC=options.index_name;
  IndC.append("/C.txt");
  std::string IndFai=options.index_name;
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

  SeqFileIn file1(toCString(options.readfile1));
  SeqFileIn file2(toCString(options.readfile2));
  file1.stream.file.seekg(0, std::ios::end);
  file2.stream.file.seekg(0, std::ios::end);
  std::streampos readfile1_size=file1.stream.file.tellg();
  std::streampos readfile2_size=file2.stream.file.tellg();
  close(file1);
  close(file2);

  /*
  Searching for all kmers of reads with the same Barcode
  */

  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  omp_lock_t lock;
  omp_init_lock(&lock);

  std::cerr << "Processing read file...";

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
    std::streampos BCI_1e;
    std::streampos BCI_2s;
    std::streampos BCI_2e;
    std::vector<std::tuple<std::string,std::streampos,std::streampos,std::streampos,std::streampos>> BCI_local; // (barcode, BCI_1s, BCI_1e, BCI_2s, BCI_2e)

    //open readfiles
    SeqFileIn file1(toCString(options.readfile1));
    SeqFileIn file2(toCString(options.readfile2));

    // calculate end and starposition in file for this thread
    std::streampos startpos=readfile1_size/options.threads*t;
    std::streampos endpos=readfile1_size/options.threads*(t+1);

    //move file 1 to start position
    if (t!=0){
      file1.stream.file.seekg(startpos);
      barcode=skipToNextBarcode(file1, id1);
    } else {
      file1.stream.file.seekg(0);
      try{
        readRecord(id1, read1, file1);
      }
      catch (Exception const & e){
        std::cerr << "ERROR: " << e.what() << std::endl;
      }
      // readRecord(id1, read1, file1);
      barcode=get10xBarcode(toCString(id1));
      file1.stream.file.seekg(0);
    }

    //align file2 with file1
    if(t!=0){
      startpos=(readfile2_size/options.threads*t)-(readfile2_size/options.threads/4);
    }
    SearchID(file2, getID(toCString(id1)), startpos, readfile2_size);

    // skip to first valid barcode
    while (barcode[0]=='*' && !atEnd(file1)) {
      skipToNextBarcode2(file1,file2,barcode);
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
      new_barcode=get10xBarcode(toCString(id1));
      if (barcode!=new_barcode){ //If Barcode changes: map kmer_list and reinitialize kmer_list
        //append Barcode Index
        BCI_1e=pos_temp;
        BCI_2e=file2.stream.file.tellg();
        BCI_local.push_back(std::make_tuple(barcode, BCI_1s, BCI_1e, BCI_2s, BCI_2e));
        BCI_1s=pos_temp;

        // map barcode and clear k_mer list
        if (!kmer_list.empty()) {
          sort(kmer_list.begin(),kmer_list.end());
          MapKmerList(kmer_list,max_window_size,max_gap_size,window_count,toCString(options.output_file),barcode, options.q, options.l, results, lookChrom);
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
        if (file1.stream.file.tellg()>endpos){
          break;
        }
        while (new_barcode[0]=='*' && !atEnd(file1)) {
          readRecord(id2, read2, file2);
          skipToNextBarcode2(file1,file2,new_barcode);
          BCI_1s=file1.stream.file.tellg();
          readRecord(id1, read1, file1);
        }
        BCI_2s=file2.stream.file.tellg();
      }

      readRecord(id2, read2, file2);
      assignValue(reads,0,read1);
      assignValue(reads,1,read2);
      barcode=new_barcode;

      for (TStringSetIterator it = begin(reads); it!=end(reads); ++it){                                            // Iterating over the reads
        std::pair <int64_t, int64_t> hash = hashkMer(infix(*it,0,k),k);                                // calculation of the hash value for the first k-mer
        int64_t minimizer_position=0;
        int64_t minimizer = InitMini(infix(*it,0,mini_window_size), k, hash, maxhash, random_seed, minimizer_position);          // calculating the minimizer of the first window
        uint_fast8_t minimizer_active_bases=1;
        if (length(*it)>mini_window_size){
          for (uint_fast32_t t=0;t<(length(*it)-1-mini_window_size);t++){               // iterating over all kmers
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
          }
          AppendPos(kmer_list, minimizer, C, dir, ref, pos, bucket_number, minimizer_active_bases,k_2);   // append last minimizer                                                                                               // if old minimizer no longer in window
        }
      }
    }
    if (!kmer_list.empty()) { // only ever happens for the last thread on the last barcode
      // BCI_local.push_back(std::make_tuple(barcode, BCI_1s, readfile1_size, BCI_2s, readfile2_size));
      sort(kmer_list.begin(),kmer_list.end());
      MapKmerList(kmer_list,max_window_size,max_gap_size,window_count,toCString(options.output_file),barcode, options.q, options.l, results, lookChrom);
    }

    close(file1);
    close(file2);

    // writing Barcode index to file
    #pragma omp ordered
    {
      std::ofstream file_bci;
      file_bci.open(options.bci_name , std::ios::app/*, std::ios::binary*/);
      for (int i=0; i<BCI_local.size(); i++){
        file_bci  << std::get<0>(BCI_local[i]) << "\t"
                  << std::get<1>(BCI_local[i]) << "\t"
                  << std::get<2>(BCI_local[i]) << "\t"
                  << std::get<3>(BCI_local[i]) << "\t"
                  << std::get<4>(BCI_local[i]) << "\n";
      }
      file_bci.close();
    }

    omp_set_lock(&lock);
    std::fstream output;
    output.open(options.output_file,std::ios::out | std::ios::app);
    output << results;
    results="";
    output.close();
    omp_unset_lock(&lock);

  }

  std::cerr << ".........done.\n";

  std::cerr << "Barcodes mapped sucessfully!\n";

  return 0;
} //map(argc,argv)

// maps k-mer list to reference genome and returns best fitting genomic windows
void MapKmerList(std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, uint_fast32_t & max_window_size, uint_fast32_t & max_gap_size, uint_fast8_t & window_count, const char* file, std::string barcode, unsigned qualityThreshold, unsigned lengthThreshold, std::string & results, std::vector<std::string> & lookChrom){

  std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrk;

  float lookQual[100]= {0,1024,6.24989, 0.624853, 0.195309, 0.0926038, 0.0541504, 0.0358415, 0.0257197, 0.0195267, 0.0154498, 0.0126139, 0.0105548, 0.00900754, 0.00781189, 0.0068662, 0.00610341, 0.00547777, 0.00495714, 0.00451843, 0.00414462, 0.003823, 0.00354385, 0.00329967, 0.00308456, 0.00289387, 0.00272383, 0.00257141, 0.00243412, 0.0023099, 0.00219705, 0.00209414, 0.00199997, 0.0019135, 0.00183386, 0.00176031, 0.0016922, 0.00162897, 0.00157012, 0.00151524, 0.00146395, 0.00141593, 0.00137087, 0.00132852, 0.00128865, 0.00125106, 0.00121556, 0.00118199, 0.00115019, 0.00112005, 0.00109142, 0.00106421, 0.00103832, 0.00101365, 0.000990122, 0.00096766, 0.000946195, 0.000925665, 0.00090601, 0.000887177, 0.000869117, 0.000851784, 0.000835136, 0.000819134, 0.000803742, 0.000788926, 0.000774656, 0.000760902, 0.000747638, 0.000734837, 0.000722477, 0.000710537, 0.000698994, 0.00068783, 0.000677027, 0.000666568, 0.000656437, 0.000646619, 0.0006371, 0.000627866, 0.000618906, 0.000610208, 0.00060176, 0.000593551, 0.000585573, 0.000577815, 0.000570269, 0.000562926, 0.000555778, 0.000548817, 0.000542037, 0.000535431, 0.000528992, 0.000522713, 0.000516589, 0.000510615, 0.000504785, 0.000499093, 0.000493536, 0.000488108};

  #define REF(X) std::get<0>(*(X))
  #define POS(X) std::get<1>(*(X))
  #define ABU(X) std::get<2>(*(X))
  #define ACT(X) std::get<3>(*(X))

  std::vector<std::tuple<double,uint_fast8_t,uint32_t,uint32_t>> best_windows(window_count,std::make_tuple(0,0,0,0)); //(maping_quality, reference, start position in referende, end position)
  std::vector<std::tuple<double,uint_fast8_t,uint32_t,uint32_t>>::iterator itrbw;

  uint_fast8_t reference=REF(kmer_list.begin());
  std::vector<std::tuple<uint_fast8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrstart=kmer_list.begin();
  uint_fast32_t start_position=POS(kmer_list.begin());
  uint_fast32_t end_position=POS(kmer_list.begin());
  double window_quality=0;
  std::tuple<double,uint_fast8_t,uint32_t,uint32_t> candidate=std::make_tuple(0,0,0,4294967295); //(maping_quality, reference, start position in referende, end position)


  window_quality+=lookQual[ABU(kmer_list.begin())]*ACT(kmer_list.begin()); // lookQual = (1/(log(abund)^5))*minimizer_active_bases

  for(itrk=kmer_list.begin()+1;itrk!=kmer_list.end();itrk++){ //iterating over kmer listed

    if (/*end position*/std::get<3>(candidate) < start_position) { // if current window no longer overlaps the qualifiing window
      ReportWindow(best_windows,candidate);
    }

    if (reference==REF(itrk) && (POS(itrk)-start_position) < max_window_size && (POS(itrk)-end_position) < max_gap_size) { //checking break criteria
      //append window by kmer_hit
      window_quality+=lookQual[ABU(itrk)]*ACT(itrk);
      end_position=POS(itrk);

    }else if (REF(itrk)!=reference || (POS(itrk)-end_position) >= max_gap_size){  // if k_mer hit from next reference or gapsize to large: report current window or candiadate window and initialize new window

      if(window_quality > std::get<0>(candidate)){ // report current window or candidate
        candidate=std::make_tuple(window_quality,reference,start_position,end_position);
        ReportWindow(best_windows,candidate);
      }else{
        ReportWindow(best_windows,candidate);
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
  candidate=std::make_tuple(window_quality,REF(itrk),POS(itrstart),POS(itrk));
  ReportWindow(best_windows,candidate); //reporting last window

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

    std::string qual=std::to_string((int)std::get<0>(*itrbw));
    std::string ref=lookChrom[std::get<1>(*itrbw)];
    std::string start=std::to_string(std::get<2>(*itrbw));
    std::string end=std::to_string(std::get<3>(*itrbw));
    results+=(ref + "\t"+ start + "\t" + end + "\t" + barcode + "\t" + qual + "\n");
  }

  return;
} //MapKmerList

//skips file1 to start of next barcode and returns the barcode
std::string skipToNextBarcode(SeqFileIn & file, CharString & id1){
  CharString id;
  Dna5String read;
  readRecord(id,read,file);
  std::string barcode=get10xBarcode(toCString(id));
  std::string new_barcode=barcode;
  std::streampos pos;
  while(barcode==new_barcode){
    pos=file.stream.file.tellg();
    readRecord(id,read,file);
    new_barcode=get10xBarcode(toCString(id));
  }
  file.stream.file.seekg(pos);
  id1=id;
  return new_barcode;
}

//skips both files till next barcode
void skipToNextBarcode2(SeqFileIn & file1, SeqFileIn & file2, std::string & barcode){
  CharString id;
  Dna5String read;
  std::streampos pos;
  pos=file1.stream.file.tellg();
  readRecord(id,read,file1);
  std::string new_barcode=get10xBarcode(toCString(id));
  while(new_barcode[0]=='*' && !atEnd(file1)){
    pos=file1.stream.file.tellg();
    readRecord(id,read,file1);
    new_barcode=get10xBarcode(toCString(id));
    readRecord(id,read,file2);
  }
  file1.stream.file.seekg(pos);
  barcode=new_barcode;
  return;
}

// searches for id in readfile and returns read and sets fileposition accordingly
void SearchID(SeqFileIn & file, CharString id, std::streampos startpos, std::streampos endpos){
  CharString new_id;
  Dna5String read;
  std::streampos pos;
  file.stream.file.seekg(startpos);
  while(new_id!=id){
    pos=file.stream.file.tellg();
    readRecord(new_id,read,file);
    new_id=getID(toCString(new_id));
  }
  file.stream.file.seekg(pos);
  return;
}
