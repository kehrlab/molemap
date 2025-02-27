# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/bam_io.h>
# include <seqan/arg_parse.h>
# include <iostream>
# include <fstream>
# include <time.h>
# include "map_linked.h"
# include "parser.h"
# include "functions.h"
# include "coverage_analysis.h"
# include "IO_functions.h"
# include "mapping_functions.h"
# include "index.h"
using namespace seqan;

std::string skipToNextBarcode(SeqFileIn & file, CharString & id1, uint8_t barcode_length);
void skipToNextBarcode2(SeqFileIn & file1, SeqFileIn & file2, std::string & barcode, uint8_t barcode_length);
void writeLastReadFileIndexEntrys(std::streampos & readfile1_size, std::streampos & readfile2_size, mapOptions & options);

int maplinked(int argc, char const ** argv){

  // parsing command line arguments
  mapOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine_map(options, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK)
      return res;
  printParseResults_map(options);

  // reading the Index
  std::cerr << "Reading in the k-mer index";

  openAddressingKmerHashtable Index;

  readKmerIndex(Index, options.kmer_index_name);

  std::cerr <<"..done.\n";
  // check for gzip compression

  std::string filenameExtension1 = options.readfile1.substr(options.readfile1.rfind('.')+1);
  std::string filenameExtension2 = options.readfile2.substr(options.readfile2.rfind('.')+1);
  if(filenameExtension1 == "gz" && filenameExtension2 == "gz"){
    if(mapLinkedZipped(Index, options, argc, argv)){
      return 1;
    }
    return 0;
  }
  if(filenameExtension1 == "gz" || filenameExtension2 == "gz"){
    std::cerr << "\nError!: Both readfiles must use same level of compression.\n";
    return 1;
  }

  // defining Parameters
  uint8_t window_count=50;   // amount of saved candidate windows

  //checking file size, barcode size and readInPairSyntax (is there /1 and /2 present)
  SeqFileIn file1(toCString(options.readfile1));
  SeqFileIn file2(toCString(options.readfile2));

  CharString (*getIdFunction)(std::string); // function pointer
  {
    std::string line;
    getline(file1.stream.file, line, ' ');
    if (line.substr(line.size()-2,line.size()) == "/1" || line.substr(line.size()-2,line.size()) == "/2"){
      getIdFunction = &getPairedID;
    }else{
      getIdFunction = &getID;
    }
  }

  file1.stream.file.seekg(0, std::ios::end);
  file2.stream.file.seekg(0, std::ios::end);
  std::streampos readfile1_size=file1.stream.file.tellg();
  std::streampos readfile2_size=file2.stream.file.tellg();
  close(file1);
  close(file2);
  uint8_t barcode_length=getBarcodeLength(options.readfile1, readfile1_size);
  if(barcode_length == 0){
    return 0;
  }

  // create folder for readfile index
  if (mkdir(toCString(options.read_index_name), 0777) == -1){
    std::cerr << "Error for readfile index target:  " << strerror(errno) << "\n";
  }

  if(checkReadfile(options.readfile1)){return 1;};
  if(checkReadfile(options.readfile2)){return 1;};

  std::cerr << "Processing read file...";

  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  omp_lock_t lock;
  omp_lock_t reslock;
  omp_init_lock(&lock);
  omp_init_lock(&reslock);

  std::vector<uint32_t> histogram(200,0);
  std::vector<result_t> globalresults;

  omp_set_num_threads(options.threads);
  uint64_t random_seed = getRandSeed(options.k);

  #pragma omp parallel for ordered
  for (int t=0; t<options.threads; t++){
    //declare variables
    std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t>> kmer_list;   // (i,j,a,m_a)   i=reference (Chromosome), j=position of matching k-mer in reference, a=abundance of k-mer in reference, m_a=minimizer_active_bases
    minimizedSequence miniSeq(options.k,options.mini_window_size,random_seed);
    StringSet<Dna5String> reads;
    resize(reads, 2, Exact());
    Dna5String read1;
    Dna5String read2;
    CharString id1;
    CharString id2;
    std::string barcode;
    std::string new_barcode;
    std::vector<result_t> results;
    // std::string results;
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
      while(line[0]!='+'){
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
        if(SearchID(file2, getIdFunction(toCString(id1)),startpos1, endpos1, getIdFunction)){
          break;
        }
        if(SearchID(file2, getIdFunction(toCString(id1)),startpos2, endpos2, getIdFunction)){
          break;
        }
        endpos1=startpos1;
        if(startpos1>100000){startpos1-=100000;}else{startpos1=0;}
        startpos2=endpos2;
        if(readfile2_size-endpos2>100000){endpos2+=100000;}else{endpos2=readfile2_size;}
      }
    }

    // skip to first valid barcode
    while (barcode[0]=='*' && !atEnd(file1)){
      skipToNextBarcode2(file1,file2,barcode,barcode_length);
    }

    BCI_1s=file1.stream.file.tellg();
    BCI_2s=file2.stream.file.tellg();

    // if skiped beyond the boundaries of threads scope: end thread
    if (file1.stream.file.tellg()>endpos){
      continue;
    }
    //proceed through readfile untill endpos
    int32_t readCount=0;
    while (!atEnd(file1)) { // proceeding through files
      pos_temp=file1.stream.file.tellg();
      readRecord(id1, read1, file1);
      new_barcode=getBarcode(toCString(id1),barcode_length);

      if (barcode!=new_barcode){ //If Barcode changes: map kmer_list and reinitialize kmer_list
        readCount=0;
        // std::cerr << "Mapping barcode: " << barcode << "\n";
        //append Barcode Index
        BCI_local.push_back(std::make_tuple(barcode, BCI_1s, BCI_2s));
        BCI_1s=pos_temp;

        // map barcode and clear k_mer list
        if (!kmer_list.empty()) {
          std::sort(kmer_list.begin(),kmer_list.end());
          MapKmerListLinked(kmer_list, window_count, barcode, results, Index.lookChrom, histogram_local, options);
          kmer_list.clear();
          if (options.Sort==0 && results.size()>1000){
            if (omp_test_lock(&lock)){
              std::fstream output;
              output.open(options.output_file,std::ios::out | std::ios::app);
              for (std::vector<result_t>::iterator itr_res=results.begin(); itr_res < results.end(); itr_res++){
                output << (*itr_res).string();
              }
              results.clear();
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
        // std::cerr << "Barcode: " << barcode << " mapped!" << "\t collecting reads of barcode: "<< new_barcode << "\n";
      }

      //process next read
      readRecord(id2, read2, file2);
      assignValue(reads,0,read1);
      assignValue(reads,1,read2);

      barcode=new_barcode;
      for (TStringSetIterator it = begin(reads); it!=end(reads); ++it){  // Iterating over the two reads of the read pair
        // initialize minimizer class
        minimizer mini;
        miniSeq.init(*it);
        // return and save all minimizers to kmer_list
        while(!miniSeq.at_end){
          mini=miniSeq.pop();
          AppendPos(kmer_list, mini.value^miniSeq.random_seed, Index, mini.active_bases, options);
        }
      }
      readCount++;
      if(readCount>=50000){
        kmer_list.clear();
        readCount=0;
      }

    }

    if (!kmer_list.empty()) { // only ever happens for the last thread on the last barcode
      // BCI_local.push_back(std::make_tuple(barcode, BCI_1s, readfile1_size, BCI_2s, readfile2_size));
      std::sort(kmer_list.begin(),kmer_list.end());
      MapKmerListLinked(kmer_list, window_count,barcode, results, Index.lookChrom, histogram_local, options);
    }

    close(file1);
    close(file2);

    // writing readfile index to files and summing up local histograms
    #pragma omp ordered
    {
      std::string readindex_file_id;
      readindex_file_id = (std::get<0>(BCI_local[0])).substr(0,2); // the first 2 characters of a barcode determine in which readindex file it is saved

      std::ofstream file_bci;
      file_bci.open(options.read_index_name + "/" + readindex_file_id , std::ios::app);

      for (int i=0; i<BCI_local.size(); i++){
        if(readindex_file_id != (std::get<0>(BCI_local[i])).substr(0,2)){
          readindex_file_id = (std::get<0>(BCI_local[i])).substr(0,2); // the first 2 characters of a barcode determine in which readindex file it is saved
          file_bci.close();
          file_bci.open(options.read_index_name + "/" + readindex_file_id , std::ios::app);
        }
        file_bci  << std::get<0>(BCI_local[i]) << "\t"
                  << std::get<1>(BCI_local[i]) << "\t"
                  << std::get<2>(BCI_local[i]) << "\n";
      }
      file_bci.close();


      // summing up local histograms
      for (int i=0; i<histogram.size(); i++){
        histogram[i]+=histogram_local[i];
      }
    }

    std::fstream output;

    if (options.Sort==0){ // write final output to file
      omp_set_lock(&lock);
      output.open(options.output_file,std::ios::out | std::ios::app);
      for (std::vector<result_t>::iterator itr_res=results.begin(); itr_res < results.end(); itr_res++){
        output << (*itr_res).string();
      }
      output.close();
      omp_unset_lock(&lock);

    }else{ // sort output and then merge into global output
      sortResults(results);
      omp_set_lock(&reslock);
      uint64_t middle = globalresults.size();
      globalresults.insert(globalresults.end(),results.begin(),results.end());
      std::inplace_merge(globalresults.begin(), globalresults.begin()+middle, globalresults.end(), compFunctionResult);
      omp_unset_lock(&reslock);
    }
    results.clear();
  } // parallel loop

  if(options.CoverageAnalysis==1){ // perform coverage analysis
    coverageAnalysis(globalresults, histogram, options);
  }

  if(options.Sort==1){ // write global results to file
    omp_set_lock(&lock);
    std::fstream output;
    output.open(options.output_file,std::ios::out | std::ios::app);
    for (std::vector<result_t>::iterator itr_res=globalresults.begin(); itr_res < globalresults.end(); itr_res++){
      output << (*itr_res).string();
    }
    output.close();
    omp_unset_lock(&lock);
  }

  // write last entrys of readfile index
  writeLastReadFileIndexEntrys(readfile1_size, readfile2_size, options);

  // write histogram to file
  writeHistogram(histogram, options.output_file);

  std::cerr << ".........done.\n";

  std::cerr << "Barcodes mapped sucessfully!\n\n";

  return 0;
} //map(argc,argv)

//skips file1 to start of next barcode and returns the barcode
std::string skipToNextBarcode(SeqFileIn & file, CharString & id1, uint8_t barcode_length){
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
void skipToNextBarcode2(SeqFileIn & file1, SeqFileIn & file2, std::string & barcode, uint8_t barcode_length){
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
bool SearchID(SeqFileIn & file, CharString id, std::streampos startpos, std::streampos endpos, CharString (*getIdFunction)(std::string)){
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
    new_id=getIdFunction(toCString(new_id));
  }
  file.stream.file.seekg(pos);
  return true;
}

// extract first entrys of readfile index subfiles and append to previous readfile index subfile
void writeLastReadFileIndexEntrys(std::streampos & readfile1_size, std::streampos & readfile2_size, mapOptions & options){
  std::ifstream file_bci_in;
  std::ofstream file_bci_out;
  std::string line; // store first line of file to append to previous file
  std::vector<std::string> files{"TG","TC","TA","GT","GG","GC","GA","CT","CG","CC","CA","AT","AG","AC"}; // TT and AA are handeled seperately

  // last entry of last file
  file_bci_out.open(options.read_index_name + "/TT" , std::ios::app);
  file_bci_out << "ZZZZZZZZZZZZZZZZ\t" << readfile1_size << "\t" << readfile2_size; // write last line of file TT
  file_bci_out.close();
  file_bci_in.open(options.read_index_name + "/TT");
  std::getline(file_bci_in, line); // extract first line of file TT
  file_bci_in.close();

  for(std::vector<std::string>::iterator file = files.begin(); file != files.end(); file++){ // iterating over files
    file_bci_out.open(options.read_index_name + "/" + *file , std::ios::app);
    file_bci_out << line;
    file_bci_out.close();
    file_bci_in.open(options.read_index_name + "/" + *file);
    std::getline(file_bci_in, line);
    file_bci_in.close();
  }

  file_bci_out.open(options.read_index_name + "/AA" , std::ios::app); // write last line of file AA
  file_bci_out << line;
  file_bci_out.close();

  return;
}
