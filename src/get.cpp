# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <iostream>
# include <fstream>
# include "functions.h"
# include "get.h"
# include <time.h>
# include <omp.h>
using namespace seqan;

void returnReads(std::vector<std::string> & BCI_BC, std::vector<std::tuple<std::streampos,std::streampos>> & BCI, std::set<std::string> & barcodes, std::ifstream & file1, std::ifstream & file2, std::string & results1, std::string & results2);

struct getReadsOptions{
  std::string readfile1;
  std::string readfile2;
  std::string read_index_name;
  std::string barcodes;
  std::string output_file;

  getReadsOptions() :
  output_file("bcmapGetOut"), read_index_name("BarcodeIndex")
  {}
};

seqan::ArgumentParser::ParseResult parseCommandLine(getReadsOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bcmap get");

    // Define arguments.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "readfile1.fastq"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "readfile2.fastq"));
    // addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Barcode_index"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Barcodes"));

    addOption(parser, seqan::ArgParseOption(
        "o", "output", "Prefix of the output files.",
        seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "o", "bcmapGetOut");
    addOption(parser, seqan::ArgParseOption(
        "r", "read_index_name", "Name of the ReadIndex.",
        seqan::ArgParseArgument::STRING, "read_index_name[IN]"));
    setDefaultValue(parser, "r", "ReadIndex");

    seqan::addUsageLine(parser, "readfile.1.fq readfile.2.fq barcodes [OPTIONS]");
    setShortDescription(parser, "Retreive all reads of a list of barcodes.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addDescription(parser,
               "Retreives all reads belonging to the given set of barcodes. "
               "The reads are quickly extracted from the readfiles using a barcode index. "
               "Barcodes can be provided in a newline seperated textfile or ',' seperated as argument.");
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;}

    // Extract argument and option values.
    getArgumentValue(options.readfile1, parser, 0);
    getArgumentValue(options.readfile2, parser, 1);
    // getArgumentValue(options.bci_name, parser, 2);
    getArgumentValue(options.barcodes, parser, 2);

    getOptionValue(options.read_index_name, parser, "r");
    getOptionValue(options.output_file, parser, "o");

    return seqan::ArgumentParser::PARSE_OK;
}


int get(int argc, char const ** argv){

  getReadsOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;
  std::cout <<'\n'
            << "readfile1        \t" << options.readfile1 << '\n'
            << "readfile2        \t" << options.readfile2 << '\n'
            << "readindex        \t" << options.read_index_name << '\n'
            << "barcodes         \t" << options.barcodes << '\n'
            << "output prefix    \t" << options.output_file << "\n\n";


  std::ifstream file1;
  file1.open(options.readfile1);
  std::ifstream file2;
  file2.open(options.readfile2);

  std::set<std::string> barcodes;
  std::string barcode;


  std::cerr << "\nobtaining barcodes...........";

  double start = omp_get_wtime();

  // get barcodes
  std::ifstream file_bc;
  file_bc.open(options.barcodes /*, std::ios::binary*/);
  if(file_bc.is_open()){      // if barcodes provided as file: read bacodes from file
    while(getline(file_bc, barcode)){
      if (barcode!=""){
        barcodes.insert(barcode);
      }
    }
    file_bc.close();
  } else {                  // if barcodes provided as argument: read barcodes from argument
    std::size_t pos_s = 0;
    std::size_t pos = options.barcodes.find(",");
    while (pos!=std::string::npos/*<options.barcodes.end()*/){
      barcodes.insert(options.barcodes.substr(pos_s,pos-pos_s));
      pos_s=pos+1;
      pos=options.barcodes.find(",",pos_s);
    }
    barcodes.insert(options.barcodes.substr(pos_s,pos-pos_s));
  }

  double end = omp_get_wtime();
  std::cerr << "done.";
  printf(" in %f seconds", end - start);
  std::cerr << "\nloading read index...........";

  start = omp_get_wtime();

  // load in barcode index
  std::ifstream file_bci;
  std::vector<std::tuple<std::streampos,std::streampos>> BCI;
  std::vector<std::string> BCI_BC;
  std::string BCI_1s;
  std::string BCI_2s;
  std::string BCI_bc;
  file_bci.open(options.read_index_name /*, std::ios::binary*/);

  while (!file_bci.eof()){
    file_bci >> BCI_bc;
    file_bci >> BCI_1s;
    file_bci >> BCI_2s;
    BCI.push_back(std::make_tuple(std::stoll(BCI_1s), std::stoll(BCI_2s)));
    BCI_BC.push_back(BCI_bc);
  }
  file_bci.close();

  end = omp_get_wtime();
  std::cerr << "done.";
  printf(" in %f seconds", end - start);
  std::cerr << "\nextract reads for barcodes...";


  start = omp_get_wtime();

  std::string output1 ="";
  std::string output2 ="";

  returnReads(BCI_BC, BCI, barcodes, file1, file2, output1, output2);

  end = omp_get_wtime();

  std::cerr << "done.";
  printf(" in %f seconds", end - start);
  std::cerr << "\nwrite output to file........";

  start = omp_get_wtime();

  std::ofstream output_file1;
  output_file1.open(options.output_file+".1.fq");
  std::ofstream output_file2;
  output_file2.open(options.output_file+".2.fq");
  output_file1 << output1;
  output_file2 << output2;
  output_file1.close();
  output_file2.close();

  // std::cerr << "\nresults:\n";
  // for (int i = 0; i < results.size(); i++){
  //   std::cerr << results[i] << "\n";
  // }
  end = omp_get_wtime();
  std::cerr << "done.";
  printf(" in %f seconds\n\n", end - start);

  file1.close();
  file2.close();

  return 0;
}

void returnReads(std::vector<std::string> & BCI_BC, std::vector<std::tuple<std::streampos,std::streampos>> & BCI, std::set<std::string> & barcodes, std::ifstream & file1, std::ifstream & file2, std::string & results1, std::string & results2){

  for (std::set<std::string>::iterator itrbc = barcodes.begin(); itrbc<barcodes.end(); itrbc++){
    // std::cerr << *itrbc << "\n";
    std::vector<std::string>::iterator itrpos = std::lower_bound(BCI_BC.begin(), BCI_BC.end(), *itrbc);
    if (itrpos!=BCI_BC.end()){
      uint_fast32_t pos = std::distance(BCI_BC.begin(), itrpos);
      file1.seekg(std::get<0>(BCI[pos]));
      file2.seekg(std::get<1>(BCI[pos]));
      std::streamoff size1=std::get<0>(BCI[pos+1])-std::get<0>(BCI[pos]);
      std::streamoff size2=std::get<1>(BCI[pos+1])-std::get<1>(BCI[pos]);

      char * buff1 = new char [size1];
      char * buff2 = new char [size2];


      file1.read(buff1,size1);
      file2.read(buff2,size2);
      std::string stringbuff1(buff1);
      std::string stringbuff2(buff2);
      stringbuff1.resize(size1);
      stringbuff2.resize(size2);

      results1+=stringbuff1;
      results2+=stringbuff2;
    }else{
      std::cerr << "\nBarcode " << *itrbc << " not present in readfile index.\n";
    }
  }
  return;
}
