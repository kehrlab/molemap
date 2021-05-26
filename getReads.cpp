# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <iostream>
# include <fstream>
# include "./src/functions.h"
# include <time.h>
using namespace seqan;


struct getReadsOptions{
  std::string readfile1;
  std::string readfile2;
  std::string bci_name;
  std::string barcodes;
  std::string output_file;
  getReadsOptions() :
  output_file("getReadsOut.fastq")
  {}
};

seqan::ArgumentParser::ParseResult parseCommandLine(getReadsOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("getReads");

    // Define arguments.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "Path to readfile1.fastq"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "Path to readfile2.fastq"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Barcode_index_name"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Barcodes to search for"));

    addOption(parser, seqan::ArgParseOption(
        "o", "output", "Path to the output file.",
        seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "o", "getReadsOut.fastq");

    setShortDescription(parser, "Retreive all reads of a list of barcodes");
    setVersion(parser, "0.1");
    setDate(parser, "May 25 2021");
    addDescription(parser,
               "Retreives all reads belonging to the given set of barcodes."
               "The reads are quickly extracted from the readfiles using a barcode index"
               "Barcodes can be provided in a newline seperated textfile or ',' seperated as argument");
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;}

    // Extract argument and option values.
    getArgumentValue(options.readfile1, parser, 0);
    getArgumentValue(options.readfile2, parser, 1);
    getArgumentValue(options.bci_name, parser, 2);
    getArgumentValue(options.barcodes, parser, 3);

    getOptionValue(options.output_file, parser, "o");

    return seqan::ArgumentParser::PARSE_OK;
}


int main(int argc, char const ** argv){

  getReadsOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;
  std::cout <<'\n'
            << "readfile1        \t" << options.readfile1 << '\n'
            << "readfile2        \t" << options.readfile2 << '\n'
            << "barcodeindex_name\t" << options.bci_name << '\n'
            << "barcodes         \t" << options.barcodes << '\n'
            << "output file      \t" << options.output_file << "\n\n";


  try { // opening read-files
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

  std::vector<std::string> barcodes;
  std::string barcode;


  // get barcodes
  std::ifstream file_bc;
  file_bc.open(options.barcodes /*, std::ios::binary*/);
  if(file_bc.is_open()){      // if barcodes provided as file: read bacodes from file
    while(getline(file_bc, barcode)){
      barcodes.push_back(barcode);
    }
    file_bc.close();
  } else {                  // if barcodes provided as argument: read barcodes from argument
    std::size_t pos_s = 0;
    std::size_t pos = options.barcodes.find(",");
    while (pos!=std::string::npos/*<options.barcodes.end()*/){
      barcodes.push_back(options.barcodes.substr(pos_s,pos-pos_s));
      pos_s=pos+1;
      pos=options.barcodes.find(",",pos_s);
    }
    barcodes.push_back(options.barcodes.substr(pos_s,pos-pos_s));
  }

  // load in barcode index
  std::ifstream file_bci;
  std::vector<std::tuple<std::string,std::streampos,std::streampos,std::streampos,std::streampos>> BCI;
  std::streampos BCI_1s;
  std::streampos BCI_1e;
  std::streampos BCI_2s;
  std::streampos BCI_2e;
  std::string BCI_bc;
  file_bci.open(options.bci_name /*, std::ios::binary*/);
  while (!file_bci.eof()){
    BCI_bc << file_bci;
    BCI_1s << file_bci;
    BCI_1e << file_bci;
    BCI_2s << file_bci;
    BCI_2e << file_bci;
    BCI.push_back(make_tuple(BCI_bc, BCI_1s, BCI_1e, BCI_2s, BCI_2e))
  }
  file_bci.close();

  for (int i=0; i<BCI.size(); i++){
    std::cerr  << std::get<0>(BCI[i]) << "\t"
              << (int)std::get<1>(BCI[i]) << "\t"
              << (int)std::get<2>(BCI[i]) << "\t"
              << (int)std::get<3>(BCI[i]) << "\t"
              << (int)std::get<4>(BCI[i]) << "\n";
  }
  // std::cerr << "\nBarcodes:\n";
  // for (int i = 0; i < barcodes.size(); i++){
  //   std::cerr << barcodes[i] << "\n";
  // }

  close(file1);
  close(file2);
}
