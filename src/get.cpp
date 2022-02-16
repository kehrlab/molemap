# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <iostream>
# include <fstream>
# include "functions.h"
# include "get.h"
# include <time.h>
using namespace seqan;

std::vector<std::string> returnReads(  std::vector<std::string> & BCI_BC, std::vector<std::tuple<std::streampos,std::streampos,std::streampos,std::streampos>> & BCI, std::vector<std::string> & barcodes, SeqFileIn & file1, SeqFileIn & file2);

struct getReadsOptions{
  std::string readfile1;
  std::string readfile2;
  std::string read_index_name;
  std::string barcodes;
  std::string output_file;

  getReadsOptions() :
  output_file("bcmapGetOut.fastq"), read_index_name("BarcodeIndex")
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
        "o", "output", "Path to the output file.",
        seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "o", "bcmapGetOut.fastq");
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
  std::vector<std::tuple<std::streampos,std::streampos,std::streampos,std::streampos>> BCI;
  std::vector<std::string> BCI_BC;
  std::string BCI_1s;
  std::string BCI_1e;
  std::string BCI_2s;
  std::string BCI_2e;
  std::string BCI_bc;
  file_bci.open(options.read_index_name /*, std::ios::binary*/);

  while (!file_bci.eof()){
    file_bci >> BCI_bc;
    file_bci >> BCI_1s;
    file_bci >> BCI_1e;
    file_bci >> BCI_2s;
    file_bci >> BCI_2e;
    BCI.push_back(std::make_tuple(std::stoll(BCI_1s), std::stoll(BCI_1e), std::stoll(BCI_2s), std::stoll(BCI_2e)));
    BCI_BC.push_back(BCI_bc);
  }
  file_bci.close();


  // lookup barcodes

  std::vector<std::string> results;
  results=returnReads(BCI_BC, BCI, barcodes, file1, file2);

  std::ofstream output;
  output.open(options.output_file);
  for (int i = 0; i < results.size(); i++){
    output << results[i] << "\n";
  }
  output.close();

  // std::cerr << "\nresults:\n";
  // for (int i = 0; i < results.size(); i++){
  //   std::cerr << results[i] << "\n";
  // }

  close(file1);
  close(file2);

  return 0;
}

std::vector<std::string> returnReads(  std::vector<std::string> & BCI_BC, std::vector<std::tuple<std::streampos,std::streampos,std::streampos,std::streampos>> & BCI, std::vector<std::string> & barcodes, SeqFileIn & file1, SeqFileIn & file2){
  std::vector<std::string> result;
  std::string read;
  std::string id;
  std::string qual;

  for (std::vector<std::string>::iterator itrbc=barcodes.begin(); itrbc<barcodes.end(); itrbc++){
    std::vector<std::string>::iterator itrpos=std::lower_bound(BCI_BC.begin(), BCI_BC.end(), *itrbc);
    if (itrpos!=BCI_BC.end()){
      uint_fast32_t pos = std::distance(BCI_BC.begin(), itrpos);
      file1.stream.file.seekg(std::get<0>(BCI[pos]));
      file2.stream.file.seekg(std::get<2>(BCI[pos]));

      while(file1.stream.file.tellg() < std::get<1>(BCI[pos])){
        readRecord(id, read, qual, file1);
        // result.push_back("\n");
        result.push_back("@"+id);
        // result.push_back("\n");
        result.push_back(read);
        result.push_back("+");
        result.push_back(qual);
        readRecord(id, read, qual, file2);
        // result.push_back("\n");
        result.push_back("@"+id);
        // result.push_back("\n");
        result.push_back(read);
        result.push_back("+");
        result.push_back(qual);      }
    }
  }
  return result;
}
