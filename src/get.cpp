# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <iostream>
# include <fstream>
# include "functions.h"
# include "get.h"
# include <time.h>
# include "parser.h"
using namespace seqan;

void returnReads(std::vector<std::string> & BCI_BC, std::vector<std::tuple<std::streampos,std::streampos>> & BCI, std::set<std::string> & barcodes, std::ifstream & file1, std::ifstream & file2, getOptions & options, std::string & fileID);

void loadReadfileIndex(std::vector<std::tuple<std::streampos,std::streampos>> & BCI, std::vector<std::string> & BCI_BC, getOptions & options, std::string & fileID){

  BCI.clear();
  BCI_BC.clear();
  std::ifstream file_bci;
  std::string BCI_1s;
  std::string BCI_2s;
  std::string BCI_bc;
  file_bci.open(options.read_index_name + "/" + fileID);

  while (!file_bci.eof()){
    file_bci >> BCI_bc;
    file_bci >> BCI_1s;
    file_bci >> BCI_2s;
    BCI.push_back(std::make_tuple(std::stoll(BCI_1s), std::stoll(BCI_2s)));
    BCI_BC.push_back(BCI_bc);
  }
  file_bci.close();

  return;
}

int get(int argc, char const ** argv){

  // parsing comand line
  getOptions options;
  seqan::ArgumentParser::ParseResult res = parseCommandLine_get(options, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;
  printParseResults_get(options);

  std::cerr << "\nobtaining barcodes...........";

  std::ifstream file1;
  file1.open(options.readfile1);
  std::ifstream file2;
  file2.open(options.readfile2);

  std::set<std::string> barcodes;
  std::string barcode;

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

  // initializing the index in readfileIndex index
  std::vector<std::tuple<std::streampos,std::streampos>> BCI;
  std::vector<std::string> BCI_BC;
  // load readfileIndex from file
  std::string fileID = barcodes.begin()->substr(0,2);
  loadReadfileIndex(BCI, BCI_BC, options, fileID);

  end = omp_get_wtime();
  std::cerr << "done.";
  printf(" in %f seconds", end - start);


  std::cerr << "\nextract reads for barcodes...";
  start = omp_get_wtime();

  returnReads(BCI_BC, BCI, barcodes, file1, file2, options, fileID);

  end = omp_get_wtime();
  std::cerr << "done.";
  printf(" in %f seconds", end - start);

  file1.close();
  file2.close();

  return 0;
}

void writeResults(std::string & results1, std::string & results2, std::ofstream & output_file1, std::ofstream & output_file2){
  output_file1 << results1;
  output_file2 << results2;
  results1 = "";
  results2 = "";
  return;
}

void returnReads(std::vector<std::string> & BCI_BC, std::vector<std::tuple<std::streampos,std::streampos>> & BCI, std::set<std::string> & barcodes, std::ifstream & file1, std::ifstream & file2, getOptions & options, std::string & fileID){
  std::ofstream output_file1;
  output_file1.open(options.output_file+".1.fq");
  std::ofstream output_file2;
  output_file2.open(options.output_file+".2.fq");
  std::string results1 = "";
  std::string results2 = "";

  for (std::set<std::string>::iterator itrbc = barcodes.begin(); itrbc!=barcodes.end(); itrbc++){ // for every barcode
    std::cerr << *itrbc << "\n";
    if(itrbc->substr(0,2)!=fileID){ // check if the correct index is loaded
      fileID=itrbc->substr(0,2);
      loadReadfileIndex(BCI, BCI_BC, options, fileID); // load correct index
    }
    std::vector<std::string>::iterator itrpos = std::lower_bound(BCI_BC.begin(), BCI_BC.end(), *itrbc);
    if (itrpos!=BCI_BC.end()){
      uint32_t pos = std::distance(BCI_BC.begin(), itrpos);
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
      if (results1.size() > 100000){
        writeResults(results1, results2, output_file1, output_file2); // write current results to output
      }
    }else{
      std::cerr << "\nBarcode " << *itrbc << " not present in readfile index.\n";
    }
  }

  writeResults(results1, results2, output_file1, output_file2); // write current results to output
  output_file1.close();
  output_file2.close();
  return;
}
