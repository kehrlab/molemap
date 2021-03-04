# include <iostream>
# include <fstream>
# include <seqan/seq_io.h>
# include <seqan/sequence.h>
using namespace std;

/*
g++ ./build_bcindex.cpp -o bcindex
*/

int main(int argc, char *argv[]){

  if(argc!=4){
    std::cerr << "Usage: ./bcindex readfile1.fastq readfile2.fastq Index_name\n\n";
    exit(-1);
  }

  std::string filename1=argv[1];
  std::string filename2=argv[2];
  std::string new_barcode;
  std::string old_barcode;
  std::string id;

  std::fstream file1;
  file1.open(filename1,std::ios::in);
  std::fstream file2;
  file2.open(filename2,std::ios::in);

  getline(file1,id);
  old_barcode=id.substr(id.find("RX:Z:")+5,16);
  std::cerr << "type: " << typeid(file1.tellg()).name() << "value: " << file1.tellg() << "\n";
  std::cerr << "barcode: " << new_barcode << "\n";

  while (!file1.eof()) {
    getline(file1,id);
    std::cerr<< id << "\n";
    new_barcode=id.substr(id.find("RX:Z:")+5,16);
    if (new_barcode!=old_barcode) {
      // SAVE POSITION
      std::cerr << "type: " << typeid(file1.tellg()).name() << "value: " << file1.tellg() << "\n";
      std::cerr << "barcode: " << new_barcode << "\n";
      old_barcode=new_barcode;
    }
    file1.ignore(100000,'\n');
    file1.ignore(100000,'\n');
    file1.ignore(100000,'\n');


  }

  // SAVE END POSITION

  // String<unsigned long long> dir;
  // resize(dir,bucket_number+1,0);

  file1.close();
  file2.close();
}
