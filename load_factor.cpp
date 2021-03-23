# include <iostream>
# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <fstream>
# include <time.h>
using namespace seqan;



int main(int argc, char *argv[]){

if(argc!=2){
  std::cerr << "Usage: ./bcmap readFile1 readFile2 k Index_name mini_window_size\n\n";
  exit(-1);
}

std::string IndC=argv[1];
IndC.append("_C.txt");
String<int long long> C;

String<int long long, External<> > extC;
if (!open(extC, IndC.c_str(), OPEN_RDONLY)){
  throw std::runtime_error("Could not open index counts file." );
}
assign(C, extC, Exact());
close(extC);

typedef Iterator<String<int long long>>::Type CIterator;

uint64_t counter=0;

for (CIterator it=begin(C);it!=end(C);++it){
  if(*it!=-1){
    counter++;
  }
}

std::cerr << "\nload factor: " << (long double)counter/length(C) << "\n";
return 0;

}
