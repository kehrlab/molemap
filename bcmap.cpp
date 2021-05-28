# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <iostream>
# include <fstream>
# include "./src/functions.h"
# include "./src/index.h"
# include "./src/get.h"
# include "./src/map.h"
# include <time.h>
using namespace seqan;

int main(int argc, char const ** argv){
  std::string arg1(argv[1]);
  --argc;
  argv++;
  if (argc==0){
    std::cerr << "\nNo command!\nCommands are: index, map, get.\nUse './bcmap [command] --help' for more information.\n";
  }
  else if(arg1=="index"){
    index(argc, argv);
  }
  else if(arg1=="map"){
    map(argc, argv);
  }
  else if(arg1=="get"){
    get(argc, argv);
  }
  else{
    std::cerr << "\nInvalid command! " << arg1 << "\nCommands are: index, map, get.\nUse './bcmap [command] --help' for more information.\n";
  }

}
