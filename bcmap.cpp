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
  if(arg1=="index"){
    index(argc, argv);
  }
  else if(arg1=="map"){
    map(argc, argv);
  }
  else if(arg1=="get"){
    get(argc, argv);
  }
  else{
    std::cerr << "\ninvalid Command " << arg1 << ".\nCommands are: index, map, get.\nUse './bcmap [Command] --help' for more information.\n"
  }

}
