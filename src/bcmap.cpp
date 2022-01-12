# include "../seqan/seq_io.h"
# include "../seqan/sequence.h"
# include "../seqan/arg_parse.h"
# include <iostream>
# include <fstream>
# include ".functions.h"
# include ".index.h"
# include ".get.h"
# include ".map.h"
# include <time.h>
using namespace seqan;

int main(int argc, char const ** argv){

  if (argc==1){
    std::cerr << "\nNo command!\nCommands are: index, map, get.\nUse './bcmap [command] --help' for more information.\n";
    return 0;
  }

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
    std::cerr << "\nInvalid command! " << arg1 << "\nCommands are: index, map, get.\nUse './bcmap [command] --help' for more information.\n";
  }

}
