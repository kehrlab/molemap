# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <iostream>
# include <fstream>
# include "index.h"
# include "get.h"
# include "map_linked.h"
# include "map_long.h"
# include <time.h>
using namespace seqan;

int main(int argc, char const ** argv){

  if (argc==1){
    std::cerr << "\nNo command:\nCommands are: index, mapLong, mapLinked, get.\nUse './molemap [command] --help' for more information.\n";
    return 1;
  }

  std::string arg1(argv[1]);
  --argc;
  argv++;

  if(arg1=="index"){
    return index(argc, argv);
  }
  if(arg1=="maplinked" || arg1=="mapLinked"){
    return maplinked(argc, argv);
  }
  if(arg1=="maplong" || arg1=="mapLong"){
    return maplong(argc, argv);
  }
  if(arg1=="get"){
    return get(argc, argv);
  }
  else{
    std::cerr << "\nInvalid command: " << arg1 << "\nCommands are: index, mapLong, mapLinked, get.\nUse './molemap [command] --help' for more information.\n";
    return 1;
  }

  return 1;
}
