# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <iostream>
# include <fstream>
# include "functions.h"
# include "parser.h"
# include <seqan/misc/interval_tree.h>

using namespace seqan;


seqan::ArgumentParser::ParseResult parseCommandLine_index(indexOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("molemap index");

    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "reference(.fastq/.fasta)"));

    // Define Options
    addOption(parser, seqan::ArgParseOption(
        "o", "kmer_index_name", "Name of the folder in which the kmer index is stored.",
        seqan::ArgParseArgument::STRING, "kmer_index_name[OUT]"));
    setDefaultValue(parser, "o", "Index");
    addOption(parser, seqan::ArgParseOption(
        "p", "preset", "Sets default parameters based on sequencing technology. Choose 'off' to manualy set parametes.",
        seqan::ArgParseArgument::STRING, "kmer_index_name[OUT]"));
    setDefaultValue(parser, "p", "off");
    setValidValues(parser, "p", "linked long off");

    addOption(parser, seqan::ArgParseOption(
        "k", "kmer_length", "Length of kmers in index.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "k", "31");
    setMinValue(parser, "k", "15");
    setMaxValue(parser, "k", "31");
    addOption(parser, seqan::ArgParseOption(
        "m", "minimizer_window", "Length of window a minimizer is chosen from.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "m", "61");

    seqan::addUsageLine(parser,"reference.fq [OPTIONS]");
    setShortDescription(parser, "Build an index of a reference genome.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addDescription(parser,"Builds an open adressing k-mer index for the given reference genome(fastq/fasta).");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine_index()
    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;}

    // Extract option and argument values.
    getOptionValue(options.k, parser, "k");
    getOptionValue(options.m, parser, "m");
    getOptionValue(options.preset, parser, "p");
    getOptionValue(options.kmer_index_name, parser, "o");
    getArgumentValue(options.reference_file, parser, 0);

    // set parameters based on preset
    if(options.preset=="linked"){
      options.k=31;
      options.m=61;
    }else if(options.preset=="long"){
      options.k=17;
      options.m=25;
    }

    return seqan::ArgumentParser::PARSE_OK;
}

void printParseResults_index(indexOptions & options){
    std::cerr <<'\n'
              << "reference        \t" << options.reference_file << '\n'
              << "kmer_index_name  \t" << options.kmer_index_name << '\n'
              << "preset           \t" << options.preset << '\n'
              << "k                \t" << options.k << '\n'
              << "minimizer_window \t" << options.m << "\n\n";
    return;
}

seqan::ArgumentParser::ParseResult parseCommandLine_map(mapOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("molemap mapLinked");

    // Define arguments.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "readfile1.fastq"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "readfile2.fastq"));
    // addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Index_name[IN]"));
    // addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Barcode_index_name[OUT]"));

    // Define Options
    addOption(parser, seqan::ArgParseOption(
        "i", "kmer_index_name", "Name of the folder in which the kmer index is stored.",
        seqan::ArgParseArgument::STRING, "kmer_index_name[IN]"));
    setDefaultValue(parser, "i", "Index");
    addOption(parser, seqan::ArgParseOption(
        "r", "Read_index_name", "Name of the folder in which the ReadIndex is stored.",
        seqan::ArgParseArgument::STRING, "Index_name[IN]"));
    setDefaultValue(parser, "r", "ReadIndex");
    addOption(parser, seqan::ArgParseOption(
        "a", "max_abundance", "Maximum abundance of minimizer in reference to be considered.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "a", "20");
    addOption(parser, seqan::ArgParseOption(
        "w", "max_window_size", "Maximum length of genomic windows.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "w", "300000");
    addOption(parser, seqan::ArgParseOption(
        "g", "max_gap_size", "Maximum gap between minimizer hits of same genomic window.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "g", "20000");
    addOption(parser, seqan::ArgParseOption(
        "o", "output", "Path to the output file.",
        seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "o", "barcode_windows.bed");
    addOption(parser, seqan::ArgParseOption(
        "s", "score_threshold", "Minimum score threshold for genomic windows.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "s", "0");
    addOption(parser, seqan::ArgParseOption(
        "l", "length", "Length threshold for genomic windows.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "l", "10000");
    setMinValue(parser, "l", "5000");
    addOption(parser, seqan::ArgParseOption(
        "t", "threads", "Number of threads available.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "t", "16");
    addOption(parser, seqan::ArgParseOption(
        "S", "sort", "Sort barcode mappings by position"));
    addOption(parser, seqan::ArgParseOption(
        "C", "coverage_analysis", "Perform coverage analysis to filter output. Only recommended for WGS data with decent coverage. (Acitvates sorting by position)"));

    seqan::addUsageLine(parser,"readfile.1.fq readfile.2.fq [OPTIONS]");
    setShortDescription(parser, "Map barcodes to reference.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addDescription(parser,
               "Barcodes will be mapped to reference genome. "
               "Returns genomic windows from which barcoded reads most likely originate. "
               "Each window is rated by a quality score. "
               "Requires readfiles to be sorted by barcode (use bcctools). "
               "Requires reference to be indexed using the 'index' command. ");
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract argument and option values.
    getArgumentValue(options.readfile1, parser, 0);
    getArgumentValue(options.readfile2, parser, 1);

    getOptionValue(options.kmer_index_name, parser, "i");
    getOptionValue(options.read_index_name, parser, "r");
    getOptionValue(options.max_abundance, parser, "a");
    getOptionValue(options.max_window_size, parser, "w");
    getOptionValue(options.max_gap_size, parser, "g");
    getOptionValue(options.output_file, parser, "o");
    getOptionValue(options.s, parser, "s");
    getOptionValue(options.l, parser, "l");
    getOptionValue(options.threads, parser, "t");
    options.Sort = isSet(parser, "S") || isSet(parser, "C");
    options.CoverageAnalysis = isSet(parser, "C");

    loadIndexParameters(options.k,options.mini_window_size,options.kmer_index_name);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    return seqan::ArgumentParser::PARSE_OK;
}

void printParseResults_map(mapOptions & options){
  std::cerr <<'\n'
  << "readfile1         \t" << options.readfile1 << '\n'
  << "readfile2         \t" << options.readfile2 << '\n'
  << "kmer_index_name   \t" << options.kmer_index_name << '\n'
  << "  k - " << options.k << '\n'
  << "  w - " << options.mini_window_size << '\n'
  << "output file       \t" << options.output_file << '\n';
  if(options.readfile1.substr(options.readfile1.rfind('.')+1)!="gz"){
    std::cerr << "read_index_name   \t" << options.read_index_name << '\n';
  }
  std::cerr << "threads           \t" << options.threads << '\n'
  << "score threshold   \t" << options.s << '\n'
  << "max kmer abundance\t" << options.max_abundance << '\n'
  << "max window size   \t" << options.max_window_size << '\n'
  << "max gap size      \t" << options.max_gap_size << '\n'
  << "length threshold  \t" << options.l << '\n'
  << "Sort by position  \t" << (options.Sort ? "true" : "false") << '\n'
  << "Coverage Analysis \t" << (options.CoverageAnalysis ? "true" : "false") << "\n\n";
  return;
}

int getReadGroupId(std::string & readGroup, std::string & readGroupId){
  std::string readGroupField;
  std::istringstream iss(readGroup);
  std::string newReadGroup="";

  std::getline(iss,readGroupField,' ');
  if(readGroupField!="@RG"){
    std::cerr << "ERROR: Readgroup header line must start with '@RG'. Please provide the Readgroup in the following format: '@RG ID:YOUR_ID ...'\n";
    std::cerr << "The provided Readgroup option will be ignored.\n";
    return 1;
  }
  newReadGroup+=readGroupField;
  std::getline(iss,readGroupField,' ');
  if(readGroupField.substr(0,3)!="ID:"){
    std::cerr << "ERROR: First entry of Readgroup header line must be 'ID'. Please provide the Readgroup in the following format: '@RG ID:YOUR_ID ...'\n";
    std::cerr << "The provided Readgroup option will be ignored.\n";
    return 1;
  }
  readGroupId=readGroupField.substr(3,std::string::npos);
  newReadGroup+="\t"+readGroupField;
  while(std::getline(iss,readGroupField,' ')){
    newReadGroup+="\t"+readGroupField;
  }

  readGroup = newReadGroup;
  return 0;
}

void getRegions(std::string regionString, longmapOptions & options){
  // chr1:123-234,chr2:412-516...
  if(regionString.substr(regionString.size()-4)==".bed"){
    // open bed-file
    std::ifstream bedFile(regionString);
    // iterate trhough bed file and insert regions
    std::string line;
    std::string chrom;
    uint32_t start;
    uint32_t end;
    std::string item;

    while(std::getline(bedFile, line)){
      std::stringstream lineS(line);
      lineS >> chrom;
      lineS >> start;
      lineS >> end;
      // if(!options.regions.count(chrom)){ // if no tree exists for this chrom: create one
      //   IntervalTree<uint32_t, bool> newTree;
      //   // options.regions[chrom]=newTree;
      //   options.regions.insert(chrom,newTree);
      // }
      seqan::addInterval(options.regions[chrom], start, end, true);
    }

    return;
  }

  region regionBuff;
  size_t pos;
  while(regionString.find(':') != std::string::npos){
    regionBuff.chrom = regionString.substr(0,regionString.find_first_of(':',0));
    regionBuff.start = std::stoi(regionString.substr(regionString.find_first_of(':',0)+1,regionString.find_first_of('-',0)-regionString.find_first_of(':',0)-1));
    regionBuff.end = std::stoi(regionString.substr(regionString.find_first_of('-',0)+1,regionString.find_first_of(',',0)-regionString.find_first_of('-',0)-1));
    // if(!options.regions.count(regionBuff.chrom)){ // if no tree exists for this chrom: create one
    //   IntervalTree<int, bool> newTree;
    //   options.regions[regionBuff.chrom]=newTree;
    // }
    seqan::addInterval(options.regions[regionBuff.chrom], regionBuff.start, regionBuff.end, true);

    if(regionString.find(',') != std::string::npos){
      regionString = regionString.substr(regionString.find_first_of(',',0)+1,std::string::npos);
    }else{
      regionString="";
    }
  }
  return;
}

seqan::ArgumentParser::ParseResult parseCommandLine_long_map(longmapOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("molemap mapLong");

    // Define arguments.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "readfile.fastq"));

    // Define Options
    addOption(parser, seqan::ArgParseOption(
        "i", "kmer_index_name", "Name of the folder in which the kmer index is stored.",
        seqan::ArgParseArgument::STRING, "kmer_index_name[IN]"));
    setDefaultValue(parser, "i", "Index");
    addOption(parser, seqan::ArgParseOption(
        "a", "max_abundance", "Maximum abundance of minimizer in reference to be considered.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "a", "20");
    addOption(parser, seqan::ArgParseOption(
        "g", "max_gap_size", "Maximum gap between minimizer hits of same genomic window.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "g", "20000");
    addOption(parser, seqan::ArgParseOption(
        "o", "output", "Path to the output file.",
        seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "o", "stdout");
    addOption(parser, seqan::ArgParseOption(
        "O", "format", "Output format.",
        seqan::ArgParseArgument::STRING, "string"));
    setDefaultValue(parser, "O", "sam");
    setValidValues(parser, "O", "sam bam");
    addOption(parser, seqan::ArgParseOption(
        "s", "score_threshold", "Minimum score threshold for genomic windows.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "s", "0");
    addOption(parser, seqan::ArgParseOption(
        "l", "min_read_len", "Minimum read length.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "l", "100");
    setMinValue(parser, "l", "0");
    addOption(parser, seqan::ArgParseOption(
        "t", "threads", "Number of threads available.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "t", "16");
    addOption(parser, seqan::ArgParseOption(
        "R", "Readgroup", "Read group header line.",
        seqan::ArgParseArgument::STRING, "string"));
    setDefaultValue(parser, "R", "");
    addOption(parser, seqan::ArgParseOption(
        "r", "Regions", "Regions to map to. Path to bed file or comma seperated list of chr:start-end",
        seqan::ArgParseArgument::STRING, "string"));
    setDefaultValue(parser, "r", "");

    seqan::addUsageLine(parser,"readfile.fq [OPTIONS]");
    setShortDescription(parser, "Map long reads to reference.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addDescription(parser,
               "Long-reads will be mapped to reference genome. "
               "Returns genomic windows from which reads most likely originate. "
               "Each window is rated by a quality score. "
               "Requires reference to be indexed using the 'index' command. ");
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract argument and option values.
    getArgumentValue(options.readfile1, parser, 0);

    getOptionValue(options.kmer_index_name, parser, "i");
    getOptionValue(options.max_abundance, parser, "a");
    getOptionValue(options.max_gap_size, parser, "g");
    getOptionValue(options.output_file, parser, "o");
    getOptionValue(options.output_format, parser, "O");
    getOptionValue(options.s, parser, "s");
    getOptionValue(options.l, parser, "l");
    getOptionValue(options.threads, parser, "t");
    getOptionValue(options.readGroup, parser, "R");
    if(isSet(parser, "R")){
      if(getReadGroupId(options.readGroup, options.readGroupId)){
        options.readGroup = "";
      }
    }
    std::string regionString;
    getOptionValue(regionString, parser, "r");
    if(isSet(parser, "r")){
      options.regionDefined=true;
      getRegions(regionString, options);
    }

    loadIndexParameters(options.k,options.mini_window_size,options.kmer_index_name);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    return seqan::ArgumentParser::PARSE_OK;
}

void printParseResults_long_map(longmapOptions & options){
    std::cerr <<'\n'
              << "readfile          \t" << options.readfile1 << '\n'
              << "kmer_index_name   \t" << options.kmer_index_name << '\n'
              << "  k - " << options.k << '\n'
              << "  w - " << options.mini_window_size << '\n'
              << "output file       \t" << options.output_file << '\n'
              << "output format     \t" << options.output_format << '\n'
              << "threads           \t" << options.threads << '\n'
              << "max kmer abundance\t" << options.max_abundance << '\n'
              << "score threshold   \t" << options.s << '\n'
              << "min read size     \t" << options.l << '\n'
              << "max gap size      \t" << options.max_gap_size << "\n\n";
    return;
}

seqan::ArgumentParser::ParseResult parseCommandLine_get(getOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("molemap get");

    // Define arguments.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "readfile1.fastq"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "readfile2.fastq"));
    // addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Barcode_index"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "Barcodes"));

    addOption(parser, seqan::ArgParseOption(
        "o", "output", "Prefix of the output files.",
        seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    setDefaultValue(parser, "o", "molemapGetOut");
    addOption(parser, seqan::ArgParseOption(
        "r", "read_index_name", "Name of the folder in which the ReadIndex is stored.",
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

void printParseResults_get(getOptions & options){
  std::cerr <<'\n'
            << "readfile1        \t" << options.readfile1 << '\n'
            << "readfile2        \t" << options.readfile2 << '\n'
            << "read_index_name  \t" << options.read_index_name << '\n'
            << "barcodes         \t" << options.barcodes << '\n'
            << "output prefix    \t" << options.output_file << "\n\n";
  return;
}

void loadIndexParameters(uint32_t & k, uint32_t & m, std::string & IndexName){
  std::string IndexInfoFile=IndexName+"/Info.txt";
  std::ifstream IndexInfo(IndexInfoFile);
  IndexInfo >> k >> m;
  IndexInfo.close();
  return;
}
