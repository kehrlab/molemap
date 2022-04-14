# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/arg_parse.h>
# include <iostream>
# include <fstream>
# include "functions.h"
# include "parser.h"
using namespace seqan;

seqan::ArgumentParser::ParseResult parseCommandLine(bcmapOptions & options, int argc, char const ** argv){
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("bcmap map");

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
        "r", "Read_index_name", "Name of the ReadIndex.",
        seqan::ArgParseArgument::STRING, "Index_name[IN]"));
    setDefaultValue(parser, "r", "ReadIndex");
    addOption(parser, seqan::ArgParseOption(
        "k", "kmer_length", "Length of kmers in index.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "k", "31");
    setMinValue(parser, "k", "8");
    setMaxValue(parser, "k", "31");
    addOption(parser, seqan::ArgParseOption(
        "m", "mini_window_size", "Length of minimizing window.",
        seqan::ArgParseArgument::INTEGER, "unsigned"));
    setDefaultValue(parser, "m", "61");
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
        "C", "coverage_analysis", "Perform coverage analysis to filter output. (Acitvates sorting by position)"));

    seqan::addUsageLine(parser,"readfile.1.fq readfile.2.fq [OPTIONS]");
    setShortDescription(parser, "Map barcodes to reference.");
    setVersion(parser, VERSION);
    setDate(parser, DATE);
    addDescription(parser,
               "Barcodes will be mapped to reference genome. "
               "Returns genomic windows from which barcoded reads most likely originate. "
               "Each window is rated by a quality score. "
               "Requires readfiles to be sorted by barcode (use bcctools). "
               "Requires reference to be indexed using the 'map' command. ");
    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract argument and option values.
    getArgumentValue(options.readfile1, parser, 0);
    getArgumentValue(options.readfile2, parser, 1);
    // getArgumentValue(options.index_name, parser, 2);
    // getArgumentValue(options.bci_name, parser, 3);

    getOptionValue(options.kmer_index_name, parser, "i");
    getOptionValue(options.read_index_name, parser, "r");
    getOptionValue(options.k, parser, "k");
    getOptionValue(options.mini_window_size, parser, "m");
    getOptionValue(options.max_window_size, parser, "w");
    getOptionValue(options.max_gap_size, parser, "g");
    getOptionValue(options.output_file, parser, "o");
    getOptionValue(options.s, parser, "s");
    getOptionValue(options.l, parser, "l");
    getOptionValue(options.threads, parser, "t");
    options.Sort = isSet(parser, "S") || isSet(parser, "C");
    options.CoverageAnalysis = isSet(parser, "C");

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    return seqan::ArgumentParser::PARSE_OK;
}

void printParseResults(bcmapOptions options){
    std::cout <<'\n'
              << "readfile1         \t" << options.readfile1 << '\n'
              << "readfile2         \t" << options.readfile2 << '\n'
              << "kmer_index_name   \t" << options.kmer_index_name << '\n'
              << "output file       \t" << options.output_file << '\n'
              << "read_index_name   \t" << options.read_index_name << '\n'
              << "threads           \t" << options.threads << '\n'
              << "score threshold   \t" << options.s << '\n'
              << "k                 \t" << options.k << '\n'
              << "minimizer window  \t" << options.mini_window_size << '\n'
              << "max window size   \t" << options.max_window_size << '\n'
              << "max gap size      \t" << options.max_gap_size << '\n'
              << "length threshold  \t" << options.l << '\n'
              << "Sort by position  \t" << (options.Sort ? "true" : "false") << '\n'
              << "Coverage Analysis \t" << (options.CoverageAnalysis ? "true" : "false") << '\n';
    return;
}
