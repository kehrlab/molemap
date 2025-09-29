# include "parser.h"
# include <seqan/sequence.h>
# include "seqan/seq_io.h"
# include "deconvolution.h"
using namespace seqan;


void filterReads(StringSet<Dna5String> & reads1, StringSet<Dna5String> & reads2, StringSet<CharString> & ids1, StringSet<CharString> & ids2, StringSet<CharString> & quals1, StringSet<CharString> & quals2, coverageIndex32_t & coverageIndex, deconOptions & options);
void writeOutput(StringSet<Dna5String> & reads1, StringSet<Dna5String> & reads2, StringSet<CharString> & ids1, StringSet<CharString> & ids2, StringSet<CharString> & quals1, StringSet<CharString> & quals2, deconOptions & options);
void fillIndex32FromReads(StringSet<Dna5String> & reads, coverageIndex32_t & coverageIndex, uint32_t k);
void fillCoverageIndex32(StringSet<Dna5String> & reads1, StringSet<Dna5String> & reads2, coverageIndex32_t & coverageIndex, uint32_t k);

int decon(int argc, char const **argv){

  // Argument parser #############################################################################
  deconOptions options;
  seqan::ArgumentParser::ParseResult res = parseDeconCommandLine(options, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK) return res;
  printDeconParseResult(options, res);

  // load readfiles
  std::cerr << "Loading readfiles...";

  StringSet<Dna5String> reads1;
  StringSet<Dna5String> reads2;
  StringSet<CharString> ids1;
  StringSet<CharString> ids2;
  StringSet<CharString> quals1;
  StringSet<CharString> quals2;

  SeqFileIn file1(toCString(options.readfile1_name));
  SeqFileIn file2(toCString(options.readfile2_name));

  readRecords(ids1, reads1, quals1, file1);
  readRecords(ids2, reads2, quals2, file2);

  close(file1);
  close(file2);

  std::cerr << "............done.\n";

  for(int loop=0; loop!=options.iterations; loop++){
    if(options.iterations!=1){
      std::cerr << "Iteration: " << loop << " of " << options.iterations << "\n";
    }
    std::cerr << "Building coverage index...";
    // double start=omp_get_wtime();

    coverageIndex32_t coverageIndex;
    uint64_t index_size = length(reads1)*(length(reads1[0])-options.k+1)*2.5;
    coverageIndex.resize(index_size);
    fillCoverageIndex32(reads1, reads2, coverageIndex, options.k);

    std::cerr << "......done.\n";

    // filter reads with low coverage (and write output)
    std::cerr << "Filtering reads...";

    filterReads(reads1, reads2, ids1, ids2, quals1, quals2, coverageIndex, options);

    std::cerr << "..............done.\n";
  }

  std::cerr << "Write Output...";

  writeOutput(reads1, reads2, ids1, ids2, quals1, quals2, options);

  std::cerr << "...............done.\n";

  return 0;
}

void writeOutput(StringSet<Dna5String> & reads1, StringSet<Dna5String> & reads2, StringSet<CharString> & ids1, StringSet<CharString> & ids2, StringSet<CharString> & quals1, StringSet<CharString> & quals2, deconOptions & options){
  typedef Iterator<StringSet<Dna5String> >::Type ReadSetIt_t;
  typedef Iterator<StringSet<CharString> >::Type IdSetIt_t;
  ReadSetIt_t itrReads1 = begin(reads1);
  ReadSetIt_t itrReads2 = begin(reads2);
  IdSetIt_t itrIds1 = begin(ids1);
  IdSetIt_t itrIds2 = begin(ids2);
  IdSetIt_t itrQuals1 = begin(quals1);
  IdSetIt_t itrQuals2 = begin(quals2);
  SeqFileOut outfile1(toCString(options.output_prefix+".1.fq"));
  SeqFileOut outfile2(toCString(options.output_prefix+".2.fq"));

  for(itrReads1; itrReads1!=end(reads1); itrReads1++){
    writeRecord(outfile1,*itrIds1,*itrReads1,*itrQuals1);
    writeRecord(outfile2,*itrIds2,*itrReads2,*itrQuals2);
    itrReads2++;
    itrIds1++;
    itrIds2++;
    itrQuals1++;
    itrQuals2++;
  }

  return;
}

// Functions for Barcode Deconvolution
void filterReads(StringSet<Dna5String> & reads1, StringSet<Dna5String> & reads2, StringSet<CharString> & ids1, StringSet<CharString> & ids2, StringSet<CharString> & quals1, StringSet<CharString> & quals2, coverageIndex32_t & coverageIndex, deconOptions & options){
  typedef Iterator<StringSet<Dna5String> >::Type ReadSetIt_t;
  typedef Iterator<StringSet<CharString> >::Type IdSetIt_t;
  uint32_t cov_buff;
  Iterator<Dna5String>::Type itrRead;

  ReadSetIt_t itrReads1 = begin(reads1);
  ReadSetIt_t itrReads2 = begin(reads2);
  IdSetIt_t itrIds1 = begin(ids1);
  IdSetIt_t itrIds2 = begin(ids2);
  IdSetIt_t itrQuals1 = begin(quals1);
  IdSetIt_t itrQuals2 = begin(quals2);

  StringSet<Dna5String> Reads1Buff, Reads2Buff;
  StringSet<CharString> Ids1Buff, Ids2Buff, Quals1Buff, Quals2Buff;

  // SeqFileOut outfile1(toCString(options.output_prefix+".1.fq"));
  // SeqFileOut outfile2(toCString(options.output_prefix+".2.fq"));

  for(itrReads1; itrReads1!=end(reads1); itrReads1++){

    uint32_t counter=0;
    // lookup first kmer
    Dna5String kmer = prefix(*itrReads1,options.k);
    hash32_t hash(kmer);
    uint64_t hashId = hash.bigger();
    coverageIndex.lookup(hashId, cov_buff);
    counter+=1*(cov_buff<=options.cutoff); // add to counter if coverage below cutoff

    //iterate over remaining kmers
    for(itrRead=begin(*itrReads1)+options.k;itrRead!=end(*itrReads1);itrRead++){
      hash.roll(*itrRead);
      hashId = hash.bigger();
      coverageIndex.lookup(hashId, cov_buff);
      counter+=1*(cov_buff<=options.cutoff); // add to counter if coverage below cutoff
    }

    kmer = prefix(*itrReads2,options.k);
    hash.rehash(kmer);
    hashId = hash.bigger();
    coverageIndex.lookup(hashId, cov_buff);
    counter+=1*(cov_buff<=options.cutoff); // add to counter if coverage below cutoff

    //iterate over remaining kmers
    for(itrRead=begin(*itrReads2)+options.k;itrRead!=end(*itrReads2);itrRead++){
      hash.roll(*itrRead);
      hashId = hash.bigger();
      coverageIndex.lookup(hashId, cov_buff);
      counter+=1*(cov_buff<=options.cutoff); // add to counter if coverage below cutoff
    }

    uint32_t kmer_count = length(*itrReads1)+length(*itrReads2)-options.k*2;
    uint32_t limit = (kmer_count)*options.cutoffPortion/100; // p% of kmers in readpair
    //write to output if p% of kmers below cutoff coverage
    if(counter < limit){
      appendValue(Reads1Buff,*itrReads1);
      appendValue(Reads2Buff,*itrReads2);
      appendValue(Ids1Buff,*itrIds1);
      appendValue(Ids2Buff,*itrIds2);
      appendValue(Quals1Buff,*itrQuals1);
      appendValue(Quals2Buff,*itrQuals2);
      // writeRecord(outfile1,*itrIds1,*itrReads1,*itrQuals1);
      // writeRecord(outfile2,*itrIds2,*itrReads2,*itrQuals2);
    }

    itrReads2++;
    itrIds1++;
    itrIds2++;
    itrQuals1++;
    itrQuals2++;
  }

  reads1=Reads1Buff;
  reads2=Reads2Buff;
  ids1=Ids1Buff;
  ids2=Ids2Buff;
  quals1=Quals1Buff;
  quals2=Quals2Buff;

  // close(outfile1);
  // close(outfile2);

}

void fillIndex32FromReads(StringSet<Dna5String> & reads, coverageIndex32_t & coverageIndex, uint32_t k){
  typedef Iterator<StringSet<Dna5String> >::Type TStringSetIterator;
  Iterator<Dna5String>::Type itread;
  uint64_t hashId;

  for (TStringSetIterator itreads = begin(reads); itreads!=end(reads); ++itreads){ //  for every read in reads

    // add first kmer to index
    Dna5String kmer = prefix(*itreads,k);
    hash32_t hash(kmer);
    hashId = hash.bigger();
    coverageIndex.insert(hashId);

    // iterate over read to extract kmers
    for(itread=begin(*itreads)+k; itread<end(*itreads);itread++){
      // roll hash value
      hash.roll(*itread);
      // add kmer to Index
      hashId = hash.bigger();
      coverageIndex.insert(hashId);
    }
  }

  return;
}

void fillCoverageIndex32(StringSet<Dna5String> & reads1, StringSet<Dna5String> & reads2, coverageIndex32_t & coverageIndex, uint32_t k){
  fillIndex32FromReads(reads1, coverageIndex, k);
  fillIndex32FromReads(reads2, coverageIndex, k);
  return;
}
