# include <seqan/bam_io.h>
# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include "parser.h"
# include "index.h"
using namespace seqan;

int loadReference(StringSet<Dna5String> & seqs, indexOptions & options){
  StringSet<CharString> ids;
  StringSet<IupacString> seqsIn;
  try {
    SeqFileIn file(toCString(options.reference_file));
    readRecords(ids, seqsIn, file);
    close(file);
  }
  catch (ParseError const & e){
    std::cerr << "ERROR: input record is badly formatted. " << e.what() << std::endl;
    return 1;
  }
  catch (IOError const & e){
    std::cerr << "ERROR: input file can not be opened. " << e.what() << std::endl;
    return 1;
  }
  seqs=seqsIn;
  return 0;
}

int loadRefContig(Dna5String & seq, SeqFileIn & file, indexOptions & options){
  CharString id;
  IupacString seqIn;
  try {
    readRecord(id, seqIn, file);
  }
  catch (ParseError const & e){
    std::cerr << "ERROR: input record is badly formatted. " << e.what() << std::endl;
    return 1;
  }
  catch (IOError const & e){
    std::cerr << "ERROR: input file can not be opened. " << e.what() << std::endl;
    return 1;
  }
  seq=seqIn;
  return 0;
}

void readKmerIndex(openAddressingKmerHashtable & Index, std::string & kmer_index_name){

  // defining index file names
  std::string IndPos=kmer_index_name;
  IndPos.append("/pos.txt");
  std::string IndRef=kmer_index_name;
  IndRef.append("/ref.txt");
  std::string IndDir=kmer_index_name;
  IndDir.append("/dir.txt");
  std::string IndC=kmer_index_name;
  IndC.append("/C.txt");
  std::string IndFai=kmer_index_name;
  IndFai.append("/fai.txt");

  //reading Index files in parallel

  #pragma omp parallel for
  for(int i=0;i<4;i++){
    if (i==0){
      String<uint32_t, External<ExternalConfigLarge<>> > extpos;
      if (!open(extpos, IndPos.c_str(), OPEN_RDONLY)){
        throw std::runtime_error("Could not open index position file." );
      }
      assign(Index.pos, extpos, Exact());
      close(extpos);
      std::cerr <<".";
    }
    if (i==1){
      String<uint8_t, External<ExternalConfigLarge<>> > extref;
      if (!open(extref, IndRef.c_str(), OPEN_RDONLY)){
        throw std::runtime_error("Could not open index position file." );
      }
      assign(Index.ref, extref, Exact());
      close(extref);

      std::ifstream input;
      input.open(toCString(IndFai), std::ios::in);
      std::string line;
      while(getline(input,line,'\t')){
        Index.lookChrom.push_back(line);
        getline(input,line);
        Index.chromLength.push_back(line);
      }
      input.close();

      std::cerr <<".";
    }
    if (i==2){
      String<uint32_t, External<> > extdir;
      if (!open(extdir, IndDir.c_str(), OPEN_RDONLY)){
        throw std::runtime_error("Could not open index directory file." );
      }
      assign(Index.dir, extdir, Exact());
      close(extdir);
      std::cerr << ".";
    }
    if (i==3){
      String<int32_t, External<> > extC;
      if (!open(extC, IndC.c_str(), OPEN_RDONLY)){
        throw std::runtime_error("Could not open index counts file." );
      }
      assign(Index.C, extC, Exact());
      close(extC);
      std::cerr << ".";
    }
  } //for omp
  Index.setBucketNumber();
  loadIndexParameters(Index.k,Index.m,kmer_index_name);
  Index.set_k_2();
  return;
}

int checkReadfile(std::string & readfileName){
  try {         // opening read-files
    SeqFileIn file1(toCString(readfileName));
    close(file1);
  }
  catch (ParseError const & e){
    std::cerr << "ERROR: input record is badly formatted. " << e.what() << std::endl;
    return 1;
  }
  catch (IOError const & e){
    std::cerr << "ERROR: input file can not be opened. " << e.what() << std::endl;
    return 1;
  }
  return 0;
}

void writeTmpHeader(longmapOptions & options, openAddressingKmerHashtable & Index, int argc, char const ** argv){
  std::fstream output;
  output.open(options.output_file+".tmp.sam",std::ios::out);

  // HD LINE
  output <<  "@HD\tVN:1.6\tSO:unsorted\n";

  // SQ Lines
  for(int i=0; i!=Index.lookChrom.size();i++){
    output << "@SQ\tSN:" << Index.lookChrom[i] << "\tLN:" << Index.chromLength[i] << "\n";
  }

  // RG Line
  if(options.readGroup!=""){
    output << options.readGroup << "\n";
  }

  // PG Line
  output << "@PG\tID:moleMap\tCL:moleMap";
  for(int i = 0; i < argc; ++i)
        output << ' ' << argv[i];
  output << "\n";
  output.close();

  return;
}

void writeSamHeader(longmapOptions & options){
  BamFileIn bamFileIn(toCString(options.output_file+".tmp.sam"));
  BamHeader header;
  readHeader(header, bamFileIn);
  if(options.output_file=="stdout"){
    BamFileOut bamFileOut(context(bamFileIn), std::cout, Sam());
    writeHeader(bamFileOut, header);
    return;
  }
  std::fstream output;
  output.open(options.output_file, std::ios::out);
  BamFileOut bamFileOut(context(bamFileIn), output, Sam());
  writeHeader(bamFileOut, header);
  return;
}

void writeBamHeader(longmapOptions & options){
  BamFileIn bamFileIn(toCString(options.output_file+".tmp.sam"));
  BamHeader header;
  readHeader(header, bamFileIn);
  if(options.output_file=="stdout"){
    BamFileOut bamFileOut(context(bamFileIn), std::cout, Bam());
    writeHeader(bamFileOut, header);
    return;
  }
  std::fstream output;
  output.open(options.output_file, std::ios::out);
  BamFileOut bamFileOut(context(bamFileIn), output, Bam());
  writeHeader(bamFileOut, header);
  return;
}

void writeSamOutput(std::vector<BamAlignmentRecord> & records, longmapOptions & options){
  BamFileIn bamFileIn(toCString(options.output_file+".tmp.sam"));
  std::fstream output;
  output.open(options.output_file, std::ios::out | std::ios::app);
  BamHeader header;
  readHeader(header, bamFileIn);
  BamFileOut bamFileOut(context(bamFileIn), output, Sam());

  for (std::vector<BamAlignmentRecord>::iterator itr_res=records.begin(); itr_res < records.end(); itr_res++){
    writeRecord(bamFileOut, (*itr_res));
  }
  records.clear();
  return;
}

void writeBamOutput(std::vector<BamAlignmentRecord> & records, longmapOptions & options){
  BamFileIn bamFileIn(toCString(options.output_file+".tmp.sam"));
  std::fstream output;
  output.open(options.output_file, std::ios::out | std::ios::app);
  BamHeader header;
  readHeader(header, bamFileIn);
  BamFileOut bamFileOut(context(bamFileIn), output, Bam());

  for (std::vector<BamAlignmentRecord>::iterator itr_res=records.begin(); itr_res < records.end(); itr_res++){
    writeRecord(bamFileOut, (*itr_res));
  }
  records.clear();
  return;
}

void writeSamCout(std::vector<BamAlignmentRecord> & records, longmapOptions & options){
  BamFileIn bamFileIn(toCString(options.output_file+".tmp.sam"));
  BamHeader header;
  readHeader(header, bamFileIn);
  BamFileOut bamFileOut(context(bamFileIn), std::cout, Sam());

  for (std::vector<BamAlignmentRecord>::iterator itr_res=records.begin(); itr_res < records.end(); itr_res++){
    writeRecord(bamFileOut, (*itr_res));
  }
  records.clear();
  return;
}

void writeBamCout(std::vector<BamAlignmentRecord> & records, longmapOptions & options){
  BamFileIn bamFileIn(toCString(options.output_file+".tmp.sam"));
  BamHeader header;
  readHeader(header, bamFileIn);
  BamFileOut bamFileOut(context(bamFileIn), std::cout, Bam());

  for (std::vector<BamAlignmentRecord>::iterator itr_res=records.begin(); itr_res < records.end(); itr_res++){
    writeRecord(bamFileOut, (*itr_res));
  }
  records.clear();
  return;
}

// void writeBedOutput(std::vector<BamAlignmentRecord> & records, std::vector<std::string> & lookChrom,  longmapOptions & options){
//   std::fstream output;
//   output.open(options.output_file,std::ios::out | std::ios::app);
//   int qual;
//   for (std::vector<BamAlignmentRecord>::iterator itr_res=records.begin(); itr_res < records.end(); itr_res++){
//     extractTagValue(qual, (*itr_res).tags, "MS");
//     output << lookChrom[(int)(*itr_res).rID] << "\t"
//            << (*itr_res).beginPos << "\t"
//            << (*itr_res).beginPos+(int)(*itr_res).tLen << "\t"
//            << (*itr_res).qName << "\t"
//            << qual << "\n";
//   }
//   output.close();
//   return;
// }

void writeHistogram(std::vector<uint32_t> & histogram, std::string & output_file){
  std::ofstream file_histogram;
  file_histogram.open(output_file+".hist");
  for (int i=0; i<histogram.size(); i++){
    file_histogram << std::to_string(histogram[i]) << "\n";
  }
  file_histogram.close();
  return;
}
