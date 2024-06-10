# include <seqan/seq_io.h>
# include <seqan/sequence.h>
# include <seqan/bam_io.h>
# include <iostream>
# include <fstream>
# include "functions.h"
# include "mapping_functions.h"
# include "parser.h"
# include "index.h"
# include "IO_functions.h"
# include "functions.h"
# include "coverage_analysis.h"
# include <zlib.h>
# include <htslib/kseq.h>


int getBarcodeLengthGz(kseq_t * seq1, kseq_t * seq2){
  std::string barcode="*";
  std::string comment;
  while(barcode=="*"){
    if(kseq_read(seq1)>=0){
      kseq_read(seq2);
      std::string id1=seq1->comment.s;
      barcode=id1.substr(id1.find("BX:Z:")+5,id1.find(' ',id1.find("BX:Z:")+5)-id1.find("BX:Z:")-5);
    }else{break;}
  }

  return barcode.size();
}

std::string getBarcodeGz(std::string id, int barcode_length){
  return id.substr(id.find("BX:Z:")+5,barcode_length);
}

void readPairedBatch(kseq_t * seq1, kseq_t * seq2, std::vector<barcodeData> & Batch, int & batchSize, int & barcode_length){
  // extracts batches of Barcodes
  int barcodeCount = 0;
  barcodeData barcodeBuff; // barcode and set of reads of this barcode
  barcodeBuff.barcode=getBarcodeGz(seq1->comment.s, barcode_length);
  barcodeBuff.append(seq1, seq2);
  std::string newBarcode;
  Batch.clear();

  while(Batch.size()<batchSize){
    if(kseq_read(seq1)<0){
      break;
    }
    kseq_read(seq2);
    newBarcode = getBarcodeGz(seq1->comment.s, barcode_length);
    if(newBarcode == barcodeBuff.barcode){
      barcodeBuff.append(seq1, seq2);
    }else{ // save barcode data and initialize new  buff
      Batch.push_back(barcodeBuff);
      barcodeBuff.Reads.clear();
      barcodeBuff.barcode=newBarcode;
      barcodeBuff.append(seq1, seq2);
    }
  }
  return;
}

int mapLinkedZipped(openAddressingKmerHashtable & Index, mapOptions & options, int argc, char const ** argv){

  std::cerr << "Processing read file...";

  omp_lock_t lock;
  omp_lock_t reslock;
  omp_init_lock(&lock);
  omp_init_lock(&reslock);

  std::vector<uint32_t> histogram(200,0);
  std::vector<result_t> globalresults;
  uint8_t window_count=50;

  omp_set_num_threads(options.threads);
  if(options.threads < 2){options.threads=2;} // make sure that the program works without parallelization

  gzFile readFile1 = gzopen(toCString(options.readfile1), "r");
  gzFile readFile2 = gzopen(toCString(options.readfile2), "r");
  kseq_t * seq1 = kseq_init(readFile1);
  kseq_t * seq2 = kseq_init(readFile2);

  int barcode_length = getBarcodeLengthGz(seq1, seq2); // determine barcode length and skip past '*' barcodes
  int batchSize = (options.threads-1)*1000;

  // read first batch
  std::vector<barcodeData> newBatch;
  std::vector<barcodeData> oldBatch;
  // read first batch

  readPairedBatch(seq1, seq2, oldBatch, batchSize, barcode_length);
  while(oldBatch.size()){
    #pragma omp parallel for
    for (int t=0; t<options.threads; t++){
      if(t==0){ // read next Batch
        readPairedBatch(seq1, seq2, newBatch, batchSize, barcode_length);
      }else{ // analyse Batch
        // set scope for thread
        minimizedSequence miniSeq(options.k,options.mini_window_size);
        int start=(t-1)*(oldBatch.size()/(options.threads-1));
        int end=(t)*(oldBatch.size()/(options.threads-1));
        if(t==options.threads-1){end=oldBatch.size();}

        std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t>> kmer_list;   // (i,j,a,m_a)   i=reference (Chromosome), j=position of matching k-mer in reference, a=abundance of k-mer in reference, m_a=minimizer_active_bases
        std::vector<result_t> results;
        std::vector<uint32_t> histogram_local(200,0);

        for(int i = start; i!=end; i++){ // for every barcode in batch
          std::vector<Dna5String>::iterator itReads = oldBatch[i].Reads.begin();
          for(itReads; itReads!=oldBatch[i].Reads.end(); itReads++){ // for ever read of the barcode
            // initialize minimizer class
            minimizer mini;
            miniSeq.init(*itReads);
            // return and save all minimizers to kmer_list
            while(!miniSeq.at_end){
              mini=miniSeq.pop();
              AppendPos(kmer_list, mini.value^miniSeq.random_seed, Index, mini.active_bases, options);
            }
          }
          // map barcode and clear k_mer list
          if (!kmer_list.empty()) {
            std::sort(kmer_list.begin(),kmer_list.end());
            MapKmerListLinked(kmer_list, window_count, oldBatch[i].barcode, results, Index.lookChrom, histogram_local, options);
            kmer_list.clear();
          } // if (!kmer_list.empty())

        } // for every barcode in batch

        if (options.Sort==0){
          omp_set_lock(&lock);
          std::fstream output;
          output.open(options.output_file,std::ios::out | std::ios::app);
          for (std::vector<result_t>::iterator itr_res=results.begin(); itr_res < results.end(); itr_res++){
            output << (*itr_res).string();
          }
          results.clear();
          output.close();
          omp_unset_lock(&lock);
        }else{ // sort output and then merge into global output
          sortResults(results);
          omp_set_lock(&lock);
          uint64_t middle = globalresults.size();
          globalresults.insert(globalresults.end(),results.begin(),results.end());
          std::inplace_merge(globalresults.begin(), globalresults.begin()+middle, globalresults.end(), compFunctionResult);
          omp_unset_lock(&lock);
        }
        results.clear();

        omp_set_lock(&reslock);
        for (int i=0; i<histogram.size(); i++){
          histogram[i]+=histogram_local[i];
        }
        omp_unset_lock(&reslock);

      }
    } // parallel for loop
    oldBatch = newBatch;
    newBatch.clear();
  } // while(oldBatch.size())

  if(options.CoverageAnalysis==1){ // perform coverage analysis
    coverageAnalysis(globalresults, histogram, options);
  }

  if(options.Sort==1){ // write global results to file
    std::fstream output;
    output.open(options.output_file,std::ios::out | std::ios::app);
    for (std::vector<result_t>::iterator itr_res=globalresults.begin(); itr_res < globalresults.end(); itr_res++){
      output << (*itr_res).string();
    }
    output.close();
  }

  // write histogram to file
  writeHistogram(histogram, options.output_file);

  std::cerr << ".........done.\n";

  std::cerr << "Barcodes mapped sucessfully!\n\n";

  return 0;
}

void MapKmerListLinked(std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, uint8_t & window_count, std::string barcode, std::vector<result_t> & results, std::vector<std::string> & lookChrom, std::vector<uint32_t> & histogram, mapOptions & options){
  unsigned qualityThreshold=50;
  std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrk;

  // std::vector<float> lookQual={0,20,3.0, 0.75, 0.38, 0.24, 0.17, 0.14, 0.11, 0.09, 0.08};
  std::vector<float> lookQual={0,5,1.44, 0.91, 0.72, 0.62, 0.56, 0.51, 0.48, 0.46, 0.43, 0.42, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.35, 0.34, 0.33}; // log(1/abundance)
  lookQual.resize(options.max_abundance, 0.33);

  #define REF(X) std::get<0>(*(X))
  #define POS(X) std::get<1>(*(X))
  #define ABU(X) std::get<2>(*(X))
  #define ACT(X) std::get<3>(*(X))

  std::vector<std::tuple<double,uint8_t,uint32_t,uint32_t>> best_windows(window_count,std::make_tuple(0,0,0,0)); //(maping_quality, reference, start position in referende, end position)
  std::vector<std::tuple<double,uint8_t,uint32_t,uint32_t>>::iterator itrbw;
  uint8_t reference=REF(kmer_list.begin());
  uint32_t start_position=POS(kmer_list.begin());
  uint32_t end_position=POS(kmer_list.begin());
  std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrstart=kmer_list.begin();
  double window_quality=0;
  std::tuple<double,uint8_t,uint32_t,uint32_t> candidate=std::make_tuple(0,0,0,4294967295); //(maping_quality, reference, start position in referende, end position)

  window_quality+=lookQual[ABU(kmer_list.begin())]*ACT(kmer_list.begin()); // lookQual = (1/(log(abund)^3))*minimizer_active_bases
  for(itrk=kmer_list.begin()+1;itrk!=kmer_list.end();itrk++){ //iterating over kmer listed

    if (/*end position*/std::get<3>(candidate) < start_position) { // if current window no longer overlaps the qualifiing window
      if (std::get<0>(candidate)>qualityThreshold){
        trimmWindowLinked(kmer_list, itrstart, itrk, candidate, lookQual);
        ReportWindow(best_windows,candidate);
      }
    }

    if (reference==REF(itrk) && (POS(itrk)-start_position) < options.max_window_size && (POS(itrk)-end_position) < options.max_gap_size) { //checking break criteria
      //append window by kmer_hit
      window_quality+=lookQual[ABU(itrk)]*ACT(itrk);
      end_position=POS(itrk);

    }else if (REF(itrk)!=reference || (POS(itrk)-end_position) >= options.max_gap_size){  // if k_mer hit from next reference or gapsize to large: report current window or candiadate window and initialize new window

      if(window_quality > std::get<0>(candidate)){ // report current window or candidate
        candidate=std::make_tuple(window_quality,reference,start_position,end_position);
        if (std::get<0>(candidate)>qualityThreshold){
          trimmWindowLinked(kmer_list, itrstart, itrk, candidate, lookQual);
          ReportWindow(best_windows,candidate);
        }
      }else{
        if (std::get<0>(candidate)>qualityThreshold){
          trimmWindowLinked(kmer_list, itrstart, itrk, candidate, lookQual);
          ReportWindow(best_windows,candidate);
        }
      }
      window_quality=lookQual[ABU(itrk)]*ACT(itrk);
      itrstart=itrk;
      reference=REF(itrk);
      start_position=POS(itrk);
      end_position=POS(itrk);

    }else if ((POS(itrk)-start_position) >= options.max_window_size){ // maximum window size criterion hurt: shrink window from start and save better one as candidate
      if (window_quality > std::get<0>(candidate)) { // check if current window better than candidate: if yes: replace candidate
        candidate=std::make_tuple(window_quality,reference,start_position,end_position);
      }
      window_quality+=lookQual[ABU(itrk)]*ACT(itrk);
      end_position=POS(itrk);
      while (POS(itrk)-POS(itrstart)>options.max_window_size){ //shrinking window untill max_window_size criterion met
        window_quality-=lookQual[ABU(itrstart)]*ACT(itrstart);
        itrstart++;
      }
      start_position=POS(itrstart);

    }else{
      std::cerr << "\nREF(itrk): " << REF(itrk) << " reference: " << reference << "\nPOS(itrk): " << POS(itrk) << " start_position: " << start_position << " end_position: " << end_position << "\n";
      std::cerr << "POS(itrstart): " << POS(itrstart) << " REF(itrstart): " << REF(itrstart) << "\n";
    }
  }
  if(window_quality > std::get<0>(candidate)){ // report last window or last candidate
    candidate=std::make_tuple(window_quality,reference,start_position,end_position);
    if (std::get<0>(candidate)>qualityThreshold){
      trimmWindowLinked(kmer_list, itrstart, itrk, candidate, lookQual);
      ReportWindow(best_windows,candidate);
    }
  }else{
    if (std::get<0>(candidate)>qualityThreshold){
      trimmWindowLinked(kmer_list, itrstart, itrk, candidate, lookQual);
      ReportWindow(best_windows,candidate);
    }
  }
  //filter low quality windows
  if (std::get<0>(*(best_windows.end()-1))!=0) {
    while(std::get<0>(*best_windows.begin())<qualityThreshold && !best_windows.empty()){
      best_windows.erase(best_windows.begin());
    }
  }else{return;}

  // filter short windows
  std::vector<int> toshort;
  for (int i = 0; i!=best_windows.size(); i++){
    if ((std::get<3>(best_windows[i])-std::get<2>(best_windows[i]))<options.l){ // options.l = lengths threshold
      toshort.push_back(i);
    }
  }
  for (int i=(toshort.size()-1);i>=0;--i) {
    best_windows.erase(best_windows.begin()+toshort[i]);
  }

  for(itrbw=best_windows.begin();itrbw!=best_windows.end(); itrbw++){
    float qual=roundf((float)(std::get<0>(*itrbw)/(std::get<3>(*itrbw)-std::get<2>(*itrbw)))*100); // quality/length_of_mapping--> QualityPerBase
    if(qual>200){
      histogram[200]++;
    }else{
      histogram[qual]++;
    }
    if(qual>=options.s){ // options.s = scoreThreshold
      results.push_back({lookChrom[std::get<1>(*itrbw)], std::get<2>(*itrbw), std::get<3>(*itrbw), barcode ,(uint16_t)qual});
    }
  }
  return;
}

void trimmWindowLinked(std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrstart, std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrk, std::tuple<double,uint8_t,uint32_t,uint32_t> & candidate, std::vector<float> & lookQual){
  double edge_quality=0;
  double min_edge_qual=100; //200
  uint32_t edge_len=2000; //2000
  uint32_t counter=0;
  double window_quality=std::get<0>(candidate);
  itrk--;

  //trimm start
  while(edge_quality<min_edge_qual){
    if(itrstart+counter>=itrk){
      candidate=std::make_tuple(0,0,0,4294967295);
      return;
    }
    edge_quality+=lookQual[ABU(itrstart+counter)]*ACT(itrstart+counter);
    while(POS(itrstart+counter)-POS(itrstart)>edge_len){ //shrinking window until edge_len criterion met
      edge_quality-=lookQual[ABU(itrstart)]*ACT(itrstart);
      window_quality-=lookQual[ABU(itrstart)]*ACT(itrstart);
      itrstart++;
      counter--;
    }
    counter++;
  }
  edge_quality=0;
  counter=0;

  //trimm end
  while(edge_quality<min_edge_qual){
    if(itrk-counter<=itrstart){ //no area in candidate has the required minimum score/base
      candidate=std::make_tuple(0,0,0,4294967295);
      return;
    }
    edge_quality+=lookQual[ABU(itrk-counter)]*ACT(itrk-counter);
    while(POS(itrk)-POS(itrk-counter)>edge_len){ //shrinking window until edge_len criterion met
      edge_quality-=lookQual[ABU(itrk)]*ACT(itrk);
      window_quality-=lookQual[ABU(itrk)]*ACT(itrk);
      itrk--;
      counter--;
    }
    counter++;
  }
  candidate=std::make_tuple(window_quality,REF(itrstart),POS(itrstart),POS(itrk));
  return;
}

void AppendPos(std::vector<std::tuple <uint8_t,uint32_t,uint32_t,uint32_t>> & kmer_list, const int64_t & hash, openAddressingKmerHashtable & Index, uint8_t & minimizer_active_bases, mapOptions & options){
      uint32_t c=Index.GetBkt(hash);
      uint32_t abundance=Index.dir[c+1]-Index.dir[c];
      if (abundance<=options.max_abundance){
        kmer_list.reserve(kmer_list.size()+abundance);
        for (uint32_t i = Index.dir[c];i!=Index.dir[c+1];i++){
          kmer_list.push_back(std::make_tuple(Index.ref[i],Index.pos[i],abundance,minimizer_active_bases));
        }
      }
      return;
}

void skipToFirstBarcode(kseq_t * seq1, kseq_t * seq2, uint8_t & barcode_length){ // skip through file until barcode[0]!="*"
  std::string barcode;
  while(kseq_read(seq1)>=0){
    kseq_read(seq2);
    barcode=getBarcode(seq1->name.s, barcode_length);
    if(barcode[0]!='*'){
      break;
    }
  }
  return;
}

uint8_t getBarcodeLength(std::string & readfile1, std::streampos & readfile1_size){
  SeqFileIn file1(toCString(readfile1));
  Dna5String read1;
  CharString id;
  readRecord(id,read1,file1);
  std::string id1;
  std::string line="X";
  uint8_t barcode_length;
  std::string barcode;
  for (uint64_t i=round((long long unsigned)readfile1_size/20); i<(uint64_t)readfile1_size;i+=round((uint64_t)readfile1_size/10)){ // check multiple positions throughout the readfile
    file1.stream.file.seekg(i);
    while(line!="+"){
      getline(file1.stream.file,line);
    }
    file1.stream.file.ignore(10000,'\n');
    readRecord(id, read1, file1);
    id1=toCString(id);
    barcode=id1.substr(id1.find("BX:Z:")+5,id1.find(' ',id1.find("BX:Z:")+5)-id1.find("BX:Z:")-5);              // determine BX:Z: entry
    if (barcode[0]=='A' || barcode[0]=='C' || barcode[0]=='T' || barcode[0]=='G'){
      close(file1);
      // std::cerr << "barcode length   \t" << barcode.size() << "\n\n";
      return barcode.size();
    }
  }
  std::cerr << "\nERROR: Incorrect BX:Z: tag in readfile. Barcodes can not be determined. Please check your input file formating.\n";
  return 0;
}

// ##################################################################################

int mapLongUnzipped(openAddressingKmerHashtable & Index, longmapOptions & options, int argc, char const ** argv){
  std::cerr << "Processing read file...";

  writeTmpHeader(options, Index, argc, argv);
  void (*writeOutput)(std::vector<BamAlignmentRecord> & records, longmapOptions & options);
  if(options.output_format=="sam"){
    writeSamHeader(options);
    if(options.output_file=="stdout"){
      writeOutput = &writeSamCout;
    }else{
      writeOutput = &writeSamOutput;
    }
  }
  if(options.output_format=="bam"){
    writeBamHeader(options);
    if(options.output_file=="stdout"){
      writeOutput = &writeBamCout;
    }else{
      writeOutput = &writeBamOutput;
    }
  }

  //checking file size
  SeqFileIn file1(toCString(options.readfile1));
  CharString (*getIdFunction)(std::string); // function pointer
  getIdFunction = &getID;
  file1.stream.file.seekg(0, std::ios::end);
  std::streampos readfile1_size=file1.stream.file.tellg();
  close(file1);

  std::vector<uint32_t> histogram(200,0);

  omp_lock_t lock; omp_lock_t reslock; omp_init_lock(&lock); omp_init_lock(&reslock);
  omp_set_num_threads(options.threads);

  #pragma omp parallel for ordered
  for (int t=0; t<options.threads; t++){
    //declare variables
    std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>> kmer_list;   // (i,j,a,m_a,o)   i=reference (Chromosome), j=position of matching k-mer in reference, a=abundance of k-mer in reference, m_a=minimizer_active_bases, o=order_of_kmers_in_seq
    Dna5String read; CharString id; CharString qual;
    BamAlignmentRecord result;
    std::vector<BamAlignmentRecord> results = {};
    std::vector<uint32_t> histogram_local(200,0);
    minimizedSequence miniSeq(options.k,options.mini_window_size);

    //open readfile and move to start position
    SeqFileIn file1(toCString(options.readfile1));
    std::streampos startpos=readfile1_size/options.threads*t;
    std::streampos endpos=readfile1_size/options.threads*(t+1);
    moveFileToStart(file1, startpos, t);
    // if skiped beyond the boundaries of threads scope: end thread
    if (file1.stream.file.tellg()>endpos){continue;}

    //proceed through readfile until endpos
    while (!atEnd(file1)) { // proceeding through files
      readRecord(id, read, qual, file1);
      // initialize minimizer class
      minimizer mini;
      miniSeq.init(read);
      // return and save all minimizers to kmer_list
      uint32_t order=0;
      while(!miniSeq.at_end){
        mini=miniSeq.pop();
        AppendPosLong(kmer_list, mini.value^miniSeq.random_seed, Index, mini.active_bases, order, options);
        order++;
      }

      // map read and clear k_mer list
      result.qName=id;
      result.seq=read;
      result.qual=qual;
      if(options.readGroupId!=""){
        appendTagValue(result.tags, "RG", options.readGroupId);
      }
      if (!kmer_list.empty()){
        std::sort(kmer_list.begin(),kmer_list.end());

        if(MapKmerListLong(kmer_list, result, histogram_local, length(read), options)){
          results.push_back(result);
        }else{//if MapKmerListLong returns 0 read is unmapped
          result.flag=4;
          results.push_back(result);
        }
        clear(result);

        kmer_list.clear();
        if (results.size()>10) {
          if (omp_test_lock(&lock)){
            writeOutput(results, options);
            omp_unset_lock(&lock);
          }else if(results.size()>100){
            omp_set_lock(&lock);
            writeOutput(results, options);
            omp_unset_lock(&lock);
          }
        }
      }else{
        result.flag=4;
        results.push_back(result);
      }

      if (file1.stream.file.tellg()>endpos){break;} // break if range for this thread exceeded
    } // while(!atEnd(file1))
    close(file1);

    // write final output and sum up local histograms
    omp_set_lock(&lock);
    writeOutput(results, options);
    for (int i=0; i<histogram.size(); i++){
      histogram[i]+=histogram_local[i];
    }
    omp_unset_lock(&lock);

  } // parallel loop

  // write histogram to file
  remove(toCString(options.output_file+".tmp.sam"));
  writeHistogram(histogram, options.output_file);

  std::cerr << ".........done.\n";
  return 0;
}

void readBatch(kseq_t * seq1, std::vector<ReadData> & Batch, int & batchSize){
  Batch.resize(batchSize);
  for(int i=0; i!=batchSize; i++){
    if(kseq_read(seq1)>=0){
      Batch[i].id=seq1->name.s;
      Batch[i].read=seq1->seq.s;
      Batch[i].qual=seq1->qual.s;
    }else{
      Batch.resize(i);
      break;
    }
  }
  return;
}

int mapLongZipped(openAddressingKmerHashtable & Index, longmapOptions & options, int argc, char const ** argv){
  std::cerr << "Processing read file...";

  writeTmpHeader(options, Index, argc, argv);
  void (*writeOutput)(std::vector<BamAlignmentRecord> & records, longmapOptions & options);
  if(options.output_format=="sam"){
    writeSamHeader(options);
    if(options.output_file=="stdout"){
      writeOutput = &writeSamCout;
    }else{
      writeOutput = &writeSamOutput;
    }
  }
  if(options.output_format=="bam"){
    writeBamHeader(options);
    if(options.output_file=="stdout"){
      writeOutput = &writeBamCout;
    }else{
      writeOutput = &writeBamOutput;
    }
  }

  std::vector<uint32_t> histogram(200,0);
  omp_lock_t lock; omp_lock_t reslock; omp_init_lock(&lock); omp_init_lock(&reslock);

  gzFile readFile = gzopen(toCString(options.readfile1), "r");
  kseq_t * seq1 = kseq_init(readFile);

  // read first batch
  int batchSize=(options.threads-1)*5000;
  std::vector<ReadData> oldBatch;
  std::vector<ReadData> newBatch;
  readBatch(seq1, oldBatch, batchSize);

  omp_set_num_threads(options.threads);
  if(options.threads < 2){options.threads=2;} // make sure that the program works without parallelization

  while(oldBatch.size()){
    #pragma omp parallel for ordered
    for (int t=0; t<options.threads; t++){
      if(t==0){
        // read new batch
        readBatch(seq1, newBatch, batchSize);
      }else{ // t!=0  // process threads part of last batch
        // set scope for thread
        minimizedSequence miniSeq(options.k,options.mini_window_size);
        int start=(t-1)*(oldBatch.size()/(options.threads-1));
        int end=(t)*(oldBatch.size()/(options.threads-1));
        if(t==options.threads-1){end=oldBatch.size();}

        //declare variables
        std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>> kmer_list;   // (i,j,a,m_a,o)   i=reference (Chromosome), j=position of matching k-mer in reference, a=abundance of k-mer in reference, m_a=minimizer_active_bases, o=order_of_kmers_in_seq
        BamAlignmentRecord result;
        std::vector<BamAlignmentRecord> results = {};
        std::vector<uint32_t> histogram_local(200,0);

        // iterate over threads scope
        for(int i = start; i!=end; i++){
          // Map read
          minimizer mini;
          miniSeq.init(oldBatch[i].read);
          //save all minimizers to kmer_list
          uint32_t order=0;
          while(!miniSeq.at_end){
            mini=miniSeq.pop();
            AppendPosLong(kmer_list, mini.value^miniSeq.random_seed, Index, mini.active_bases, order, options);
            order++;
          }

          // map read and clear k_mer list
          result.qName=oldBatch[i].id;
          result.seq=oldBatch[i].read;
          result.qual=oldBatch[i].qual;
          if(options.readGroupId!=""){
            appendTagValue(result.tags, "RG", options.readGroupId);
          }

          if (!kmer_list.empty()){
            std::sort(kmer_list.begin(),kmer_list.end());

            if(MapKmerListLong(kmer_list, result, histogram_local, length(oldBatch[i].read), options)){
              results.push_back(result);
            }else{//if MapKmerList returns 0 read is unmapped
              result.flag=4;
              results.push_back(result);
            }
            clear(result);

            kmer_list.clear();
            if (results.size()>10) {
              if (omp_test_lock(&lock)){
                writeOutput(results, options);
                omp_unset_lock(&lock);
              }else if(results.size()>100){
                omp_set_lock(&lock);
                writeOutput(results, options);
                omp_unset_lock(&lock);
              }
            }
          }else{
            result.flag=4;
            results.push_back(result);
          }
        } // for(int i = start; i!=end; i++) For read in scope
        // write final output and sum up local histograms
        omp_set_lock(&lock);
        writeOutput(results, options);
        for (int i=0; i<histogram.size(); i++){
          histogram[i]+=histogram_local[i];
        }
        omp_unset_lock(&lock);
      } // if not thread 0
    } // for thread in threads
    oldBatch=newBatch;
  } //while not at end of file

  // write histogram to file
  remove(toCString(options.output_file+".tmp.sam"));
  writeHistogram(histogram, options.output_file);
  std::cerr << ".........done.\n";

  return 0;
}

void moveFileToStart(SeqFileIn & file1, std::streampos & startpos, int & t){
  //move file 1 to start position
  Dna5String read;
  CharString id;
  if (t!=0){
    std::string line="-";
    std::streampos temppos = startpos;
    file1.stream.file.seekg(startpos);
    while(true){
      temppos = file1.stream.file.tellg();
      getline(file1.stream.file,line);
      if(line[0]=='@'){
        getline(file1.stream.file,line);
        getline(file1.stream.file,line);
        if(line[0]=='+'){
          file1.stream.file.seekg(temppos);
          break;
        }else{
          file1.stream.file.seekg(temppos);
          getline(file1.stream.file,line);
        }
      }
    } // while(true)
    try{
      readRecord(id, read, file1);
    }
    catch (Exception const & e){
      std::cerr << "ERROR: " << e.what() << std::endl;
    }
    file1.stream.file.seekg(temppos);

  } else { // if t=0
    file1.stream.file.seekg(0);
    try{
      readRecord(id, read, file1);
    }
    catch (Exception const & e){
      std::cerr << "ERROR: " << e.what() << std::endl;
    }
    file1.stream.file.seekg(0);
  }
  return;
}

int MapKmerListLong(std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>> & kmer_list, BamAlignmentRecord & result, std::vector<uint32_t> & histogram, uint32_t readLength, longmapOptions & options){
  unsigned qualityThreshold=0;
  uint32_t max_window_size=readLength*5;
  uint32_t window_count=10;
  std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t, uint32_t>>::const_iterator itrk=kmer_list.begin();
  // std::vector<float> lookQual={0,5,1.44, 0.91, 0.72, 0.62, 0.56, 0.51, 0.48, 0.46, 0.43, 0.42, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.35, 0.34, 0.33}; // log(1/abundance)
  std::vector<float> lookQual={0,4,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}; // log(1/abundance)
  lookQual.resize(options.max_abundance, 1);

  #define REF(X) std::get<0>(*(X)) // ref contig
  #define POS(X) std::get<1>(*(X)) // position of kmer in ref
  #define ABU(X) std::get<2>(*(X)) // abundance of kmer in ref
  #define ACT(X) std::get<3>(*(X)) // active time of kmer in seq
  #define ORD(X) std::get<4>(*(X)) // order within sequence

  // std::vector<std::tuple<double,uint8_t,uint32_t,uint32_t>> bestWindows(window_count, std::make_tuple(0,0,0,0)); //(maping_quality, reference, start position in referende, end position)
  std::tuple<double,uint8_t,uint32_t,uint32_t,bool> best_window = std::make_tuple(0,0,0,0,true); //(maping_quality, reference, start position in referende, end position)
  std::tuple<double,uint8_t,uint32_t,uint32_t,bool> candidate=std::make_tuple(0,0,0,4294967295,true); //(maping_quality, reference, start position in referende, end position)
  uint8_t reference=REF(itrk);
  uint32_t start_position=POS(itrk);
  uint32_t end_position=POS(itrk);
  std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrstart=kmer_list.begin();
  double window_quality=0;

  window_quality+=lookQual[ABU(kmer_list.begin())]*ACT(kmer_list.begin()); // lookQual = (1/(log(abund)^3))*minimizer_active_bases

  for(itrk;itrk!=kmer_list.end();itrk++){ //iterating over kmer listed

    if (/*end position*/std::get<3>(candidate) < start_position) { // if current window no longer overlaps the qualifiing window
      if (std::get<0>(candidate)>qualityThreshold){
        trimmWindowLong(kmer_list, itrstart, itrk, candidate, lookQual);
        if(std::get<0>(candidate) > std::get<0>(best_window)){best_window=candidate;}
      }
    }

    if (reference==REF(itrk) && (POS(itrk)-start_position) < max_window_size && (POS(itrk)-end_position) < options.max_gap_size) { //checking break criteria
      //append window by kmer_hit
      window_quality+=lookQual[ABU(itrk)]*ACT(itrk);
      end_position=POS(itrk);
      // order+=1-2*(old_ord<ORD(itrk));

    }else if (REF(itrk)!=reference || (POS(itrk)-end_position) >= options.max_gap_size){  // if k_mer hit from next reference or gapsize to large: report current window or candiadate window and initialize new window

      if(window_quality > std::get<0>(candidate)){ // report current window or candidate
        candidate=std::make_tuple(window_quality,reference,start_position,end_position,true);
        if (std::get<0>(candidate)>qualityThreshold){
          trimmWindowLong(kmer_list, itrstart, itrk, candidate, lookQual);
          if(std::get<0>(candidate) > std::get<0>(best_window)){best_window=candidate;}
        }
      }else{
        if (std::get<0>(candidate)>qualityThreshold){
          trimmWindowLong(kmer_list, itrstart, itrk, candidate, lookQual);
          if(std::get<0>(candidate) > std::get<0>(best_window)){best_window=candidate;}
        }
      }
      window_quality=lookQual[ABU(itrk)]*ACT(itrk);
      itrstart=itrk;
      reference=REF(itrk);
      start_position=POS(itrk);
      end_position=POS(itrk);
      // order=0;

    }else if ((POS(itrk)-start_position) >= max_window_size){ // maximum window size criterion hurt: shrink window from start and save better one as candidate
      if (window_quality > std::get<0>(candidate)) { // check if current window better than candidate: if yes: replace candidate
        candidate=std::make_tuple(window_quality,reference,start_position,end_position,true);
      }
      window_quality+=lookQual[ABU(itrk)]*ACT(itrk);
      end_position=POS(itrk);
      while (POS(itrk)-POS(itrstart)>max_window_size){ //shrinking window untill max_window_size criterion met
        window_quality-=lookQual[ABU(itrstart)]*ACT(itrstart);
        itrstart++;
      }
      start_position=POS(itrstart);

    }else{
      std::cerr << "\nREF(itrk): " << REF(itrk) << " reference: " << reference << "\nPOS(itrk): " << POS(itrk) << " start_position: " << start_position << " end_position: " << end_position << "\n";
      std::cerr << "POS(itrstart): " << POS(itrstart) << " REF(itrstart): " << REF(itrstart) << "\n";
    }
  } // for kmer hit in kmer list

  if(window_quality > std::get<0>(candidate)){ // report last window or last candidate
    candidate=std::make_tuple(window_quality,reference,start_position,end_position,true);
    if (std::get<0>(candidate)>qualityThreshold){
      trimmWindowLong(kmer_list, itrstart, itrk, candidate, lookQual);
      if(std::get<0>(candidate) > std::get<0>(best_window)){best_window=candidate;}
    }
  }else{
    if (std::get<0>(candidate)>qualityThreshold){
      trimmWindowLong(kmer_list, itrstart, itrk, candidate, lookQual);
      if(std::get<0>(candidate) > std::get<0>(best_window)){best_window=candidate;}
    }
  }

  //filter low quality windows
  if (std::get<0>(best_window)<qualityThreshold || std::get<0>(best_window)==0) {
    return 0;
  }

  // filter short windows
  if ((std::get<3>(best_window)-std::get<2>(best_window))<options.l){
    return 0;
  }

  // report best window
  // float qual=roundf((float)(std::get<0>(best_window)/(std::get<3>(best_window)-std::get<2>(best_window)))*100); // quality/length_of_mapping--> QualityPerBase
  int qual=roundf((float)(std::get<0>(best_window))*10); //
  // float qual=roundf((float)(std::get<0>(best_window)/(readLength))*100); // quality/length_of_mapping--> QualityPerBase
  result.rID=std::get<1>(best_window);
  result.beginPos=std::get<2>(best_window);
  result.tLen=(int32_t)(std::get<3>(best_window)-std::get<2>(best_window));
  result.flag=16*(std::get<4>(best_window));
  result.cigar=seqan::CigarElement<char>('M',readLength);

  appendTagValue(result.tags, "MS", std::abs(qual));
  if(options.readGroupId!=""){
    appendTagValue(result.tags, "RG", options.readGroupId);
  }
  // result.mapQ=(uint8_t)qual;
  // MS TAG (MAPPING SCORE)

  if(result.mapQ < options.s){ // options.s = scoreThreshold
    return 0;
  }

  // add to quality histogram
  if(qual>200){
    histogram[200]++;
  }else{
    histogram[qual]++;
  }

  return 1;
}

void trimmWindowLong(std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>> & kmer_list, std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrstart, std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrk, std::tuple<double,uint8_t,uint32_t,uint32_t,bool> & candidate, std::vector<float> & lookQual){
  double edge_quality=0;
  double min_edge_qual=20; //100
  uint32_t edge_len=4000; //2000
  uint32_t counter=0;
  double window_quality=std::get<0>(candidate);
  itrk--;

  // std::tuple<double,uint8_t,uint32_t,uint32_t> backupCandidate=candidate; // backup of candidate in case the trimming deletes the candidate

  //trimm start
  while(edge_quality<min_edge_qual){
    if(itrstart+counter>=itrk){
      candidate=std::make_tuple(0,0,0,4294967295,true);
      return;
    }
    edge_quality+=lookQual[ABU(itrstart+counter)]*ACT(itrstart+counter);
    while(POS(itrstart+counter)-POS(itrstart)>edge_len){ //shrinking window until edge_len criterion met
      edge_quality-=lookQual[ABU(itrstart)]*ACT(itrstart);
      window_quality-=lookQual[ABU(itrstart)]*ACT(itrstart);
      itrstart++;
      counter--;
    }
    counter++;
  }
  edge_quality=0;
  counter=0;

  //trimm end
  while(edge_quality<min_edge_qual){
    if(itrk-counter<=itrstart){ //no area in candidate has the required minimum score/base
      candidate=std::make_tuple(0,0,0,4294967295,true);
      return;
    }
    edge_quality+=lookQual[ABU(itrk-counter)]*ACT(itrk-counter);
    while(POS(itrk)-POS(itrk-counter)>edge_len){ //shrinking window until edge_len criterion met
      edge_quality-=lookQual[ABU(itrk)]*ACT(itrk);
      window_quality-=lookQual[ABU(itrk)]*ACT(itrk);
      itrk--;
      counter--;
    }
    counter++;
  }

  float orderScore;
  int32_t orderCount=0;
  std::vector<std::tuple<uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>>::const_iterator itrOrd=itrstart;
  uint32_t old_ord=ORD(itrOrd);
  itrOrd++;
  for(itrOrd; itrOrd<=itrk; itrOrd++){
    orderCount+=1-2*(old_ord<ORD(itrOrd));
    old_ord=ORD(itrOrd);
  }
  orderScore=(float)std::abs(orderCount)/(float)std::distance(itrstart, itrk);
  // orderScore=(float)orderCount/(float)std::distance(itrstart, itrk);
  window_quality*=orderScore;
  bool reverse=true;
  if(orderCount<0){
    reverse=false;
  }

  candidate=std::make_tuple(window_quality,REF(itrstart),POS(itrstart),POS(itrk),reverse);
  return;
}

void AppendPosLong(std::vector<std::tuple <uint8_t,uint32_t,uint32_t,uint32_t,uint32_t>> & kmer_list, const int64_t & hash, openAddressingKmerHashtable & Index, uint8_t & minimizer_active_bases, uint32_t order, longmapOptions & options){
      uint32_t c=Index.GetBkt(hash);
      uint32_t abundance=Index.dir[c+1]-Index.dir[c];
      if (abundance<=options.max_abundance){
        kmer_list.reserve(kmer_list.size()+abundance);
        for (uint32_t i = Index.dir[c];i!=Index.dir[c+1];i++){
          kmer_list.push_back(std::make_tuple(Index.ref[i],Index.pos[i],abundance,minimizer_active_bases,order));
        }
      }
      return;
}
