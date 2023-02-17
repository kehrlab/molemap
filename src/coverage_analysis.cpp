# include "coverage_analysis.h"
# include "parser.h"
# include "functions.h"


void filterBarcodeMappings(std::vector<result_t> & Barcode_Mappings, CoverageProfiles & CovProf){
  std::vector<result_t> filtered_Mappings;
  filtered_Mappings.reserve(Barcode_Mappings.size());
  // filtered_file.open("filteredMappings.bed");
  uint16_t scoreThreshold;
  // std::string res = "";
  for (std::map<std::string, std::pair<uint32_t, uint32_t> >::iterator itr_contig = CovProf.contigPositions.begin(); itr_contig != CovProf.contigPositions.end(); itr_contig++){
    std::string contig = itr_contig->first;

    std::vector<result_t>::iterator itr_start = Barcode_Mappings.begin() + std::get<0>(CovProf.contigPositions.at(contig));
    std::vector<result_t>::iterator itr_end = Barcode_Mappings.begin() + std::get<1>(CovProf.contigPositions.at(contig));
    for (std::vector<result_t>::iterator itr = itr_start; itr < itr_end; itr++){
      scoreThreshold = CovProf.scoreMap.at(contig)[(*itr).start/1000];
      if((*itr).quality >= scoreThreshold){
        filtered_Mappings.push_back(*itr);
        // res+=(*itr).contig+"\t"+std::to_string((*itr).start)+"\t"+std::to_string((*itr).end)+"\t"+(*itr).barcode+"\t"+std::to_string((*itr).score)+"\n";
      }
    }
    // filtered_file << res;
    // res = "";
  }

  Barcode_Mappings = filtered_Mappings;
  return;
}

std::vector<uint16_t> getScoreProfile(std::vector<std::vector<uint32_t>> & coverageProfile, std::vector<uint16_t> & scoreThresholds, uint32_t median){
  std::vector<uint16_t> scoreThresholdProfile(coverageProfile.size(),0);
  std::vector<uint16_t>::iterator itr_thresh = scoreThresholdProfile.begin();
  std::vector<std::vector<uint32_t>>::iterator itr_cov;
  uint32_t coverage_threshold = (median/10)*8;
  uint16_t size = scoreThresholds.size()-1;
  for (itr_cov = coverageProfile.begin(); itr_cov < coverageProfile.end(); itr_cov++){
    for (int i = size; i != -1; i--){
      if ((*itr_cov)[i] > coverage_threshold){
        *itr_thresh = scoreThresholds[i];
        break;
      }
    }
    itr_thresh++;
  }

  return scoreThresholdProfile;
}

uint32_t getCoverageMedian(std::vector<std::vector<uint32_t>> coverageProfile){
  // calculation of a sparse median for the coverage profile at trustworthyScore
  std::vector<std::vector<uint32_t>>::iterator itr;
  std::vector<uint32_t> sparseCoverage(coverageProfile.size()/10+1,0);
  uint32_t counter = 0;
  uint16_t trustScore = coverageProfile[0].size()-1;
  for (itr = coverageProfile.begin(); itr < coverageProfile.end();itr+=10){
    sparseCoverage[counter] = (*itr)[trustScore];
    counter++;
  }
  std::sort(sparseCoverage.begin(),sparseCoverage.end());
  uint32_t median = sparseCoverage[(int)(sparseCoverage.size()/2)];
  return median;
}

std::vector<std::vector<uint32_t>> getCoverageProfile(std::vector<result_t> & Barcode_Mappings, std::string contig, std::vector<uint16_t> & scoreThresholds, std::map<std::string, std::pair<uint32_t, uint32_t>> & contigPositions, uint32_t max_mapping_size){
  std::vector<result_t>::iterator itr_start = Barcode_Mappings.begin() + std::get<0>(contigPositions.at(contig));
  std::vector<result_t>::iterator itr_end = Barcode_Mappings.begin() + std::get<1>(contigPositions.at(contig));
  uint32_t contigBuckets = ((*(itr_end-1)).end / 1000) + 1 + (max_mapping_size / 1000);
  std::vector<std::vector<uint32_t>> coverageProfile(contigBuckets, std::vector<uint32_t>(scoreThresholds.size(), 0));
  std::vector<std::vector<uint32_t>> startProfile(contigBuckets, std::vector<uint32_t>(scoreThresholds.size(), 0)); //stores how many mappings start in every bucket
  std::vector<std::vector<uint32_t>> endProfile(contigBuckets, std::vector<uint32_t>(scoreThresholds.size(), 0)); //stores how many mappings end in every bucket
  uint32_t start_pos;
  uint32_t end_pos;
  // determine how many barcodes start and end in every reference bucket
  for(itr_start;itr_start < itr_end;itr_start++){ // for every barcode mapping on contig
    for(uint16_t score_index = 0; score_index < scoreThresholds.size(); score_index++){
      if((*itr_start).quality >= scoreThresholds[score_index]){
        start_pos = (*itr_start).start/1000;
        end_pos = (*itr_start).end/1000;
        startProfile[start_pos][score_index]++;
        // std::cerr << "start: " << (*itr_start).start << " end: " << (*itr_start).end << " contig: " << (*itr_start).contig << "\n";
        // std::cerr << __LINE__ << "\t" << end_pos << "\t" << endProfile.size() << "\t" << score_index << "\n";
        endProfile[end_pos][score_index]++;
        // std::cerr << __LINE__ << "\n";
      }else{ // all other score thresholds will also not be met
        break;
      }
    }
  }
  // determine coverage from balance between starting and ending barcodes
  std::vector<uint32_t> starts_sums(scoreThresholds.size(),0);
  std::vector<uint32_t> ends_sums(scoreThresholds.size(),0);

  for(uint32_t i=0; i<coverageProfile.size();i++){
    for(uint16_t s=0; s<scoreThresholds.size();s++){
      starts_sums[s]+=startProfile[i][s];
      ends_sums[s]+=endProfile[i][s];
      // std::cerr << starts[i] << "\t" << ends[i] << "\t" << starts_sum << "," << ends_sum << "\n";
      coverageProfile[i][s]=starts_sums[s]-ends_sums[s];
    }
  }

  return coverageProfile;
}

bool SmallerContig(const result_t & BM1 , const result_t & BM2){
  return BM1.chrom < BM2.chrom;
}

std::map<std::string, std::pair<uint32_t, uint32_t> > getContigPositions(std::vector<result_t> & Barcode_Mappings){
  std::map<std::string, std::pair<uint32_t, uint32_t> > contigPositions;
  std::vector<result_t>::iterator itr_start = Barcode_Mappings.begin();
  std::vector<result_t>::iterator itr_end;
  result_t target;
  target.start=0;

  while(itr_start < Barcode_Mappings.end()){
    target.chrom=(*itr_start).chrom;
    itr_end = std::upper_bound(itr_start, Barcode_Mappings.end(), target, SmallerContig);
    // std::cerr << (*(itr_end-1)).contig << "\t" << (*(itr_end-1)).start << "\t" << (*(itr_end-1)).end << "\t" << (*(itr_end-1)).barcode << "\n";
    contigPositions.insert({(*itr_start).chrom, {itr_start-Barcode_Mappings.begin(), itr_end-Barcode_Mappings.begin()}});
    itr_start = itr_end;
  }
  return contigPositions;
}

std::vector<uint16_t> getScoreThresholds(uint16_t trustworthyScore){
  std::vector<uint16_t> scoreThresholds = {0};
  float num_thresholds = 5;
  for (int i = 1; i<=num_thresholds; i++){
    scoreThresholds.push_back( (uint16_t)((float)((float)trustworthyScore/num_thresholds)*i) );
  }

  return scoreThresholds;
}

std::vector<uint32_t> smooththeHistogram(std::vector<uint32_t> scoreHistogram){
  int len = scoreHistogram.size();
  std::vector<uint32_t> smoothedHistogram(len,0);
  // smoothe first and second value
  smoothedHistogram[0]=(scoreHistogram[0]+scoreHistogram[1]+scoreHistogram[2])/3;
  smoothedHistogram[1]=(scoreHistogram[0]+scoreHistogram[1]+scoreHistogram[2]+scoreHistogram[3])/4;
  // smoothe body of histogram
  for (int i = 2; i<len-2; i++){
    smoothedHistogram[i]=(scoreHistogram[i-2]+scoreHistogram[i-1]+scoreHistogram[i]+scoreHistogram[i+1]+scoreHistogram[i+2])/5;
  }
  // smoothe last and second last value
  smoothedHistogram[len-2]=(scoreHistogram[len-1]+scoreHistogram[len-2]+scoreHistogram[len-3]+scoreHistogram[len-4])/4;
  smoothedHistogram[len-1]=(scoreHistogram[len-1]+scoreHistogram[len-2]+scoreHistogram[len-3])/3;
  return smoothedHistogram;
}

std::vector<std::pair<uint16_t,uint32_t>> getlocalMax(std::vector<uint32_t> & scoreHist){
  std::vector<std::pair<uint16_t,uint32_t>> localMax;
  if(scoreHist[0]>scoreHist[1]){
    localMax.push_back(std::make_pair(0,scoreHist[0]));
  }
  for(int i = 1; i < scoreHist.size()-1; i++){
    if (scoreHist[i] > scoreHist[i-1] && scoreHist[i] > scoreHist[i+1]){
      localMax.push_back(std::make_pair(i,scoreHist[i]));
    }
  }
  if(scoreHist[scoreHist.size()-1] > scoreHist[scoreHist.size()-2]){
    localMax.push_back(std::make_pair(scoreHist.size()-1,scoreHist[scoreHist.size()-1]));
  }
  return localMax;
}

std::vector<std::pair<uint16_t,uint32_t>> getlocalMin(std::vector<uint32_t> & scoreHist){
  std::vector<std::pair<uint16_t,uint32_t>> localMin;
  if(scoreHist[0]<scoreHist[1]){
    localMin.push_back(std::make_pair(0,scoreHist[0]));
  }
  for(int i = 1; i < scoreHist.size()-1; i++){
    if (scoreHist[i] < scoreHist[i-1] && scoreHist[i] < scoreHist[i+1]){
      localMin.push_back(std::make_pair(i,scoreHist[i]));
    }
  }
  if(scoreHist[scoreHist.size()-1] < scoreHist[scoreHist.size()-2]){
    localMin.push_back(std::make_pair(scoreHist.size()-1,scoreHist[scoreHist.size()-1]));
  }
  return localMin;
}

uint16_t getTrustworthyScore(std::vector<std::pair<uint16_t,uint32_t>> locMax, std::vector<std::pair<uint16_t,uint32_t>> locMin){

  uint16_t start;
  uint16_t end;
  std::pair<uint16_t,uint32_t> max1 = {0,0};
  std::pair<uint16_t,uint32_t> max2 = {0,0};
  std::pair<uint16_t,uint32_t> min = {0,4294967295};
  // determine 2 biggeset local maxima
  for(int i = 0; i < locMax.size(); i++){
    if (std::get<1>(locMax[i])>std::get<1>(max1)){
      max2=max1;
      max1=locMax[i];
    } else if(std::get<1>(locMax[i])>std::get<1>(max2)){
      max2=locMax[i];
    }
  }
  // constrain position of interest to inbetween the 2 biggest maxima
  if(std::get<0>(max1) > std::get<0>(max2)){
    start = std::get<0>(max2);
    end = std::get<0>(max1);
  }else{
    start = std::get<0>(max1);
    end = std::get<0>(max2);
  }
  // choose minimum between 2 biggest maxima
  for(int i = 0; i < locMin.size(); i++){
    if(std::get<0>(locMin[i]) < end){
      if(std::get<0>(locMin[i]) > start && std::get<1>(locMin[i]) < std::get<1>(min)){
        min = locMin[i];
      }
    }else{
      break;
    }
  }

  return std::get<0>(min);
}

uint16_t computeTrustworthyScore(std::vector<uint32_t> scoreHistogram){
  std::vector<std::pair<uint16_t,uint32_t>> locMax = getlocalMax(scoreHistogram);
  std::vector<std::pair<uint16_t,uint32_t>> locMin = getlocalMin(scoreHistogram);
  uint16_t trustworthyScore = getTrustworthyScore(locMax, locMin);
  return trustworthyScore;
}

int parseScoreThreshold(uint32_t trustworthyScore, std::vector<uint32_t> hist){
  if(trustworthyScore != 0){
    return trustworthyScore;
  }else{
    hist = smooththeHistogram(hist);
    trustworthyScore = computeTrustworthyScore(hist);
    return trustworthyScore;
  }
}

void coverageAnalysis(std::vector<result_t> & barcodeMappings, std::vector<uint32_t> & histogram, mapOptions & options){
  uint32_t trustworthyScore = parseScoreThreshold(options.s, histogram);
  std::vector<uint16_t> scoreThresholds = getScoreThresholds(trustworthyScore);

  CoverageProfiles CovProf(barcodeMappings);

  for (std::map<std::string, std::pair<uint32_t, uint32_t> >::iterator itr_contig = CovProf.contigPositions.begin(); itr_contig != CovProf.contigPositions.end(); itr_contig++){
    CovProf.addCoverageProfile(barcodeMappings, itr_contig->first, scoreThresholds, options.max_window_size);
  }

  for (std::map<std::string, std::pair<uint32_t, uint32_t> >::iterator itr_contig = CovProf.contigPositions.begin(); itr_contig != CovProf.contigPositions.end(); itr_contig++){
    CovProf.addScoreProfile(itr_contig->first, scoreThresholds);
  }

  filterBarcodeMappings(barcodeMappings, CovProf);

  return;
}
