#ifndef COVERAGE_ANALYSIS_H
#define COVERAGE_ANALYSIS_H

# include <map>
# include <vector>
# include "functions.h"
# include "parser.h"

void coverageAnalysis(std::vector<result_t> & barcodeMappings, std::vector<uint32_t> & histogram, bcmapOptions & options);
std::map<std::string, std::pair<uint32_t, uint32_t> > getContigPositions(std::vector<result_t> & Barcode_Mappings);
std::vector<std::vector<uint32_t>> getCoverageProfile(std::vector<result_t> & Barcode_Mappings, std::string contig, std::vector<uint16_t> & scoreThresholds, std::map<std::string, std::pair<uint32_t, uint32_t>> & contigPositions, uint32_t max_mapping_size);
std::vector<uint16_t> getScoreProfile(std::vector<std::vector<uint32_t>> & coverageProfile, std::vector<uint16_t> & scoreThresholds, uint32_t median);
uint32_t getCoverageMedian(std::vector<std::vector<uint32_t>> coverageProfile);

class CoverageProfiles{
  public:
    std::map<std::string, std::vector<std::vector<uint32_t>> > coverageMap; // contains coverageProfiles for every contig
    std::map<std::string, std::vector<uint16_t>> scoreMap; // contains score Profiles for every contig, the score threshold provided marks the breakof for optimal coverage at this position
    std::map<std::string, std::pair<uint32_t, uint32_t> > contigPositions; // positions of barcode maooings for each contig within the barcodeMappings vector

    CoverageProfiles(std::vector<result_t> & Barcode_Mappings){ //constuctor
      contigPositions = getContigPositions(Barcode_Mappings);
    }

    void addCoverageProfile(std::vector<result_t> & Barcode_Mappings, std::string contig, std::vector<uint16_t> & scoreThresholds, uint32_t max_mapping_size){ // ever profile contains the coverage for all reference widnows of size 1000 for every score threshold
      coverageMap.insert({contig,getCoverageProfile(Barcode_Mappings, contig, scoreThresholds, contigPositions, max_mapping_size)});
      return;
    }

    void addScoreProfile(std::string contig, std::vector<uint16_t> & scoreThresholds){
      uint32_t medianCoverage = getCoverageMedian(coverageMap.at(contig));
      scoreMap.insert({contig,getScoreProfile(coverageMap.at(contig), scoreThresholds, medianCoverage)});
      return;
    }
};


#endif
