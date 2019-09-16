//
// RapMap - Rapid and accurate mapping of short reads to transcriptomes using
// quasi-mapping.
// Copyright (C) 2015-2018 Rob Patro, Avi Srivastava, Hirak Sarkar
//
// This file is part of RapMap.
//
// RapMap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RapMap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RapMap.  If not, see <http://www.gnu.org/licenses/>.
//

#include <fstream>
#include <map>
#include <utility>
#include <sstream>
#include "sparsepp/spp.h"
#include "SegmentMappingInfo.hpp"
#include "RapMapUtils.hpp"

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>


SegmentMappingInfo::SegmentMappingInfo() {}

bool SegmentMappingInfo::loadFromFile(const std::string& fname,
                                      const std::vector<std::string>& segmentNames,
                                      std::shared_ptr<spdlog::logger> log
                                      ) {
  std::ifstream ifile(fname);
  std::string line;
  // The first line is the header
  std::getline(ifile, line);
  auto header = rapmap::utils::tokenize(line, '\t');
  std::map<std::string, int32_t> hmap;
  for (size_t i = 0; i < header.size(); ++i) {
    hmap[header[i]] = i;
  }

  spp::sparse_hash_map<std::string, uint32_t> segmentID;
  segmentID.reserve(segmentNames.size());

  spp::sparse_hash_map<std::string, uint32_t> txpID;
  spp::sparse_hash_map<std::string, uint32_t> geneNameIDMap;

  std::vector<bool> seenSegment(segmentNames.size(), false);

  uint32_t ctr{0};
  for (const auto& sn : segmentNames) {
    segmentID[sn] = ctr;
    ++ctr;
  }

  // number of segments we care about
  //txpListRanges_.resize(segmentNames.size());

  int32_t segID = hmap["segID"];
  int32_t txAnnIDs = hmap["txAnnIDs"];
  int32_t posInTx = hmap["posInTxs"];
  int32_t geneID = hmap["geneID"];

  std::vector<std::pair<uint32_t, uint32_t>> segmentIntervals;

  uint64_t lastIdx{0};
  while (std::getline(ifile, line)) {
    auto tokens = rapmap::utils::tokenize(line, '\t');
    std::string sid = tokens[segID];
    std::string txlist = tokens[txAnnIDs];
    std::string posInTxlist = tokens[posInTx];
    std::string gid = tokens[geneID];

    auto sidIt = segmentID.find(sid);
    uint32_t segIdx{std::numeric_limits<uint32_t>::max()};
    if (sidIt == segmentID.end()) {
      //log->error("Segment mapping {} referred to a segment {} not in the sequence file.", fname, sid);
      continue;
      //return false;
    } else {
      segIdx = sidIt->second;
    }

    if (seenSegment[segIdx]) {
      log->error("Segment {} was seen more than once in the segment mapping {}.", sidIt->first, fname);
      return false;
    }
    seenSegment[segIdx] = true;

    auto txtokens = rapmap::utils::tokenize(txlist, ',');
    auto posInTxtokens = rapmap::utils::tokenize(posInTxlist, ',');
    auto firstIdx = lastIdx;
    uint32_t len{0};
    auto posInTx =  posInTxtokens.begin();
    for (auto& tx : txtokens) {
      // index for this transcript
      auto txIt = txpID.find(tx);
      uint32_t txpIdx{std::numeric_limits<uint32_t>::max()};
      if (txIt == txpID.end()) {
        txpIdx = txpID.size();
        txpID[tx] = txpIdx;
        txpNames_.push_back(tx);
      } else {
        txpIdx = txIt->second;
      }
      txpList_.push_back(txpIdx);
      segTxPos_[std::make_pair(segIdx, txpIdx)] = std::stoi(*posInTx);
      ++lastIdx;
      ++len;
      ++posInTx;
    }
    txpListRanges_.push_back(std::make_pair(firstIdx, len));
    auto& r = txpListRanges_.back();
    std::sort(txpList_.begin() + r.first, txpList_.begin() + r.first + r.second);
    //std::cerr << "sid : " << sid << " :: txlist : " << txlist << "\n";

    auto geneIt = geneNameIDMap.find(gid);
    uint32_t geneIdx{std::numeric_limits<uint32_t>::max()};
    if (geneIt == geneNameIDMap.end()) {
      geneIdx = geneNameIDMap.size();
      geneNameIDMap[gid] = geneIdx;
      geneNames_.push_back(gid);
    } else {
      geneIdx = geneIt->second;
    }
    genesList_.push_back(geneIdx);
  }

  size_t siIdx{0};
  /*
  txpListRanges_.resize(segmentIntervals.size());
  for (auto& si : segmentIntervals) {
    txpListRanges_[siIdx] = nonstd::span<int32_t>(txpList_.data() + si.first, si.second);
    std::sort(txpListRanges_[siIdx].begin(), txpListRanges_[siIdx].end());
    ++siIdx;
  }
  */

  log->info("Total txps {} for {} segments.", txpNames_.size(), txpListRanges_.size());
  log->info("Total txps list length = {}", txpList_.size());

  /*
  siIdx = 0;
  for (auto p : txpListRanges_) {
    std::stringstream s;
    s << "transcripts for " << siIdx << " : {";
    nonstd::span<int32_t> r(txpList_.data() + p.first, p.second);
    for (auto tid : r) {
      s << tid << ',';
    }
    s << "}\n";
    log->info(s.str());
    siIdx++;
  }
  */


  ifile.close();
  return true;
}

void SegmentMappingInfo::serialize(const std::string& outDir) {
  std::ofstream mappingStream(outDir + "segmentMappingInfo.bin", std::ios::binary);
  {
    //ScopedTimer timer;
    cereal::BinaryOutputArchive segArchive(mappingStream);
    segArchive(txpListRanges_);
    segArchive(txpList_);
    segArchive(txpNames_);
    segArchive(genesList_);
    segArchive(geneNames_);
    segArchive(segTxPos_);
  }
  mappingStream.close();
}

bool SegmentMappingInfo::writeSegmentOutput(const std::string& segFile, const std::vector<std::string>& segNames) {
  std::ofstream ofile(segFile+".tsv");
//std::cerr << segNames[410295] << " " << segNames[492716] << "\n";

  // number of segments
  //ofile << segNames.size() << '\n';
  //ofile << txpNames_.size() << '\n';

  auto lt = countMap_.lock_table();

  //ofile << lt.size() << '\n';

  //for (auto& sn : segNames) {
  //  ofile << sn << '\n';
  //}

  for (auto& kv : lt) {
    auto sp = kv.first;
    ofile << segNames[sp.first] << '\t' ;

    if (sp.second < std::numeric_limits<decltype(sp.second)>::max()) {
      ofile << segNames[sp.second];
    } else {
      ofile << "NA";
    }

    for (size_t i = 0; i < 8; ++i) {
      ofile << '\t' << kv.second.typeCounts[i];
    }
    ofile << '\n';
  }

  ofile.close();

  // Write multimappings
  std::ofstream ofile2(segFile+".tsv.multimappings");
  auto lt2 = multimaps_.lock_table();
  for (auto& kv : lt2) {
     ofile2 << kv.first << '\t' << kv.second.typeCounts[0] << '\n' ;
  }
  ofile2.close();
  std::cerr << "\n Found " << multimaps_.size() << " Multimapping segment-Pair Classes (MPCs)\n";

  // Write New Junctions
  std::ofstream ofile3(segFile+".tsv.njpcounts");
  auto lt3 = newJuncsJPCountMap_.lock_table();
  for (auto& kv : lt3) {
    auto sp = kv.first;
    ofile3 << segNames[sp.first] << '\t' ;

    if (sp.second < std::numeric_limits<decltype(sp.second)>::max()) {
      ofile3 << segNames[sp.second];
    } else {
      ofile3 << "NA";
    }

    for (size_t i = 0; i < 8; ++i) {
      ofile3 << '\t' << kv.second.typeCounts[i];
    }
    ofile3 << '\n';
  }
  ofile3.close();


  return true;
}

void SegmentMappingInfo::load(const std::string& indDir) {
  std::ifstream mappingStream(indDir + "segmentMappingInfo.bin", std::ios::binary);
  {
    //ScopedTimer timer;
    cereal::BinaryInputArchive segArchive(mappingStream);
    segArchive(txpListRanges_);
    segArchive(txpList_);
    segArchive(txpNames_);
    segArchive(genesList_);
    segArchive(geneNames_);
    segArchive(segTxPos_);
  }
  mappingStream.close();
}
