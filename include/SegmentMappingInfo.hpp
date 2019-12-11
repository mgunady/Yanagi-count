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

#ifndef __RAPMAP_SEGMENT_MAPPING_INFO_HPP__
#define __RAPMAP_SEGMENT_MAPPING_INFO_HPP__

#include <vector>
#include <map>
#include <utility>
#include "nonstd/span.hpp"
#include "cuckoohash_map.hh"
#include "spdlog/spdlog.h"
#include "metro/metrohash64.h"

using SegmentIDType = uint32_t;

struct SegmentCountValue {
  std::array<uint32_t, 8> typeCounts{{0, 0, 0, 0, 0, 0, 0, 0}};
};


// Parse segment information from yanagi (https://github.com/mgunady/yanagi)
class SegmentMappingInfo {
  using TranscriptIDType = int32_t;
  using GeneIDType = int32_t;
public:
  SegmentMappingInfo();
  bool loadFromFile(const std::string& fname,
                    const std::vector<std::string>& segmentNames,
                    std::shared_ptr<spdlog::logger> log);

  void serialize(const std::string& outDir);

  void load(const std::string& outDir);

  size_t numSegments() const { return txpListRanges_.size(); }

  std::string getSegID(int64_t sID) const { return segNames[sID];}

  nonstd::span<TranscriptIDType> transcriptsForSegment(int64_t segmentID) const {
    auto& p = txpListRanges_[segmentID];
    return nonstd::span<TranscriptIDType>(static_cast<int32_t*>(const_cast<int32_t*>(txpList_.data()) + p.first), static_cast<size_t>(p.second));
  }
  
	/**
   * modified from : https://en.cppreference.com/w/cpp/algorithm/set_intersection
   * return true of the two sorted ranges have a non-null intersection, false otherwise.
   **/
  template<class InputIt1, class InputIt2, class OutputIt>
  int32_t getIntersection(InputIt1 first1, InputIt1 last1,
					   InputIt2 first2, InputIt2 last2, OutputIt d_first) {
	int32_t len = 0;
	while (first1 != last1 && first2 != last2) {
	  if (*first1 < *first2) {
		++first1;
	  } else  {
		if (!(*first2 < *first1)) {
		  *d_first++ = *first1++;
		  len ++;
		}
		++first2;
	  }
	}
	return len;
  }
  
  int32_t transcriptsForPair(SegmentIDType s1, SegmentIDType s2, std::vector<TranscriptIDType>& txp_intersect) {
	auto p = std::make_pair(s1, s2);
	//nonstd::span<TranscriptIDType> txpList;
	//std::vector<int> txp_intersect;
	if(pairTxMap_.contains(p)) {
		txp_intersect = pairTxMap_.find(p);
	} else {
		auto txps1 = transcriptsForSegment(s1);
		if(s1 == s2) {
			txp_intersect = std::vector<TranscriptIDType>(txps1.begin(), txps1.end());
		} else {
			auto txps2 = transcriptsForSegment(s2);
			txp_intersect.reserve(std::min(txps1.size(), txps2.size()));
			getIntersection(txps1.begin(), txps1.end(),
							  txps2.begin(), txps2.end(),
							  std::back_inserter(txp_intersect));
			//txpList = nonstd::span<TranscriptIDType>(static_cast<int32_t*>(const_cast<int32_t*>(txp_intersect.data())), static_cast<size_t>(interLen));
		}
		pairTxMap_.insert(p, txp_intersect);
	}
    return txp_intersect.size();
  }

  size_t tableSize() const { return countMap_.size(); }

  void addHit(SegmentIDType s1, SegmentIDType s2, uint8_t mappingType) {
    auto p = std::make_pair(s1, s2);
    auto upfn = [mappingType](SegmentCountValue& x) -> void {
      x.typeCounts[mappingType] += 1;
    };
    SegmentCountValue v; v.typeCounts[mappingType] = 1;
    countMap_.upsert(p, upfn, v);
  }

  void addMultimap(std::string segMultiClass) { 
    auto upfn = [](SegmentCountValue& x) -> void { x.typeCounts[0] += 1; };
    SegmentCountValue v; v.typeCounts[0] = 1;
    multimaps_.upsert(segMultiClass, upfn, v);  
  }

  void addJPHit(SegmentIDType s1, SegmentIDType s2, int32_t pos1, int32_t pos2) {
    auto p = std::make_pair(s1, s2);
    auto upfn = [pos1, pos2](SegmentCountValue& x) -> void { 
       auto sumpos1 = x.typeCounts[0]*x.typeCounts[1] + pos1;
       auto sumpos2 = x.typeCounts[0]*x.typeCounts[2] + pos2;
       x.typeCounts[0] += 1;
       x.typeCounts[1] = sumpos1/x.typeCounts[0];
       x.typeCounts[2] = sumpos2/x.typeCounts[0]; 
    };
    SegmentCountValue v; v.typeCounts[0] = 1; v.typeCounts[1] = pos1; v.typeCounts[2] = pos2;
    newJuncsJPCountMap_.upsert(p, upfn, v);
  }


  bool writeSegmentOutput(const std::string& segFile, const std::vector<std::string>& segNames);

  uint32_t getSegTxPosition(int32_t segID, int32_t txID) const { return segTxPos_.find(std::make_pair(segID, txID))->second; }
  
  int32_t getGeneOfSeg(int32_t sid) const { return genesList_[sid];}

private:
    struct pairhash {
    public:
      template <typename T>
      std::size_t operator()(const std::pair<T, T> &x) const
      {
        T d[2] = {x.first, x.second};
        uint64_t hashKey{0};
        MetroHash64::Hash(reinterpret_cast<uint8_t*>(d), 2*sizeof(T), reinterpret_cast<uint8_t*>(&hashKey), 0);
        return hashKey;
      }
    };

    // The segment id provides an
    // interval into the transcript list vector of
    // all transcripts corresponding to this segment.
    std::vector<std::pair<uint32_t, uint32_t>> txpListRanges_;
    std::vector<TranscriptIDType> txpList_;
    std::vector<std::string> txpNames_;

    std::vector<GeneIDType> genesList_;
    std::vector<std::string> geneNames_;

    std::vector<std::string> segNames;

  std::map<std::pair<SegmentIDType, TranscriptIDType>, uint32_t> segTxPos_;
  cuckoohash_map<std::pair<SegmentIDType, SegmentIDType>, SegmentCountValue, pairhash> countMap_;
  cuckoohash_map<std::pair<SegmentIDType, SegmentIDType>, std::vector<TranscriptIDType>, pairhash> pairTxMap_;
  
  cuckoohash_map<std::string, SegmentCountValue> multimaps_;
  cuckoohash_map<std::pair<SegmentIDType, SegmentIDType>, SegmentCountValue, pairhash> newJuncsJPCountMap_;

};

#endif // __RAPMAP_SEGMENT_MAPPING_INFO_HPP__
