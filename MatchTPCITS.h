// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file MatchTPCITS.h
/// \brief Class to perform TPC ITS matching
/// \author ruben.shahoyan@cern.ch

#ifndef ALICEO2_GLOBTRACKING_MATCHTPCITS_
#define ALICEO2_GLOBTRACKING_MATCHTPCITS_

#include <Rtypes.h>
#include <array>
#include <vector>
#include <string>
#include <TStopwatch.h>
#include "DetectorsBase/Track.h"

class TChain;

namespace o2
{

namespace ITS {  class CookedTrack; }
namespace TPC {  class TrackTPC; }
 
 class MCCompLabel;
 
namespace globaltracking
{

///< flags to tell the status of TPC-ITS tracks comparison
enum TrackRejFlag : int {
  Accept = 0,
  RejectOnY,   // rejected comparing DY difference of tracks
  RejectOnZ,
  RejectOnSnp,
  RejectOnTgl,
  RejectOnQ2Pt,
  RejectOnChi2,
  NSigmaShift = 10  
};
  
///< TPC track parameters propagated to reference X, with time bracket and index of
///< original track in the currently loaded TPC reco output
struct TrackLocTPC {
  o2::track::TrackParCov track;
  int trOrigID = -1;   ///< original index of the TPC track in the input packet (tree entry)
  float timeMin = 0.f; ///< min. possible time (in TPC time-bin units)
  float timeMax = 0.f; ///< max. possible time (in TPC time-bin units)
  TrackLocTPC(const o2::track::TrackParCov& src, int trid) : track(src),trOrigID(trid) {}
  TrackLocTPC() = default;
  ClassDefNV(TrackLocTPC,1);
};

///< ITS track outward parameters propagated to reference X, with time bracket and index of
///< original track in the currently loaded ITS reco output
struct TrackLocITS {
  o2::track::TrackParCov track;
  int trOrigID = -1;    ///< original index of the ITS track in the input packet (tree entry)
  int eventID = -1;     ///< packet (tree entry) this track belongs to
  int roFrame = -1;     ///< ITS readout frame assigned to this track
  float timeMin = 0.f;  ///< min. possible time (in TPC time-bin units)
  float timeMax = 0.f;  ///< max. possible time (in TPC time-bin units)
  TrackLocITS(const o2::track::TrackParCov& src, int trid, int evid) : track(src),trOrigID(trid),eventID(evid) {}
  TrackLocITS() = default;
  ClassDefNV(TrackLocITS,1);
};

///< record of single ITS track matching to given TPC track and reference on
///< the match record of the same TPC track 
struct MatchRecord {
  static const int Dummy; ///< flag dummy index, must be negative
  float chi2 = -1.f;           ///< matching chi2
  int trOrigID = Dummy;        ///< original index of the track in its packet (tree entry)
  int eventID = Dummy;         ///< packet (tree entry) the track belongs to
  int nextRecID = Dummy;       ///< index of eventual next record
  MatchRecord(int trid, int evid, float chi2match) :
    trOrigID(trid), eventID(evid), chi2(chi2match) {}
  MatchRecord(int trid, int evid, float chi2match, int nxt) :
    trOrigID(trid), eventID(evid), chi2(chi2match), nextRecID(nxt) {}
  MatchRecord() = default;
  ClassDefNV(MatchRecord,1);
};

class MatchTPCITS {
  
 public:

  ///< perform matching for provided input
  void run();
  
  ///< perform all initializations
  void init();

  ///< set ITS ROFrame duration in microseconds
  void setITSROFrame(float fums) { mITSROFrame = fums; }

  ///< get number of matching records for TPC track
  int getNMatchRecords(int tpcTrackID) const;


#ifdef _ALLOW_DEBUG_TREES_
  enum DebugFlagTypes : UInt_t {
    MatchTreeAll      = 0x1<<1,    ///< produce matching candidates tree for all candidates
    MatchTreeAccOnly  = 0x1<<2     ///< fill the matching candidates tree only once the cut is passed
  };
  ///< check if partucular flags are set
  bool isDebugFlag(UInt_t flags) const { return mDBGFlags & flags; }

  ///< get debug trees flags
  UInt_t getDebugFlags() const { return mDBGFlags; }

  ///< set or unset debug stream flag
  void setDebugFlag(UInt_t flag, bool on=true);

  ///< fill matching debug tree
  void fillITSTPCmatchTree(int itsID,int tpcID, int rejFlag, float chi2=-1.);
#endif
  
 private:
  
  void attachInputChains();
  bool prepareTPCData();
  bool prepareITSData();
  bool loadTPCData();
  bool loadITSData();
  void doMatching(int sec);
  int compareITSTPCTracks(const TrackLocITS& tITS,const TrackLocTPC& tTPC, float& chi2) const;
  bool registerMatchRecord(const TrackLocITS& tITS,const TrackLocTPC& tTPC, float& chi2);
  float getPredictedChi2NoZ(const o2::track::TrackParCov& tr1, const o2::track::TrackParCov& tr2) const;
  bool propagateToRefX(o2::track::TrackParCov &trc);
  void addTrackCloneForNeighbourSector(const TrackLocITS& src, int sector);
  void printCandidates() const; // temporary

  
  ///< convert TPC time bin to ITS ROFrame units
  int tpcTimeBin2ITSROFrame(float tbin) const {
    int rof = tbin*mTPCBin2ITSROFrame;
    return rof<0 ? 0 : rof;
  }

  ///< convert ITS ROFrame to TPC time bin units
  float itsROFrame2TPCTimeBin(int rof) const {
    return rof*mITSROFrame2TPCBin;
  }

  ///< convert Z interval to TPC time-bins
  float z2TPCBin(float z) const {
    return z*mZ2TPCBin;
  }

  ///< convert TPC time-bins to Z interval
  float tpcBin2Z(float t) const {
    return t*mTPCBin2Z;
  }

  ///< rought check of 2 track params difference, return -1,0,1 if it is <,within or > than tolerance
  int roughCheckDif(float delta, float toler, int rejFlag) const {
    return delta>toler ? rejFlag : (delta<-toler ? -rejFlag : Accept);
  }

  bool mInitDone = false;   ///< flag init already done
  
  int mCurrTPCTreeEntry=-1; ///< current TPC tree entry loaded to memory
  int mCurrITSTreeEntry=-1; ///< current ITS tree entry loaded to memory
  float mXTPCInnerRef = 80.0; ///< reference radius at which TPC provides the tracks 
  float mXRef = 70.0;         ///< reference radius to propage tracks for matching
  float mYMaxAtXRef = 0.;     ///< max Y in the sector at reference X

  bool mMCTruthON = false;   ///< flag availability of MC truth

  ///< do we use track Z difference to reject fake matches? makes sense for triggered mode only
  bool mCompareTracksDZ = false;
  
  ///<tolerance on abs. different of ITS/TPC params
  std::array<float,o2::track::kNParams> mCrudeAbsDiff = {2.f, 2.f, 0.2f, 0.2f, 4.f};
  ///<tolerance on per-component ITS/TPC params NSigma
  std::array<float,o2::track::kNParams> mCrudeNSigma = {49.f,49.f,49.f,49.f,49.f};
  float mCutChi2TPCITS = 200.f;  ///< cut on matching chi2
  float mSectEdgeMargin2 = 0.;   ///< crude check if ITS track should be matched also in neighbouring sector

  int mMaxMatchCandidates = 5; ///< max allowed matching candidates per TPC track

  
  ///< safety margin (in TPC time bins) for ITS-TPC tracks time (in TPC time bins!) comparison 
  float mITSTPCTimeBinSafeMargin = 1.f; 

  ///< safety margin in cm when estimating TPC track tMin and tMax from assigned time0 and its
  ///< track Z position
  float mTPCTimeEdgeZSafeMargin = 20.f;
  
  ///< safety margin in TPC time bins when estimating TPC track tMin and tMax from
  ///< assigned time0 and its track Z position (converted from mTPCTimeEdgeZSafeMargin)
  float mTPCTimeEdgeTSafeMargin = 0.f;

  float mTimeBinTolerance = 10.; ///<tolerance in time-bin for ITS-TPC time bracket matching

  float mITSROFrame = -1.; ///< ITS RO frame in \mus
  float mTPCVDrift0 = -1.; ///< TPC nominal drift speed in cm/microseconds
  float mITSROFrame2TPCBin = 0.; ///< conversion coeff from ITS ROFrame units to TPC time-bin
  float mTPCBin2ITSROFrame = 0.; ///< conversion coeff from TPC time-bin to ITS ROFrame units
  float mZ2TPCBin = 0.; ///< conversion coeff from Z to TPC time-bin
  float mTPCBin2Z = 0.; ///< conversion coeff from TPC time-bin to Z
  float mNTPCBinsFullDrift = 0.; ///< max time bin for full drift
  float mTPCZMax = 0.;

  TChain *mChainITS=nullptr; ///< input chain for ITS tracks
  TChain *mChainTPC=nullptr; ///< input chain for TPC tracks
 
  ///>>>------ these are input arrays which should not be modified by the matching code
  //           since this info is provided by external device
  std::vector<o2::ITS::CookedTrack> *mITSTracksArrayInp = nullptr;             ///<input tracks
  std::vector<o2::TPC::TrackTPC> *mTPCTracksArrayInp = nullptr;                ///<input tracks

  o2::dataformats::MCTruthContainer<o2::MCCompLabel> *mITSTrkLabels = nullptr; ///< input ITS Track MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> *mTPCTrkLabels = nullptr; ///< input TPC Track MC labels
  /// <<<-----
  
  std::vector<int> mMatchRecordID; ///< refs on 1st matchRecord in mMatchRecords of each TPC track
  std::vector<MatchRecord> mMatchRecords; ///< match records pool

  std::vector<TrackLocTPC> mTPCWork; ///<TPC track params prepared for matching
  std::vector<TrackLocITS> mITSWork; ///<ITS track params prepared for matching
  std::vector<o2::MCCompLabel> mTPCLblWork; ///<TPC track labels
  std::vector<o2::MCCompLabel> mITSLblWork; ///<ITS track labels

  ///< per sector indices of TPC track entry in mTPCWork 
  std::array<std::vector<int>,o2::constants::math::NSectors> mTPCSectIndexCache;
  ///< per sector indices of ITS track entry in mITSWork 
  std::array<std::vector<int>,o2::constants::math::NSectors> mITSSectIndexCache;

  ///<indices of 1st entries with time-bin above the value
  std::array<std::vector<int>,o2::constants::math::NSectors> mTPCTimeBinStart;
  ///<indices of 1st entries of ITS tracks with givem ROframe
  std::array<std::vector<int>,o2::constants::math::NSectors> mITSTimeBinStart; 

  std::string mITSTrackBranchName = "ITSTrack"; ///< name of branch containing input ITS tracks
  std::string mTPCTrackBranchName = "Tracks";    ///< name of branch containing input TPC tracks
  std::string mITSMCTruthBranchName = "ITSTrackMCTruth"; ///< name of branch containing ITS MC labels
  std::string mTPCMCTruthBranchName = "TracksMCTruth"; ///< name of branch containing input TPC tracks
  
#ifdef _ALLOW_DEBUG_TREES_
  std::unique_ptr<TreeStreamRedirector> mDBGOut;
  UInt_t mDBGFlags = 0;
#endif
 

  ///----------- aux stuff --------------///
  static constexpr float TolerSortTime = 0.1; ///<tolerance for comparison of 2 tracks times
  static constexpr float TolerSortTgl = 1e-4; ///<tolerance for comparison of 2 tracks tgl
  static constexpr float Tan70 = 2.74747771e+00; // tg(70 degree): std::tan(70.*o2::constants::math::PI/180.);
  static constexpr float Cos70I2 = 1.+Tan70*Tan70; // 1/cos^2(70) = 1 + tan^2(70)
  static constexpr float MaxSnp = 0.85;
  static constexpr float MaxTgp = 1.61357f; // max tg corresponting to MaxSnp = MaxSnp/std::sqrt((1.-MaxSnp)*(1.+MaxSnp));

  TStopwatch mTimer;
  
  ClassDefNV(MatchTPCITS,1);
};
 
}
}


#endif
