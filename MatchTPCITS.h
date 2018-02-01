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

#define _ALLOW_DEBUG_TREES_ // to allow debug and control tree output

#include <Rtypes.h>
#include <array>
#include <vector>
#include <string>
#include <TStopwatch.h>
#include "ReconstructionDataFormats/Track.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "CommonUtils/TreeStreamRedirector.h"

class TChain;

namespace o2
{

namespace dataformats
{
template <typename TruthElement>
class MCTruthContainer;
}
 
namespace ITS
{
class CookedTrack;
}

namespace TPC
{
class TrackTPC;
}
 
namespace globaltracking
{

constexpr int DummyID = -1;
  
enum TimerLevel : int
{ // how often stop/start timer
  TimingOff, TimingTotal, TimingLoadData, TimingFillDebugTree, TimingAll
};

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

///< timing (in TPC time-bins) bracket assumed for the track
struct timeBracket {
  float tmin = 0.f;     ///< min possible time(bin) 
  float tmax = 0.f;     ///< max possible time(bin)
  timeBracket() = default;
  timeBracket(float mn,float mx) : tmin(mn),tmax(mx) {}
  void set(float tmn,float tmx) {
    tmin = tmn;
    tmax = tmx;
  }
};
 
///< identifier for the track entry
struct trackOrigin {
  int id    = DummyID;  ///< entry in the chunk
  int chunk = DummyID;  ///< data chunk (tree entry, message packet..)
  trackOrigin(int tid,int tch) : id(tid), chunk(tch) {}
  trackOrigin() = default;
  void set(int tid,int tch) {
    id = tid;
    chunk = tch;
  }
};
 
///< TPC track parameters propagated to reference X, with time bracket and index of
///< original track in the currently loaded TPC reco output
struct TrackLocTPC {
  o2::track::TrackParCov track;
  trackOrigin source;  ///< track origin id
  timeBracket timeBins;      ///< bracketing time-bins
  int matchID = DummyID;     ///< entry (non if DummyID) of its matchTPC struct in the mMatchesTPC
  TrackLocTPC(const o2::track::TrackParCov& src, int tid, int tch) : track(src),source(tid,tch) {}
  TrackLocTPC() = default;
  ClassDefNV(TrackLocTPC,1);
};

///< ITS track outward parameters propagated to reference X, with time bracket and index of
///< original track in the currently loaded ITS reco output
struct TrackLocITS {
  o2::track::TrackParCov track;
  trackOrigin source;        ///< track origin id
  timeBracket timeBins;      ///< bracketing time-bins
  int roFrame = DummyID;     ///< ITS readout frame assigned to this track
  int matchID = DummyID;     ///< entry (non if DummyID) of its matchITS struct in the mMatchesITS
  TrackLocITS(const o2::track::TrackParCov& src, int tid, int tch) : track(src),source(tid,tch) {}
  TrackLocITS() = default;
  ClassDefNV(TrackLocITS,1);
};

/// RS: at the moment matchTPC and matchITS are similar, but may diverge in future
 
///< each TPC track having at least 1 matching ITS candidate records in matchTPC the
///< the ID of the 1st (best) matchRecord in the mMatchRecordsTPC container
struct matchTPC {
  trackOrigin source;        ///< track origin id
  int first = DummyID;  ///< 1st match for this track in the mMatchRecordsTPC
  matchTPC(const trackOrigin& src) : source(src) {}
  matchTPC() = default;
};

///< each ITS track having at least 1 matching TPC candidate records in matchITS the TODO
///< the ID of its 1st matchLink in the mMatchLinksITS container
struct matchITS {
  trackOrigin source;        ///< track origin id
  int first = DummyID;  ///< 1st match for this track in the mMatchRecordsITS
  matchITS(const trackOrigin& src) : source(src) {}
  matchITS() = default;
};
 
///< record TPC or ITS track associated with single ITS or TPC track and reference on
///< the next (worse chi2) matchRecord of the same TPC or ITS track 
struct matchRecord {
  float chi2 = -1.f;             ///< matching chi2
  int matchID = DummyID;      ///< id of matchITS struct in mMatchesITS container
  int nextRecID = DummyID;       ///< index of eventual next record

  matchRecord(int itsMID, float chi2match) : matchID(itsMID), chi2(chi2match) {}
  matchRecord(int itsMID, float chi2match, int nxt) : matchID(itsMID), chi2(chi2match), nextRecID(nxt) {}
  matchRecord() = default;
};
 
class MatchTPCITS {
  
 public:
  
  ///< perform matching for provided input
  void run();
  
  ///< perform all initializations
  void init();

  ///< set ITS ROFrame duration in microseconds
  void setITSROFrame(float fums) { mITSROFrame = fums; }

  ///< set chain containing ITS tracks 
  void setInputChainITS(TChain* chain) { mChainITS = chain;}

  ///< set chain containing TPC tracks 
  void setInputChainTPC(TChain* chain) { mChainTPC = chain;}

  ///< print settings
  void print() const;
  void printCandidatesTPC() const;
  void printCandidatesITS() const;


  ///< set timing level
  void setTimingLevel(int v) { mTimingLevel = v; }
  ///< get timing level
  bool getTimingLevel() const { return mTimingLevel; }

  //>>> ====================== cuts ================================>>>
  
  ///< set cuts on absolute difference of ITS vs TPC track parameters
  void setCrudeAbsDiffCut(const std::array<float,o2::track::kNParams> &vals) {
    mCrudeAbsDiffCut = vals;
  }
  ///< get cuts on absolute difference of ITS vs TPC track parameters
  const std::array<float,o2::track::kNParams> & getCrudeAbsDiffCut() const {
    return mCrudeAbsDiffCut;
  }

  ///< set cuts on absolute difference of ITS vs TPC track parameters
  void setCrudeNSigmaCut(const std::array<float,o2::track::kNParams> &vals) {
    mCrudeNSigmaCut = vals;
  }
  ///< get cuts on absolute difference of ITS vs TPC track parameters
  const std::array<float,o2::track::kNParams> & getCrudeNSigmaCut() const {
    return mCrudeNSigmaCut;
  }

  ///< set cut matching chi2
  void setCutMatchingChi2(float val) { mCutMatchingChi2 = val; }
  ///< get cut on matching chi2
  float getCutMatchingChi2() const { return mCutMatchingChi2; }

  ///< set max number of matching candidates to consider
  void setMaxMatchCandidates(int n) { mMaxMatchCandidates = n>1 ? 1:n; }
  ///< get max number of matching candidates to consider
  int  getMaxMatchCandidates() const { return mMaxMatchCandidates; }

  ///< set tolerance (TPC time bins) on ITS-TPC times comparison
  void setTimeBinTolerance(float val) {  mTimeBinTolerance = val; }
  ///< get tolerance (TPC time bins) on ITS-TPC times comparison
  float getTimeBinTolerance() const {  return mTimeBinTolerance; }

  ///< set tolerance on TPC time-bins estimate from highest cluster Z
  void setTPCTimeEdgeZSafeMargin(float val) {  mTPCTimeEdgeZSafeMargin = val; }
  ///< get tolerance on TPC time-bins estimate from highest cluster Z
  float getTPCTimeEdgeZSafeMargin() const {  return mTPCTimeEdgeZSafeMargin; }
  
  //<<< ====================== cuts ================================<<<
  
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

  ///< set the name of output debug file
  void setDebugTreeFileName(std::string name) {
    if (!name.empty()) {
      mDebugTreeFileName = name;
    }
  }

  ///< get the name of output debug file
  const std::string& getDebugTreeFileName() const {return mDebugTreeFileName;}
  
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
  float getPredictedChi2NoZ(const o2::track::TrackParCov& tr1, const o2::track::TrackParCov& tr2) const;
  bool propagateToRefX(o2::track::TrackParCov &trc);
  void addTrackCloneForNeighbourSector(const TrackLocITS& src, int sector);

  ///------------------- manipulations with matches records ----------------------
  bool registerMatchRecordTPC(TrackLocITS& tITS, TrackLocTPC& tTPC, float chi2);
  void registerMatchRecordITS(TrackLocITS& tITS, int matchTPCID, float chi2);
  void suppressMatchRecordITS(int matchITSID, int matchTPCID);
  matchTPC& getTPCMatchEntry(TrackLocTPC& tTPC);
  matchITS& getITSMatchEntry(TrackLocITS& tITS);

  ///< get number of matching records for TPC track referring to this matchTPC
  int getNMatchRecordsTPC(const matchTPC& tpcMatch) const;
  
  ///< get number of matching records for ITS track referring to this matchITS
  int getNMatchRecordsITS(const matchITS& itsMatch) const;
  
  ///< get number of matching records for TPC track referring to matchTPS struct with matchTPCID
  int getNMatchRecordsTPC(int matchTPCID) const {
    return matchTPCID==DummyID ? 0 : getNMatchRecordsTPC(mMatchesTPC[matchTPCID]);
  }
  ///< get number of matching records for ITS track referring to matchITS struct with matchITSID
  int getNMatchRecordsITS(int matchITSID) const {    
    return matchITSID==DummyID ? 0 : getNMatchRecordsITS(mMatchesITS[matchITSID]);
  }     
  
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

  ///< start timing if current level timing is allowed
  void timingOn(TimerLevel lev) {
    if (mTimingLevel >= lev) {
      mTimer.Start(false);
    }
  }

  ///< stop timing if current level timing is allowed
  void timingOff(TimerLevel lev) {
    if (mTimingLevel >= lev) {
      mTimer.Stop();
    }
  }
  
  //================================================================
  
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
  std::array<float,o2::track::kNParams> mCrudeAbsDiffCut = {2.f, 2.f, 0.2f, 0.2f, 4.f};
  
  ///<tolerance on per-component ITS/TPC params NSigma
  std::array<float,o2::track::kNParams> mCrudeNSigmaCut = {49.f,49.f,49.f,49.f,49.f};

  float mCutMatchingChi2 = 200.f;  ///< cut on matching chi2

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

  float mTimeBinTolerance = 10.f; ///<tolerance in time-bin for ITS-TPC time bracket matching

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

  ///< container for matchTPC structures (1 per TPCtrack with some matches to ITS)
  std::vector<matchTPC> mMatchesTPC;
  ///< container for matchITS structures (1 per ITStrack with some matches to TPC)
  std::vector<matchITS> mMatchesITS; 

  ///< container for record the match of TPC track to single ITS track
  std::vector<matchRecord> mMatchRecordsTPC; 
  ///< container for reference to matchRecord involving particular ITS track
  std::vector<matchRecord> mMatchRecordsITS; 
  
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
  std::unique_ptr<o2::utils::TreeStreamRedirector> mDBGOut;
  UInt_t mDBGFlags = 0;
  std::string mDebugTreeFileName = "dbg_match.root"; ///< name for the debug tree file
#endif
 

  ///----------- aux stuff --------------///
  static constexpr float TolerSortTime = 0.1; ///<tolerance for comparison of 2 tracks times
  static constexpr float TolerSortTgl = 1e-4; ///<tolerance for comparison of 2 tracks tgl
  static constexpr float Tan70 = 2.74747771e+00; // tg(70 degree): std::tan(70.*o2::constants::math::PI/180.);
  static constexpr float Cos70I2 = 1.+Tan70*Tan70; // 1/cos^2(70) = 1 + tan^2(70)
  static constexpr float MaxSnp = 0.85;
  static constexpr float MaxTgp = 1.61357f; // max tg corresponting to MaxSnp = MaxSnp/std::sqrt((1.-MaxSnp)*(1.+MaxSnp));

  int mTimingLevel = TimingAll;
  TStopwatch mTimer;
  
  ClassDefNV(MatchTPCITS,1);
};

 

//______________________________________________
inline matchTPC& MatchTPCITS::getTPCMatchEntry(TrackLocTPC& tTPC)
{
  ///< return the matchTPC entry referred by the tTPC track,
  ///< create if neaded
  if (tTPC.matchID == DummyID) { // does this TPC track already have any match? If not, create matchTPC entry 
    tTPC.matchID = mMatchesTPC.size();
    mMatchesTPC.emplace_back(tTPC.source);
    return mMatchesTPC.back();
  }
  return mMatchesTPC[tTPC.matchID];
}

//______________________________________________
inline matchITS& MatchTPCITS::getITSMatchEntry(TrackLocITS& tITS)
{
  ///< return the matchITS entry referred by the tITS track,
  ///< create if neaded
  if (tITS.matchID == DummyID) { // does this ITS track already have any match? If not, create matchITS entry 
    tITS.matchID = mMatchesITS.size();
    mMatchesITS.emplace_back(tITS.source);
    return mMatchesITS.back();
  }
  return mMatchesITS[tITS.matchID];
}


 
}
}


#endif
