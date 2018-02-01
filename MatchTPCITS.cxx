// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <TChain.h>
#include <cassert>
#include <TChain.h>

#include "FairLogger.h"
#include "Field/MagneticField.h"
#include "Field/MagFieldFast.h"
#include "ITSReconstruction/CookedTrack.h"
#include "TPCReconstruction/TrackTPC.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

#include "DetectorsBase/Propagator.h"
#include "CommonUtils/TreeStream.h"

#include "TPCBase/Defs.h"
#include "TPCBase/ParameterElectronics.h"
#include "TPCBase/ParameterDetector.h"
#include "TPCBase/ParameterGas.h"
#include "MathUtils/Cartesian3D.h"
#include "MathUtils/Utils.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DetectorsBase/GeometryManager.h"

#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include "DataFormatsParameters/GRPObject.h"

#include "MatchTPCITS.h"

using namespace o2::globaltracking;

using MatrixDSym4 = ROOT::Math::SMatrix<double, 4, 4, ROOT::Math::MatRepSym<double, 4>>;
using MatrixD4 = ROOT::Math::SMatrix<double, 4, 4, ROOT::Math::MatRepStd<double, 4>>;

//______________________________________________
void MatchTPCITS::run()
{
  ///< perform matching for provided input
  if (!mInitDone) {
    LOG(FATAL)<<"init() was not done yet"<<FairLogger::endl;
  }

  timingOn(TimingTotal);
  
  while(prepareTPCData()) {
    while(prepareITSData()) {
      for (int  sec=o2::constants::math::NSectors;sec--;) {
	doMatching(sec);
      }
    }
    printCandidatesTPC();
    printCandidatesITS();
  }

  timingOff(TimingTotal);

  if (mTimingLevel>=TimingTotal) {
    mTimer.Print();
  }
#ifdef _ALLOW_DEBUG_TREES_
  mDBGOut.reset();
#endif
}

//______________________________________________
void MatchTPCITS::init()
{
  ///< perform initizalizations, precalculate what is needed
  if (mInitDone) {
    LOG(ERROR)<<"Initialization was already done"<<FairLogger::endl;
    return;
  }
  mYMaxAtXRef = mXRef*std::tan(o2::constants::math::SectorSpanRad*0.5); ///< max Y in the sector at reference X
  mSectEdgeMargin2 = mCrudeAbsDiffCut[o2::track::kY]*mCrudeAbsDiffCut[o2::track::kY]; ///< precalculated ^2

  const auto & gasParam = o2::TPC::ParameterGas::defaultInstance();
  const auto & elParam = o2::TPC::ParameterElectronics::defaultInstance();
  const auto & detParam = o2::TPC::ParameterDetector::defaultInstance();
  float tpcTBin = elParam.getZBinWidth();
  mTPCVDrift0 = gasParam.getVdrift();
  mTPCZMax = detParam.getTPClength();
  
  assert(mITSROFrame>0.f);

  mITSROFrame2TPCBin = mITSROFrame/tpcTBin;
  mTPCBin2ITSROFrame = 1./mITSROFrame2TPCBin;
  mTPCBin2Z = tpcTBin*mTPCVDrift0;
  mZ2TPCBin = 1./mTPCBin2Z;
  mNTPCBinsFullDrift = mTPCZMax*mZ2TPCBin;

  mTPCTimeEdgeTSafeMargin = z2TPCBin(mTPCTimeEdgeZSafeMargin);

  attachInputChains();
  
#ifdef _ALLOW_DEBUG_TREES_
  // debug streamer
  if (mDBGFlags) {
    mDBGOut = std::make_unique<o2::utils::TreeStreamRedirector>(mDebugTreeFileName.data(),"recreate");
  }
#endif

  mMatchRecordsTPC.clear();
  mMatchRecordsITS.clear();
  mMatchesTPC.clear();
  mMatchesITS.clear();
  
  mInitDone = true;

  print();
  
}


//______________________________________________
int MatchTPCITS::getNMatchRecordsTPC(int matchTPCID) const
{
  ///< get number of matching records for TPC track reffering to matchTPS struct with matchTPCID
  int count = 0, recID = mMatchesTPC[matchTPCID].first;
  while (recID!=DummyID) {
    recID = mMatchRecordsTPC[recID].nextRecID;
    count++;
  }
  return count;
}

//______________________________________________
int MatchTPCITS::getNMatchRecordsITS(int matchITSID) const
{
  ///< get number of matching records for ITS track reffering to matchITS struct with matchITSID
  int count = 0, recID = mMatchesITS[matchITSID].first;
  while (recID!=DummyID) {
    auto& itsRecord = mMatchRecordsITS[recID];
    recID = itsRecord.nextRecID;
    if (itsRecord.matchTPCID!=DummyID) { // some match records might be disabled
      count++;
    }
  }
  return count;
}

//______________________________________________
void MatchTPCITS::attachInputChains()
{
  if (!mChainITS) {
    LOG(FATAL) <<"ITS data input chain is not set"<<FairLogger::endl;
  }
  if (!mChainTPC) {
    LOG(FATAL) <<"TPC data input chain is not set"<<FairLogger::endl;
  }

  if (!mChainITS->GetBranch(mITSTrackBranchName.data())) {
    LOG(FATAL) <<"Did not find ITS tracks branch "<<mITSTrackBranchName<<" in the input chain"<<FairLogger::endl;
  }
  mChainITS->SetBranchAddress(mITSTrackBranchName.data(),&mITSTracksArrayInp);
  LOG(INFO)<<"Attached ITS tracks "<<mITSTrackBranchName<<" branch with "
	   <<mChainITS->GetEntries()<<" entries"<<FairLogger::endl;

  if (!mChainTPC->GetBranch(mTPCTrackBranchName.data())) {
    LOG(FATAL) <<"Did not find TPC tracks branch "<<mTPCTrackBranchName<<" in the input chain"<<FairLogger::endl;
  }
  mChainTPC->SetBranchAddress(mTPCTrackBranchName.data(),&mTPCTracksArrayInp);
  LOG(INFO)<<"Attached TPC tracks "<<mTPCTrackBranchName<<" branch with "
	   <<mChainTPC->GetEntries()<<" entries"<<FairLogger::endl;
  
  // is there MC info available ?
  if (mChainITS->GetBranch(mITSMCTruthBranchName.data())) {
    mChainITS->SetBranchAddress(mITSMCTruthBranchName.data(),&mITSTrkLabels);
    LOG(INFO)<<"Found ITS Track MCLabels branch "<<mITSMCTruthBranchName<<FairLogger::endl;
  }
  // is there MC info available ?
  if (mChainTPC->GetBranch(mTPCMCTruthBranchName.data())) {
    mChainTPC->SetBranchAddress(mTPCMCTruthBranchName.data(),&mTPCTrkLabels);
    LOG(INFO)<<"Found TPC Track MCLabels branch "<<mTPCMCTruthBranchName<<FairLogger::endl;
  }

  mMCTruthON = (mITSTrkLabels && mTPCTrkLabels);  
  mCurrTPCTreeEntry = -1;
  mCurrITSTreeEntry = -1;

}

//______________________________________________
bool MatchTPCITS::prepareTPCData()
{
  ///< load next chunk of TPC data and prepare for matching
  mMatchRecordsTPC.clear();

  if (!loadTPCData()) {
    return false;
  }

  int ntr = mTPCTracksArrayInp->size();
  
  mMatchesTPC.reserve(mMatchesTPC.size() + ntr);
  // number of records might be actually more than N tracks!
  mMatchRecordsTPC.reserve(mMatchRecordsTPC.size() + ntr);

  // copy the track params, propagate to reference X and build sector tables
  mTPCWork.clear();
  mTPCWork.reserve(ntr);
  if (mMCTruthON) {
    mTPCLblWork.clear();
    mTPCLblWork.reserve(ntr);    
  }
  for (int sec=o2::constants::math::NSectors;sec--;) {
    mTPCSectIndexCache[sec].clear();
    mTPCTimeBinStart[sec].clear();
  }
  
  for (int it=0;it<ntr;it++) {
    
    o2::TPC::TrackTPC& trcOrig = (*mTPCTracksArrayInp)[it];

    // make sure the track was propagated to inner TPC radius at the ref. radius
    if (trcOrig.getX()>mXTPCInnerRef+0.1) continue; // failed propagation to inner TPC radius, cannot be matched

    mTPCWork.emplace_back(static_cast<o2::track::TrackParCov&>(trcOrig),it,mCurrTPCTreeEntry); // working copy of track param
    auto & trc = mTPCWork.back();    
    // propagate to matching Xref
    if (!propagateToRefX(trc.track)) {
      mTPCWork.pop_back(); // discard track whose propagation to mXRef failed
      continue;
    }
    if (mMCTruthON) {
      mTPCLblWork.emplace_back(  mTPCTrkLabels->getLabels(it)[0] );
    }
    // max time increment (in time bins) to have last cluster at the Z of endcap
    float dtZEdgeTPC = z2TPCBin( mTPCZMax-fabs(trcOrig.getLastClusterZ()) );
    // max time decrement (in time bins) to have last cluster at the Z=0 (CE)
    float dtZCETPC = z2TPCBin( fabs(trcOrig.getLastClusterZ()) );
    // RS: consider more effective narrowing
    float time0 = trcOrig.getTimeVertex(mTPCBin2Z);
    trc.timeMin = time0 - dtZCETPC - mTPCTimeEdgeTSafeMargin; 
    trc.timeMax = time0 + dtZEdgeTPC + mTPCTimeEdgeTSafeMargin;

    // cache work track index
    mTPCSectIndexCache[o2::utils::Angle2Sector( trc.track.getAlpha() )].push_back( mTPCWork.size()-1 ); 
  }

  // sort tracks in each sector according to their timeMax
  for (int sec=o2::constants::math::NSectors; sec--;) {
    auto &indexCache = mTPCSectIndexCache[sec];
    LOG(INFO) <<"Sorting "<<sec<<" | "<<indexCache.size()<<" TPC tracks"<<FairLogger::endl;
    if (!indexCache.size()) continue;
    std::sort(indexCache.begin(), indexCache.end(),
	      [this](int a, int b) {	       
		auto &trcA = mTPCWork[a];
		auto &trcB = mTPCWork[b];		
		return (trcA.timeMax - trcB.timeMax) < 0.;
	      });

    // build array of 1st entries with tmax corresponding to each ITS RO cycle
    float tmax = mTPCWork[ indexCache.back() ].timeMax;
    int nbins = 1 +  tpcTimeBin2ITSROFrame(tmax);
    auto &tbinStart = mTPCTimeBinStart[sec];
    tbinStart.resize( nbins>1 ? nbins : 1, -1);
    tbinStart[0] = 0;
    for (int itr=0;itr<(int)indexCache.size();itr++) {
      auto &trc = mTPCWork[ indexCache[itr] ];
      int bTrc = tpcTimeBin2ITSROFrame(trc.timeMax);
      if (bTrc<0) {
	continue;
      }
      if (tbinStart[bTrc]==-1) {
	tbinStart[bTrc] = itr;
      }
    }
    for (int i=1;i<nbins;i++) {   
      if (tbinStart[i]==-1) { // fill gaps with preceding indices
	tbinStart[i] = tbinStart[i-1];
      }
    }
  } // loop over tracks of single sector
  
  return true;
}

//_____________________________________________________
bool MatchTPCITS::prepareITSData()
{
  // load next chunk of ITS data and prepare for matching
  
  if (!loadITSData()) {
    return false;
  }

  int ntr = mITSTracksArrayInp->size();
  
  mMatchesITS.reserve(mMatchesITS.size() + ntr);
  // number of records might be actually more than N tracks!
  mMatchRecordsITS.reserve(mMatchRecordsITS.size() + ntr); 

  mITSWork.clear();
  mITSWork.reserve(ntr*1.3); // tracks close to sector edge may require extra copies
  if (mMCTruthON) {
    mITSLblWork.clear();
    mITSLblWork.reserve(ntr*1.3);
  }
  for (int sec=o2::constants::math::NSectors;sec--;) {
    mITSSectIndexCache[sec].clear();
  }

  for (int it=0;it<ntr;it++) {
    auto& trcOrig = (*mITSTracksArrayInp)[it];

    if (trcOrig.getParamOut().getX()<1.) {
      continue; // backward refit failed
    }
    // working copy of outer track param
    mITSWork.emplace_back(static_cast<o2::track::TrackParCov&>(trcOrig.getParamOut()),it,mCurrITSTreeEntry);
    auto & trc = mITSWork.back();

    // TODO: why I did this?
    if ( !trc.track.rotate( o2::utils::Angle2Alpha(trc.track.getPhiPos()) ) ) {
      mITSWork.pop_back(); // discard failed track
      continue;
    } 
    // make sure the track is at the ref. radius
    if (!propagateToRefX(trc.track)) {
      mITSWork.pop_back(); // discard failed track
      continue; // add to cache only those ITS tracks which reached ref.X and have reasonable snp
    }
    if (mMCTruthON) {
      mITSLblWork.emplace_back( mITSTrkLabels->getLabels(it)[0] );
    }

    trc.timeMin = itsROFrame2TPCTimeBin(trcOrig.getROFrame());
    trc.timeMax = trc.timeMin + mITSROFrame2TPCBin;
    trc.roFrame = trcOrig.getROFrame();

    // cache work track index
    int sector = o2::utils::Angle2Sector( trc.track.getAlpha() );
    mITSSectIndexCache[sector].push_back( mITSWork.size()-1 );

    // If the ITS track is very close to the sector edge, it may match also to a TPC track in the neighbouring sector.
    // For a track with Yr and Phir at Xr the distance^2 between the poisition of this track in the neighbouring sector
    // when propagated to Xr (in this neighbouring sector) and the edge will be (neglecting the curvature)
    // [(Xr*tg(10)-Yr)/(tgPhir+tg70)]^2  / cos(70)^2  // for the next sector
    // [(Xr*tg(10)+Yr)/(tgPhir-tg70)]^2  / cos(70)^2  // for the prev sector
    // Distances to the sector edges in neighbourings sectors (at Xref in theit proper frames)
    float tgp = trc.track.getSnp();
    tgp /= std::sqrt((1.f-tgp)*(1.f+tgp)); // tan of track direction XY

    // sector up
    float dy2Up = (mYMaxAtXRef-trc.track.getY())/(tgp + Tan70);
    if ( (dy2Up*dy2Up*Cos70I2)<mSectEdgeMargin2) { // need to check this track for matching in sector up
      addTrackCloneForNeighbourSector(trc, sector<(o2::constants::math::NSectors-1) ? sector+1 : 0);
    }
    // sector down
    float dy2Dn = (mYMaxAtXRef+trc.track.getY())/(tgp - Tan70);
    if ( (dy2Dn*dy2Dn*Cos70I2)<mSectEdgeMargin2) { // need to check this track for matching in sector down
      addTrackCloneForNeighbourSector(trc, sector>1 ? sector-1 : o2::constants::math::NSectors-1); 
    }
  }
  
  // sort tracks in each sector according to their time, then tgl
  for (int sec=o2::constants::math::NSectors; sec--;) {
    auto &indexCache = mITSSectIndexCache[sec];
    LOG(INFO) <<"Sorting "<<sec<<" | "<<indexCache.size()<<" ITS tracks"<<FairLogger::endl;
    if (!indexCache.size()) {
      continue;
    }
    std::sort(indexCache.begin(), indexCache.end(),
	      [this](int a, int b) {
		auto &trackA = mITSWork[a];
		auto &trackB = mITSWork[b];		
		if (trackA.roFrame < trackB.roFrame) { // ITS tracks have the same time coverage
		  return true;
		}
		else if (trackA.roFrame > trackB.roFrame) { 
		  return false;
		}
		return trackA.track.getTgl() < trackB.track.getTgl();
	      });
    
    // build array of 1st entries with of each ITS RO cycle
    int nbins = 1 + mITSWork[ indexCache.back() ].roFrame;
    auto &tbinStart = mITSTimeBinStart[sec];
    tbinStart.resize( nbins>1 ? nbins : 1, -1);
    tbinStart[0] = 0;
    for (int itr=0;itr<(int)indexCache.size();itr++) {
      auto &trc = mITSWork[ indexCache[itr] ];
      if (tbinStart[trc.roFrame]==-1) {
	tbinStart[trc.roFrame] = itr;
      }
    }
    for (int i=1;i<nbins;i++) {   
      if (tbinStart[i]==-1) { // fill gaps with preceding indices
	tbinStart[i] = tbinStart[i-1];
      }
    }

    
  } // loop over tracks of single sector

  return true;
}

//_____________________________________________________
bool MatchTPCITS::loadITSData()
{
  ///< load next chunk of ITS data
  timingOff(TimingLoadData);
  
  while(++mCurrITSTreeEntry < mChainITS->GetEntries()) {
    mChainITS->GetEntry(mCurrITSTreeEntry);
    LOG(INFO)<<"Starting ITS entry "<<mCurrITSTreeEntry<<" -> "
	     <<mITSTracksArrayInp->size()<<" tracks"<<FairLogger::endl;
    if (!mITSTracksArrayInp->size()) {
      continue;
    }
    return true;
    timingOn(TimingLoadData);
  }
  --mCurrITSTreeEntry;
  timingOn(TimingLoadData);
  return false;
}

//_____________________________________________________
bool MatchTPCITS::loadTPCData()
{
  ///< load next chunk of TPC data
  timingOff(TimingLoadData);
  
  while(++mCurrTPCTreeEntry < mChainTPC->GetEntries()) {
    mChainTPC->GetEntry(mCurrTPCTreeEntry);
    LOG(INFO)<<"Starting TPC entry "<<mCurrTPCTreeEntry<<" -> "
	     <<mTPCTracksArrayInp->size()<<" tracks"<<FairLogger::endl;
    if (mTPCTracksArrayInp->size()<1) {
      continue;
    }
    timingOn(TimingLoadData);
    return true;
  }
  --mCurrTPCTreeEntry;

  timingOn(TimingLoadData);
  return false;
}

//_____________________________________________________
void MatchTPCITS::doMatching(int sec)
{
  ///< run matching for currently cached ITS data for given TPC sector
  LOG(INFO)<<"Matching sector "<<sec<<FairLogger::endl;
  
  auto &cacheITS = mITSSectIndexCache[sec]; // array of cached ITS track indices for this sector
  auto &cacheTPC = mTPCSectIndexCache[sec]; // array of cached ITS track indices for this sector
  auto &tbinStartTPC = mTPCTimeBinStart[sec];   // array of 1st TPC track with timeMax in ITS ROFrame
  auto &tbinStartITS = mITSTimeBinStart[sec];
  int nTracksTPC = cacheTPC.size(), nTracksITS = cacheITS.size();
  LOG(INFO)<<"N tracks TPC="<<nTracksTPC<<" ITS="<<nTracksITS<<" in sector "<<sec<<FairLogger::endl;
  if (!nTracksTPC || !nTracksITS) {
    return;
  }
  
  /// full drift time + safety margin
  float maxTDriftSafe = (mNTPCBinsFullDrift+mITSTPCTimeBinSafeMargin+mTPCTimeEdgeTSafeMargin);

  // get min ROFrame (in TPC time-bins) of ITS tracks currently in cache
  auto minROFITS = mITSWork[ cacheITS.front() ].roFrame;

  if (minROFITS>=int(tbinStartTPC.size())) {
    LOG(INFO)<<"ITS min ROFrame "<<minROFITS<<" exceeds all cached TPC track ROF eqiuvalent "<<
	     cacheTPC.size()-1<<FairLogger::endl;
    return;
  }
  
  int nCheckTPCControl = 0, nCheckITSControl = 0, nMatchesControl = 0; // temporary

  int idxMinTPC = tbinStartTPC[minROFITS]; // index of 1st cached TPC track within cached ITS ROFrames
  for (int itpc=idxMinTPC;itpc<nTracksTPC; itpc++) {
    auto & trefTPC = mTPCWork[ cacheTPC[itpc] ];
    // estimate ITS 1st ROframe bin this track may match to: TPC track are sorted according to their
    // timeMax, hence the timeMax - MaxmNTPCBinsFullDrift are non-decreasing
    int itsROBin = tpcTimeBin2ITSROFrame(trefTPC.timeMax - maxTDriftSafe);
    if (itsROBin>=int(tbinStartITS.size())) { // time of TPC track exceeds the max time of ITS in the cache
      break;
    }
    int iits0 = tbinStartITS[itsROBin];
    nCheckTPCControl++;
    for (auto iits=iits0;iits<nTracksITS;iits++) {
      auto &trefITS = mITSWork[ cacheITS[iits] ];
      // compare if the ITS and TPC tracks may overlap in time
      if (trefTPC.timeMax<trefITS.timeMin) {
	// since TPC tracks are sorted in timeMax and ITS tracks are sorted in timeMin
	//all following ITS tracks also will not match
	break;
      }
      nCheckITSControl++;
      float chi2 = -1;
      int rejFlag = compareITSTPCTracks(trefITS,trefTPC,chi2);

#ifdef _ALLOW_DEBUG_TREES_      
      if ( mDBGOut && ((rejFlag==Accept && isDebugFlag(MatchTreeAccOnly)) || isDebugFlag(MatchTreeAll)) ) {
	fillITSTPCmatchTree(cacheITS[iits],cacheTPC[itpc],rejFlag, chi2);
      }
#endif
      
      if (rejFlag == RejectOnTgl) {
	// ITS tracks in each ROFrame are ordered in Tgl, hence if this check failed on Tgl check
	// (i.e. tgl_its>tgl_tpc+tolerance), tnem all other ITS tracks in this ROFrame will also have tgl too large.
	// Jump on the 1st ITS track of the next ROFrame
	int nextROF = trefITS.roFrame+1;
	if ( nextROF >= int(tbinStartITS.size()) ) { // no more ITS ROFrames in cache
	  break;
	}
	iits = tbinStartITS[nextROF]-1;  // next track to be checked -1
	continue;
      }
      if (rejFlag!=Accept) continue;
      registerMatchRecordTPC(trefITS,trefTPC,chi2);// register matching candidate
      nMatchesControl++;
    }    
  }

  // RS: this is temporary dump
  printf("TPC track checked: %d (starting from %d), total checks: %d, total matches: %d\n",
	 nCheckTPCControl,idxMinTPC,nCheckITSControl,nMatchesControl);
}

//______________________________________________
matchTPC& MatchTPCITS::getTPCMatchEntry(TrackLocTPC& tTPC)
{
  ///< return the matchTPC entry referred by tje tTPC track,
  ///< create if neaded
  if (tTPC.matchID == DummyID) { // does this TPC track already have any match? If not, create matchTPC entry 
    tTPC.matchID = mMatchesTPC.size();
    mMatchesTPC.emplace_back(tTPC.source);
    //TODO fill other info
    return mMatchesTPC.back();
  }
  return mMatchesTPC[tTPC.matchID];
}

//______________________________________________
matchITS& MatchTPCITS::getITSMatchEntry(TrackLocITS& tITS)
{
  if (tITS.matchID == DummyID) { // does this ITS track already have any match? If not, create matchITS entry 
    tITS.matchID = mMatchesITS.size();
    mMatchesITS.emplace_back(tITS.source);
    //TODO fill other info
    return mMatchesITS.back();
  }
  return mMatchesITS[tITS.matchID];
}

//______________________________________________
void MatchTPCITS::suppressMatchRecordITS(int matchITSID, int matchTPCID)
{
  ///< suppress the reference on the matchTPC with id=matchTPCID in
  ///< the list of matches recorded by for matchITS with id matchITSID
  auto& itsMatch = mMatchesITS[matchITSID];
  int recordID = itsMatch.first;     // 1st entry in mMatchRecordsITS
  while(recordID!=DummyID) { // navigate over records for given ITS track
    if (mMatchRecordsITS[recordID].matchTPCID == matchTPCID) {
      mMatchRecordsITS[recordID].matchTPCID = DummyID;
      return;
    }
    recordID = mMatchRecordsITS[recordID].nextRecID; // check next record
  }
}

//______________________________________________
bool MatchTPCITS::registerMatchRecordTPC(TrackLocITS& tITS, TrackLocTPC& tTPC, float& chi2)
{
  ///< record matching candidate, making sure that number of ITS candidates per TPC track, sorted
  ///< in matching chi2 does not exceed allowed number
  const int OverrideExisting = DummyID-100;

  auto & mtcTPC = getTPCMatchEntry(tTPC); // get matchTPC structure of this TPC track, create if none
  int nextID = mtcTPC.first;  // get 1st matchRecordTPC this matchTPC refers to
  if (nextID == DummyID) { // no matches yet, just add new record 
    registerMatchRecordITS(tITS, tTPC.matchID); // register matchTPC entry in the ITS records
    mtcTPC.first = mMatchRecordsTPC.size(); // new record will be added in the end
    mMatchRecordsTPC.emplace_back(tITS.matchID, chi2); // create new record with empty reference on next match
    return true;
  }
  
  int count=0, topID=DummyID;
  do {
    auto& nextMatchRec = mMatchRecordsTPC[nextID];
    count++;
    if (chi2<nextMatchRec.chi2) { // need to insert new record before nextMatchRec?
      if (count<mMaxMatchCandidates) {
	break; // will insert in front of nextID
      }
      else { // max number of candidates reached, will overwrite the last one
	nextMatchRec.chi2 = chi2;
	suppressMatchRecordITS(nextMatchRec.matchITSID, tTPC.matchID); // flag as disabled the overriden ITS match
	registerMatchRecordITS(tITS,tTPC.matchID); // register matchTPC entry in the ITS records
	nextMatchRec.matchITSID = tITS.matchID; // reuse the record of suppressed ITS match to store better one
	nextID = OverrideExisting; // to flag overriding existing (last) candidate
	break; // return true?
      }
    }
    topID = nextID; // check next match record
    nextID = nextMatchRec.nextRecID;
  } while (nextID!=DummyID);

  // if count == mMaxMatchCandidates, the max number of candidates was already reached, and the
  // new candidated was either discarded (if its chi2 is worst one) or has overwritten worst
  // existing candidate. Otherwise, we need to add new entry
  if (count<mMaxMatchCandidates) {
    if (topID==DummyID) { // the new match is top candidate
      mtcTPC.first = mMatchRecordsTPC.size(); // register new record as top one
    }
    else { // there are better candidates
      mMatchRecordsTPC[topID].nextRecID =  mMatchRecordsTPC.size(); // register to his parent
    }
    // nextID==-1 will mean that the while loop run over all candidates->the new one is the worst (goes to the end)
    registerMatchRecordITS(tITS,tTPC.matchID);  // register matchTPC entry in the ITS records
    mMatchRecordsTPC.emplace_back(tITS.matchID, chi2, nextID); // create new record with empty reference on next match
    return true; 
  }
  return nextID==OverrideExisting; // unless nextID was assigned OverrideExisting, new candidate was discarded
}

//______________________________________________
void MatchTPCITS::registerMatchRecordITS(TrackLocITS& tITS, int matchTPCID)
{
  ///< re
  auto & itsMatch = getITSMatchEntry(tITS); // if needed, create new entry
  int nextRecord = itsMatch.first;  // entry of 1st match record in mMatchRecordsITS
  int idnew = mMatchRecordsITS.size();
  mMatchRecordsITS.emplace_back(matchTPCID); // associate inded of matchTPD with this record
  if (nextRecord==DummyID) { // this is the 1st match for this TPC track
    itsMatch.first = idnew;
    return;
  }
  while (mMatchRecordsITS[nextRecord].nextRecID != DummyID) { // navigate to last record
    nextRecord = mMatchRecordsITS[nextRecord].nextRecID; 
  };
  mMatchRecordsITS[nextRecord].nextRecID = idnew; // register new link
  
}

//______________________________________________
int MatchTPCITS::compareITSTPCTracks(const TrackLocITS& tITS,const TrackLocTPC& tTPC, float& chi2) const
{
  ///< compare pair of ITS and TPC tracks
  auto & trackTPC = tTPC.track;
  auto & trackITS = tITS.track;
  chi2 = -1.f;
  int rejFlag = Accept;
  float diff;   // make rough check differences and their nsigmas
  
  if (std::abs(trackTPC.getAlpha()-trackITS.getAlpha())>1e-5) {
    LOG(ERROR)<<"alpha mistmatch"<<FairLogger::endl;
    trackTPC.Print();
    trackITS.Print();
  }

  // start with check on Tgl, since rjection on it will allow to profit from sorting
  diff = trackITS.getParam(o2::track::kTgl)-trackTPC.getParam(o2::track::kTgl);
  if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiffCut[o2::track::kTgl], RejectOnTgl)) ) {
    return rejFlag;
  }
  diff *= diff/(trackITS.getDiagError2(o2::track::kTgl)+trackTPC.getDiagError2(o2::track::kTgl));
  if ( (rejFlag=roughCheckDif(diff,mCrudeNSigmaCut[o2::track::kTgl], RejectOnTgl+NSigmaShift)) ) {
    return rejFlag;
  }

  diff = trackITS.getParam(o2::track::kY)-trackTPC.getParam(o2::track::kY);
  if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiffCut[o2::track::kY], RejectOnY)) ) {
    return rejFlag;
  }
  diff *= diff/(trackITS.getDiagError2(o2::track::kY)+trackTPC.getDiagError2(o2::track::kY));
  if ( (rejFlag=roughCheckDif(diff,mCrudeNSigmaCut[o2::track::kY], RejectOnY+NSigmaShift)) ) {
    return rejFlag;
  }
  
  if (mCompareTracksDZ) { // in continuous mode we usually don't use DZ
    diff = trackITS.getParam(o2::track::kZ)-trackTPC.getParam(o2::track::kZ);
    if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiffCut[o2::track::kZ],RejectOnZ)) ) {
      return rejFlag;
    }
    diff *= diff/(trackITS.getDiagError2(o2::track::kZ)+trackTPC.getDiagError2(o2::track::kZ));
    if ( (rejFlag=roughCheckDif(diff,mCrudeNSigmaCut[o2::track::kZ], RejectOnZ+NSigmaShift)) ) {
      return rejFlag;
    }    
  }

  diff = trackITS.getParam(o2::track::kSnp)-trackTPC.getParam(o2::track::kSnp);
  if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiffCut[o2::track::kSnp], RejectOnSnp)) ) {
    return rejFlag;
  }
  diff *= diff/(trackITS.getDiagError2(o2::track::kSnp)+trackTPC.getDiagError2(o2::track::kSnp));
  if ( (rejFlag=roughCheckDif(diff,mCrudeNSigmaCut[o2::track::kSnp], RejectOnSnp+NSigmaShift)) ) {
    return rejFlag;
  }
  
  diff = trackITS.getParam(o2::track::kQ2Pt)-trackTPC.getParam(o2::track::kQ2Pt);
  if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiffCut[o2::track::kQ2Pt], RejectOnQ2Pt)) ) {
    return rejFlag;
  }
  diff *= diff/(trackITS.getDiagError2(o2::track::kQ2Pt)+trackTPC.getDiagError2(o2::track::kQ2Pt));
  if ( (rejFlag=roughCheckDif(diff,mCrudeNSigmaCut[o2::track::kQ2Pt], RejectOnQ2Pt+NSigmaShift)) ) {
    return rejFlag;
  }

  // calculate mutual chi2 excluding Z in continuos mode
  chi2 = getPredictedChi2NoZ(tITS.track,tTPC.track);
  if (chi2>mCutMatchingChi2) return RejectOnChi2;

  return Accept;
}

//______________________________________________
void MatchTPCITS::printCandidatesTPC() // temporary
{
  ///< print mathing records
  
  timingOff(TimingPrintout);

  printf("\n\nPrinting all TPC -> ITS matches\n");
  
  for (auto& tTPC : mTPCWork) {
    int matchTPCID = tTPC.matchID;
    int nrec = matchTPCID==DummyID ? 0 : getNMatchRecordsTPC(matchTPCID);
    printf("*** trackTPC#%5d(%4d) : Ncand = %d\n",tTPC.source.id,tTPC.source.chunk,nrec);
    if (!nrec) {
      continue;
    }
    auto& tpcMatch = mMatchesTPC[matchTPCID];
    int count=0, recID = tpcMatch.first;
    while (recID!=DummyID) {
      const auto & recTPC = mMatchRecordsTPC[recID];
      const auto & itsMatch = mMatchesITS[recTPC.matchITSID];
      printf("  * cand %2d : ITS track %5d(%4d) Chi2: %f\n",count,itsMatch.source.id,itsMatch.source.chunk,recTPC.chi2);
      recID = recTPC.nextRecID;
      count++;
    }
  }
  timingOn(TimingPrintout);  
}

//______________________________________________
void MatchTPCITS::printCandidatesITS() // temporary
{
  ///< print mathing records
  
  timingOff(TimingPrintout);
  
  printf("\n\nPrinting all ITS -> TPC matches\n");

  for (auto & tITS : mITSWork) {
    int matchITSID = tITS.matchID;
    int nrec = matchITSID==DummyID ? 0 : getNMatchRecordsITS(matchITSID);
    printf("*** trackITS#%5d(%4d) : Ncand = %d\n",tITS.source.id,tITS.source.chunk,nrec);
    if (!nrec) {
      continue;
    }
    auto& itsMatch = mMatchesITS[matchITSID];
    int count=0, recID = itsMatch.first;
    while (recID!=DummyID) {
      const auto & recITS = mMatchRecordsITS[recID];
      const auto & tpcMatch = mMatchesTPC[recITS.matchTPCID];
      printf("  * cand %2d : TPC track %5d(%4d)\n",count,tpcMatch.source.id,tpcMatch.source.chunk);
      recID = recITS.nextRecID;
      count++;
    }
  }

  timingOn(TimingPrintout);
  
}

//______________________________________________
float MatchTPCITS::getPredictedChi2NoZ(const o2::track::TrackParCov& tr1, const o2::track::TrackParCov& tr2) const
{
  /// get chi2 between 2 tracks, neglecting Z parameter.
  /// 2 tracks must be defined at the same parameters X,alpha (check is currently commented)

  //  if (std::abs(tr1.getAlpha() - tr2.getAlpha()) > FLT_EPSILON) {
  //    LOG(ERROR) << "The reference Alpha of the tracks differ: "
  //	       << tr1.getAlpha() << " : " << tr2.getAlpha() << FairLogger::endl;
  //    return 2. * o2::track::HugeF;
  //  }
  //  if (std::abs(tr1.getX() - tr2.getX()) > FLT_EPSILON) {
  //    LOG(ERROR) << "The reference X of the tracks differ: "
  //	       << tr1.getX() << " : " << tr2.getX() << FairLogger::endl;
  //    return 2. * o2::track::HugeF;
  //  }
  MatrixDSym4 covMat;
  covMat(0, 0) = static_cast<double>(tr1.getSigmaY2())     + static_cast<double>(tr2.getSigmaY2());
  covMat(1, 0) = static_cast<double>(tr1.getSigmaSnpY())   + static_cast<double>(tr2.getSigmaSnpY());
  covMat(1, 1) = static_cast<double>(tr1.getSigmaSnp2())   + static_cast<double>(tr2.getSigmaSnp2());
  covMat(2, 0) = static_cast<double>(tr1.getSigmaTglY())   + static_cast<double>(tr2.getSigmaTglY());
  covMat(2, 1) = static_cast<double>(tr1.getSigmaTglSnp()) + static_cast<double>(tr2.getSigmaTglSnp());
  covMat(2, 2) = static_cast<double>(tr1.getSigmaTgl2())   + static_cast<double>(tr2.getSigmaTgl2());
  covMat(3, 0) = static_cast<double>(tr1.getSigma1PtY())   + static_cast<double>(tr2.getSigma1PtY());
  covMat(3, 1) = static_cast<double>(tr1.getSigma1PtSnp()) + static_cast<double>(tr2.getSigma1PtSnp());
  covMat(3, 2) = static_cast<double>(tr1.getSigma1PtTgl()) + static_cast<double>(tr2.getSigma1PtTgl());
  covMat(3, 3) = static_cast<double>(tr1.getSigma1Pt2())   + static_cast<double>(tr2.getSigma1Pt2());
  if (!covMat.Invert()) {
    LOG(ERROR) << "Cov.matrix inversion failed: " << covMat << FairLogger::endl;
    return 2. * o2::track::HugeF;
  }
  double chi2diag = 0., chi2ndiag = 0., diff[o2::track::kNParams-1] = {
    tr1.getParam(o2::track::kY)    - tr2.getParam(o2::track::kY),
    tr1.getParam(o2::track::kSnp)  - tr2.getParam(o2::track::kSnp),
    tr1.getParam(o2::track::kTgl)  - tr2.getParam(o2::track::kTgl),
    tr1.getParam(o2::track::kQ2Pt) - tr2.getParam(o2::track::kQ2Pt)
  };
  for (int i = o2::track::kNParams-1; i--;) {
    chi2diag += diff[i] * diff[i] * covMat(i, i);
    for (int j = i; j--;) {
      chi2ndiag += diff[i] * diff[j] * covMat(i, j);
    }
  }
  return chi2diag + 2. * chi2ndiag;
}

//______________________________________________
void MatchTPCITS::addTrackCloneForNeighbourSector(const TrackLocITS& src, int sector)
{
  // add clone of the src ITS track cashe, propagate it to ref.X in requested sector
  // and register its index in the sector cache. Used for ITS tracks which are so close
  // to their setctor edge that their matching should be checked also in the neighbouring sector
  
  mITSWork.push_back(src); // clone the track defined in given sector
  auto &trc = mITSWork.back().track;
  if ( trc.rotate(o2::utils::Sector2Angle(sector)) &&
       o2::Base::Propagator::Instance()->PropagateToXBxByBz(trc,mXRef)) { //TODO: use faster prop here, no 3d field, materials
    mITSSectIndexCache[sector].push_back( mITSWork.size()-1 ); // register track CLONE
    if (mMCTruthON) {
      mITSLblWork.emplace_back( mITSTrkLabels->getLabels(src.source.id)[0] );
    }
  }
  else {
    mITSWork.pop_back(); // rotation / propagation failed
  }
}

//______________________________________________
bool MatchTPCITS::propagateToRefX(o2::track::TrackParCov &trc)
{
  // propagate track to matching reference X, making sure its assigned alpha
  // is consistent with TPC sector
  bool refReached = false;
  refReached = mXRef<10.; // RS: tmp, to cover mXRef~0
  while ( o2::Base::Propagator::Instance()->PropagateToXBxByBz(trc,mXRef) ) {
    if (refReached) break; // RS: tmp
    // make sure the track is indeed within the sector defined by alpha
    if ( fabs(trc.getY()) < mXRef*tan(o2::constants::math::SectorSpanRad/2) ) {
      refReached = true;
      break; // ok, within
    }
    auto alphaNew = o2::utils::Angle2Alpha(trc.getPhiPos());
    if ( !trc.rotate( alphaNew ) != 0 ) {
      break; // failed (RS: check effect on matching tracks to neighbouring sector)
    }
  }
  return refReached && std::abs(trc.getSnp())<MaxSnp;
  
}

//______________________________________________
void MatchTPCITS::print() const
{
  ///< print all settings
  printf("\n******************** TPC-ITS matching component ********************\n");
  if (!mInitDone) {
    printf("init is not done yet\n");
    return;
  }

  printf("MC truth: %s\n",mMCTruthON ? "on":"off");
  printf("Matching reference X: %.3f\n",mXRef);
  printf("Account Z dimension: %s\n",mCompareTracksDZ ? "on":"off");
  printf("Cut on matching chi2: %.3f\n",mCutMatchingChi2);
  printf("Max number ITS candidates per TPC track: %d\n",mMaxMatchCandidates);
  printf("Crude cut on track params: ");
  for (int i=0;i<o2::track::kNParams;i++) {
    printf(" %.3e",mCrudeAbsDiffCut[i]);
  }
  printf("\n");

  printf("N Sigma cut on track params: ");
  for (int i=0;i<o2::track::kNParams;i++) {
    printf(" %6.2f",mCrudeNSigmaCut[i]);
  }
  printf("\n");

  printf("TPC-ITS time(bins) bracketing safety margin: %6.2f\n",mTimeBinTolerance);
  printf("TPC Z->time(bins) bracketing safety margin: %6.2f\n",mTPCTimeEdgeZSafeMargin);

#ifdef _ALLOW_DEBUG_TREES_
  
  printf("Output debug tree (%s) file: %s\n",mDBGFlags ? "on":"off",mDebugTreeFileName.data());
  if (getDebugFlags()) {
    printf("Debug stream flags:\n");
    if ( isDebugFlag(MatchTreeAll|MatchTreeAccOnly) ) {
      printf("* matching canditate pairs: %s\n", isDebugFlag(MatchTreeAccOnly) ? "accepted":"all");
    }
  }
#endif

  
  printf("**********************************************************************\n");
}

#ifdef _ALLOW_DEBUG_TREES_
//______________________________________________
void MatchTPCITS::setDebugFlag(UInt_t flag, bool on)
{
  ///< set debug stream flag
  if (on) {
    mDBGFlags |= flag;
  }
  else {
    mDBGFlags &= ~flag;
  }
}

//_________________________________________________________
void MatchTPCITS::fillITSTPCmatchTree(int itsID, int tpcID, int rejFlag, float chi2)
{
  ///< fill debug tree for ITS TPC tracks matching check

  timingOff(TimingFillDebugTree);
  
  auto &trackITS = mITSWork[ itsID ];
  auto &trackTPC = mTPCWork[ tpcID ];
  if (chi2<0.) { // need to recalculate
    chi2 = getPredictedChi2NoZ(trackITS.track,trackTPC.track);
  }
  o2::MCCompLabel lblITS,lblTPC;
  (*mDBGOut)<<"match"<<"chi2Match="<<chi2<<"its="<<trackITS<<"tpc="<<trackTPC;
  if (mMCTruthON) {
    lblITS = mITSLblWork[itsID];
    lblTPC = mTPCLblWork[tpcID];
    (*mDBGOut)<<"match"<<"itsLbl="<<lblITS<<"tpcLbl="<<lblTPC;
  }
  (*mDBGOut)<<"match"<<"rejFlag="<<rejFlag<<"\n";

  timingOn(TimingFillDebugTree);
}
#endif
