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
#include "CommonUtils/TreeStreamRedirector.h"

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

const int MatchRecord::Dummy = -1;

//______________________________________________
void MatchTPCITS::run()
{
  ///< perform matching for provided input
  if (!mInitDone) {
    LOG(FATAL)<<"init() was not done yet"<<FairLogger::endl;
  }
  mTimer.Start();
  
  while(prepareTPCData()) {
    while(prepareITSData()) {
      for (int  sec=o2::constants::math::NSectors;sec--;) {
	doMatching(sec);
      }
    }
    printCandidates();
  }

  mTimer.Stop();
  mTimer.Print();

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
  mSectEdgeMargin2 = mCrudeAbsDiff[o2::track::kY]*mCrudeAbsDiff[o2::track::kY]; ///< precalculated ^2

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
    mDBGOut = std::make_unique<TreeStreamRedirector>("dbg_match.root","recreate");
  }
#endif

  mInitDone = true;
  
}


//______________________________________________
int MatchTPCITS::getNMatchRecords(int tpcTrackID) const
{
  ///< get number of matching records for TPC track
  int count = 0, recID = mMatchRecordID[tpcTrackID];
  while (recID!=MatchRecord::Dummy) {
    recID = mMatchRecords[recID].nextRecID;
    count++;
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
  mMatchRecordID.clear();
  mMatchRecords.clear();
  
  mTimer.Stop();
  if (!loadTPCData()) {
    mTimer.Start(kFALSE);
    return false;
  }
  mTimer.Start(kFALSE);

  // prepare empty matching records for each TPC track
  mMatchRecordID.resize(mTPCTracksArrayInp->size(),MatchRecord::Dummy);
  // number of records might be actually more than N tracks!
  mMatchRecords.reserve(mTPCTracksArrayInp->size()); 

  // copy the track params, propagate to reference X and build sector tables
  int ntr = mTPCTracksArrayInp->size();
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

    mTPCWork.emplace_back(static_cast<o2::track::TrackParCov&>(trcOrig),it); // working copy of track param
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
  mTimer.Stop();
  if (!loadITSData()) {
    mTimer.Start(kFALSE);
    return false;
  }
  mTimer.Start(kFALSE);
  
  int ntr = mITSTracksArrayInp->size();
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
  while(++mCurrITSTreeEntry < mChainITS->GetEntries()) {
    mChainITS->GetEntry(mCurrITSTreeEntry);
    LOG(INFO)<<"Starting ITS entry "<<mCurrITSTreeEntry<<" -> "
	     <<mITSTracksArrayInp->size()<<" tracks"<<FairLogger::endl;
    if (!mITSTracksArrayInp->size()) {
      continue;
    }
    return true;
  }
  --mCurrITSTreeEntry;
  return false;
}

//_____________________________________________________
bool MatchTPCITS::loadTPCData()
{
  ///< load next chunk of TPC data
  while(++mCurrTPCTreeEntry < mChainTPC->GetEntries()) {
    mChainTPC->GetEntry(mCurrTPCTreeEntry);
    LOG(INFO)<<"Starting TPC entry "<<mCurrTPCTreeEntry<<" -> "
	     <<mTPCTracksArrayInp->size()<<" tracks"<<FairLogger::endl;
    if (mTPCTracksArrayInp->size()<1) {
      continue;
    }
    return true;
  }
  --mCurrTPCTreeEntry;
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
    for (auto iits=iits0;iits<nTracksITS;iits++) {
      auto &trefITS = mITSWork[ cacheITS[iits] ];
      // compare if the ITS and TPC tracks may overlap in time
      if (trefTPC.timeMax<trefITS.timeMin) {
	// since TPC tracks are sorted in timeMax and ITS tracks are sorted in timeMin
	//all following ITS tracks also will not match
	break;
      }
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
      registerMatchRecord(trefITS,trefTPC,chi2);// register matching candidate
    }    
  }  
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
  if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiff[o2::track::kTgl], RejectOnTgl)) ) {
    return rejFlag;
  }
  diff *= diff/(trackITS.getDiagError2(o2::track::kTgl)+trackTPC.getDiagError2(o2::track::kTgl));
  if ( (rejFlag=roughCheckDif(diff,mCrudeNSigma[o2::track::kTgl], RejectOnTgl+NSigmaShift)) ) {
    return rejFlag;
  }

  diff = trackITS.getParam(o2::track::kY)-trackTPC.getParam(o2::track::kY);
  if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiff[o2::track::kY], RejectOnY)) ) {
    return rejFlag;
  }
  diff *= diff/(trackITS.getDiagError2(o2::track::kY)+trackTPC.getDiagError2(o2::track::kY));
  if ( (rejFlag=roughCheckDif(diff,mCrudeNSigma[o2::track::kY], RejectOnY+NSigmaShift)) ) {
    return rejFlag;
  }
  
  if (mCompareTracksDZ) { // in continuous mode we usually don't use DZ
    diff = trackITS.getParam(o2::track::kZ)-trackTPC.getParam(o2::track::kZ);
    if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiff[o2::track::kZ],RejectOnZ)) ) {
      return rejFlag;
    }
    diff *= diff/(trackITS.getDiagError2(o2::track::kZ)+trackTPC.getDiagError2(o2::track::kZ));
    if ( (rejFlag=roughCheckDif(diff,mCrudeNSigma[o2::track::kZ], RejectOnZ+NSigmaShift)) ) {
      return rejFlag;
    }    
  }

  diff = trackITS.getParam(o2::track::kSnp)-trackTPC.getParam(o2::track::kSnp);
  if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiff[o2::track::kSnp], RejectOnSnp)) ) {
    return rejFlag;
  }
  diff *= diff/(trackITS.getDiagError2(o2::track::kSnp)+trackTPC.getDiagError2(o2::track::kSnp));
  if ( (rejFlag=roughCheckDif(diff,mCrudeNSigma[o2::track::kSnp], RejectOnSnp+NSigmaShift)) ) {
    return rejFlag;
  }
  
  diff = trackITS.getParam(o2::track::kQ2Pt)-trackTPC.getParam(o2::track::kQ2Pt);
  if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiff[o2::track::kQ2Pt], RejectOnQ2Pt)) ) {
    return rejFlag;
  }
  diff *= diff/(trackITS.getDiagError2(o2::track::kQ2Pt)+trackTPC.getDiagError2(o2::track::kQ2Pt));
  if ( (rejFlag=roughCheckDif(diff,mCrudeNSigma[o2::track::kQ2Pt], RejectOnQ2Pt+NSigmaShift)) ) {
    return rejFlag;
  }

  // calculate mutual chi2 excluding Z in continuos mode
  chi2 = getPredictedChi2NoZ(tITS.track,tTPC.track);
  if (chi2>mCutChi2TPCITS) return RejectOnChi2;

  return Accept;
}

//______________________________________________
void MatchTPCITS::printCandidates() const // temporary
{
  for (int itr=0;itr<int(mTPCTracksArrayInp->size());itr++) {
    int nrec = getNMatchRecords(itr);
    printf("*** trackTPC#%5d : Ncand = %d\n",itr,nrec);
    int count=0, recID = mMatchRecordID[itr];
    while (recID!=MatchRecord::Dummy) {
      auto& rec = mMatchRecords[recID];
      printf("  * cand %2d : ITS track %5d(%4d) Chi2: %f\n",count,rec.trOrigID,rec.eventID,rec.chi2);
      recID = mMatchRecords[recID].nextRecID;
      count++;
    }
  }
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


#ifdef _ALLOW_DEBUG_TREES_
//______________________________________________
void setDebugFlag(UInt_t flag, bool on)
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
  mTimer.Stop();
  // Note: we cannot use (auto &) here since need non-const versions of the tracks
  auto trackITSnc = mITSWork[ itsID ]; // ITS track copy
  auto trackTPCnc = mTPCWork[ tpcID ]; // TPC track copy
  if (chi2<0.) { // need to recalculate
    chi2 = getPredictedChi2NoZ(trackITSnc.track,trackTPCnc.track);
  }
  o2::MCCompLabel lblITS,lblTPC;
  (*mDBGOut)<<"match"<<"chi2Match="<<chi2<<"its="<<trackITSnc<<"tpc="<<trackTPCnc;
  if (mMCTruthON) {
    lblITS = mITSLblWork[itsID];
    lblTPC = mTPCLblWork[tpcID];
    (*mDBGOut)<<"match"<<"itsLbl="<<lblITS<<"tpcLbl="<<lblTPC;
  }
  (*mDBGOut)<<"match"<<"rejFlag="<<rejFlag<<"\n";

  mTimer.Start(kFALSE);
}
#endif
