#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <sstream>

#include <TStopwatch.h>
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
#include "DetectorsBase/GeometryManager.h"

#include "CommonConstants/PhysicsConstants.h"

#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TFile.h>
#include <TGeoGlobalMagField.h>
#include "DataFormatsParameters/GRPObject.h"
#include <string>
#include <array>
#include <vector>
#endif

#define _ALLOW_DEBUG_TREES_ // uncomment this to produce debug/tuning trees

TStopwatch timer;


using o2field = o2::field::MagneticField;
using o2fieldFast = o2::field::MagFieldFast;
using GRP = o2::parameters::GRPObject;

using MatrixDSym4 = ROOT::Math::SMatrix<double, 4, 4, ROOT::Math::MatRepSym<double, 4>>;
using MatrixD4 = ROOT::Math::SMatrix<double, 4, 4, ROOT::Math::MatRepStd<double, 4>>;

using namespace o2::utils;
using namespace o2::Base;
using namespace o2::TPC;
using namespace o2::ITS;

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
  ClassDefNV(MatchRecord,1);
};

const int MatchRecord::Dummy = -1;

int mCurrTPCTreeEntry=-1; ///< current TPC tree entry loaded to memory
int mCurrITSTreeEntry=-1; ///< current ITS tree entry loaded to memory

const float mXTPCInnerRef = 80.0; ///< reference radius at which TPC provides the tracks 
const float mXRef = 70.0;
const float mYMaxAtXRef = mXRef*std::tan(o2::constants::math::SectorSpanRad*0.5); ///< max Y in the sector at reference X
//const float mXRef = 70.0;

float mCrudeAbsDiff[o2::track::kNParams] = {2., 2., 0.2, 0.2, 4.}; ///<tolerance on abs. different of ITS/TPC params
float mCrudeNSigma[o2::track::kNParams] = {49.,49.,49.,49.,49.}; ///<tolerance on per-component ITS/TPC params NSigma
float mCutChi2TPCITS = 200;
float mSectEdgeMargin2 = 0; // crude check if ITS track should be matched also in neighbouring sector

///< safety margin (in TPC time bins) for ITS-TPC tracks time (in TPC time bins!) comparison 
float mITSTPCTimeBinSafeMargin = 1.; 

///< safety margin in cm when estimating TPC track tMin and tMax from assigned time0 and its
///< track Z position
float mTPCTimeEdgeZSafeMargin = 20.;

///< safety margin in TPC time bins when estimating TPC track tMin and tMax from
///< assigned time0 and its track Z position (converted from mTPCTimeEdgeZSafeMargin)
float mTPCTimeEdgeTSafeMargin = 0.;

bool mCompareTracksDZ = false; // do we use track Z difference to reject fake matches?

constexpr float TolerSortTime = 0.1; ///<tolerance for comparison of 2 tracks times
constexpr float TolerSortTgl = 1e-4; ///<tolerance for comparison of 2 tracks tgl

constexpr float Tan70 = 2.74747771e+00; // tg(70 degree): std::tan(70.*o2::constants::math::PI/180.);
constexpr float Cos70I2 = 1.+Tan70*Tan70; // 1/cos^2(70) = 1 + tan^2(70)
constexpr float MaxSnp = 0.85;
constexpr float MaxTgp = 1.61357f; // max tg corresponting to MaxSnp = MaxSnp/std::sqrt((1.-MaxSnp)*(1.+MaxSnp));
//constexpr float Max;

TChain *mChainITS=nullptr, *mChainTPC=nullptr;

float mITSROFrame = 0.; ///< ITS RO frame in \mus
float mTPCVDrift0 = 0.; ///< TPC nominal drift speed in cm/microseconds
float mITSROFrame2TPCBin = 0.; ///< conversion coeff from ITS ROFrame units to TPC time-bin
float mTPCBin2ITSROFrame = 0.; ///< conversion coeff from TPC time-bin to ITS ROFrame units
float mZ2TPCBin = 0.; ///< conversion coeff from Z to TPC time-bin
float mTPCBin2Z = 0.; ///< conversion coeff from TPC time-bin to Z
float mNTPCBinsFullDrift = 0.; ///< max time bin for full drift
float mTPCZMax = 0.;


float mTimeBinTolerance = 10.; ///<tolerance in time-bin for ITS-TPC time bracket matching

bool mMCTruthON = false;

///>>>------ these are input arrays which should not be modified by the matching code
//           since this info is provided by external device
std::vector<o2::ITS::CookedTrack> *mITSTracksArrayInp = nullptr;             ///<input tracks
std::vector<o2::TPC::TrackTPC> *mTPCTracksArrayInp = nullptr;                ///<input tracks

o2::dataformats::MCTruthContainer<o2::MCCompLabel> *mITSTrkLabels = nullptr; ///< input ITS Track MC labels
o2::dataformats::MCTruthContainer<o2::MCCompLabel> *mTPCTrkLabels = nullptr; ///< input TPC Track MC labels
/// <<<-----

std::vector<int> mMatchRecordID; ///< refs on 1st matchRecord in mMatchRecords of each TPC track
std::vector<MatchRecord> mMatchRecords; ///< match records pool
int mMaxMatchCandidates = 5; ///< max allowed matching candidates per TPC track

std::vector<TrackLocTPC> mTPCWork; ///<TPC track params prepared for matching
std::vector<TrackLocITS> mITSWork; ///<ITS track params prepared for matching
std::vector<o2::MCCompLabel> mTPCLblWork; ///<TPC track labels
std::vector<o2::MCCompLabel> mITSLblWork; ///<ITS track labels

std::array<std::vector<int>,o2::constants::math::NSectors> mTPCSectIndexCache;
std::array<std::vector<int>,o2::constants::math::NSectors> mITSSectIndexCache;

 ///<indices of 1st entries with time-bin above the value
std::array<std::vector<int>,o2::constants::math::NSectors> mTPCTimeBinStart;
///<indices of 1st entries of ITS tracks with givem ROframe
std::array<std::vector<int>,o2::constants::math::NSectors> mITSTimeBinStart; 

#ifdef _ALLOW_DEBUG_TREES_

enum DebugFlagTypes : UInt_t {
  MatchTreeAll      = 0x1<<1,    ///< produce matching candidates tree for all candidates
  MatchTreeAccOnly  = 0x1<<2     ///< fill the matching candidates tree only once the cut is passed
};

std::unique_ptr<TreeStreamRedirector> mDBGOut;
UInt_t mDBGFlags = 0;

///< check if partucular flags are set
bool isDebugFlag(UInt_t flags) { return mDBGFlags & flags; }

///< set or unset debug stream flag
void setDebugFlag(UInt_t flag, bool on=true);

///< get debug trees flags
UInt_t getDebugFlags() { return mDBGFlags; }

#endif

///< get number of matching records for TPC track
int getNMatchRecords(int tpcTrackID)
{
  int count = 0, recID = mMatchRecordID[tpcTrackID];
  while (recID!=MatchRecord::Dummy) {
    recID = mMatchRecords[recID].nextRecID;
    count++;
  }
  return count;
}

///< set ITS RO Frame, TPC time bin duration in microseconds and TPC nominal VDrift
void setTPCITSParams(float itsROFrame, float tpcTBin, float tpcVDrift, float tpcZMax)
{
  mTPCZMax = tpcZMax;
  mITSROFrame = itsROFrame;
  mITSROFrame2TPCBin = mITSROFrame/tpcTBin;
  mTPCBin2ITSROFrame = 1./mITSROFrame2TPCBin;
  mTPCVDrift0 = tpcVDrift; 
  mTPCBin2Z = tpcTBin*tpcVDrift;
  mZ2TPCBin = 1./mTPCBin2Z;
  mNTPCBinsFullDrift = mTPCZMax*mZ2TPCBin;
}

///< convert TPC time bin to ITS ROFrame units
int tpcTimeBin2ITSROFrame(float tbin) {
  int rof = tbin*mTPCBin2ITSROFrame;
  return rof<0 ? 0 : rof;
}

///< convert ITS ROFrame to TPC time bin units
float itsROFrame2TPCTimeBin(int rof) {
  return rof*mITSROFrame2TPCBin;
}

///< convert Z interval to TPC time-bins
float z2TPCBin(float z)
{
  return z*mZ2TPCBin;
}

///< convert TPC time-bins to Z interval
float tpcBin2Z(float t)
{
  return t*mTPCBin2Z;
}

///< rought check of 2 track params difference, return -1,0,1 if it is <,within or > than tolerance
int roughCheckDif(float delta, float toler, int rejFlag)
{
  return delta>toler ? rejFlag : (delta<-toler ? -rejFlag : Accept);
}

int initSimGeomAndField(std::string geomFileName="O2geometry.root",
			std::string grpFileName="o2sim_grp.root",
			std::string geomName="FAIRGeom",
			std::string grpName="GRP");
void initFieldFromGRP(const o2::parameters::GRPObject* grp);


TChain* initInputITS(std::string &inputTracks);
TChain* initInputTPC(std::string &inputTracks);
bool prepareTPCData();
bool prepareITSData();
bool loadTPCData();
bool loadITSData();
void doMatching();
void doMatching(int sec);
int compareITSTPCTracks(const TrackLocITS& tITS,const TrackLocTPC& tTPC, float& chi2);
bool registerMatchRecord(const TrackLocITS& tITS,const TrackLocTPC& tTPC, float& chi2);
float getPredictedChi2NoZ(const o2::track::TrackParCov& tr1, const o2::track::TrackParCov& tr2); //const; //RS make it const
bool propagateToRefX(o2::track::TrackParCov &trc);
void addTrackCloneForNeighbourSector(const TrackLocITS& src, int sector);
void fillITSTPCmatchTree(int itsID,int tpcID, int rejFlag, float chi2=-1.);

void printCandidates(); // temporary


const o2fieldFast* field = nullptr;



void run_match(float rate = 0. // continuous or triggered mode
	       ,std::string outputfile="o2match_itstpc.root"
	       ,std::string inputTracksITS="o2track_its.root"
	       ,std::string inputTracksTPC="tracksFromNative.root"
	       ,std::string inputGeom="O2geometry.root"
	       ,std::string inputGRP="o2sim_grp.root"
	       )
{ 
  ///<Initialize logger
  FairLogger *logger = FairLogger::GetLogger();
  logger->SetLogVerbosityLevel("LOW");
  logger->SetLogScreenLevel("INFO");

  // request debug
#ifdef _ALLOW_DEBUG_TREES_
  setDebugFlag(MatchTreeAll,true);
  //setDebugFlag(MatchTreeAccOnly,true);
#endif
  // Setup timer
  timer.Start();

#ifdef _ALLOW_DEBUG_TREES_
  // debug streamer
  if (mDBGFlags) {
    mDBGOut = std::make_unique<TreeStreamRedirector>("dbg_match.root","recreate");
  }
#endif
  
  // load geometry and field
  initSimGeomAndField(inputGeom,inputGRP);
  o2field *fieldSlow = (o2field*)TGeoGlobalMagField::Instance()->GetField();
  fieldSlow->AllowFastField(true);
  field = fieldSlow->getFastField();

  // ================= DIFFERENT INITS HERE =======================
  
  // set ITS and TPC time bin units in \mus
  const ParameterGas &gasParam = ParameterGas::defaultInstance();
  const ParameterElectronics &elParam = ParameterElectronics::defaultInstance();
  const ParameterDetector& detParam = ParameterDetector::defaultInstance();
  // TPC time-bin and ITS RO frame in \mus and max drift length
  setTPCITSParams(10., elParam.getZBinWidth(), gasParam.getVdrift(), detParam.getTPClength()); 
  mTPCTimeEdgeTSafeMargin = z2TPCBin(mTPCTimeEdgeZSafeMargin);
  // the ITS track in given sector might actually be the partner of TPC track in sector
  // above or below if it is very close to the sector edge. Here is the square of the
  // threshold for estimated distance of the track from neighboring sector edge when
  // rotated to this sector and propagated to Xref 
  mSectEdgeMargin2 = mCrudeAbsDiff[o2::track::kY]*mCrudeAbsDiff[o2::track::kY];
  
  mChainITS = initInputITS(inputTracksITS);
  mChainTPC = initInputTPC(inputTracksTPC);
  mMCTruthON = (mITSTrkLabels && mTPCTrkLabels);
  
  while(prepareTPCData()) {
    while(prepareITSData()) {    
      doMatching();
    }
    printCandidates();
  }

  timer.Stop();
  timer.Print();
  
#ifdef _ALLOW_DEBUG_TREES_
  mDBGOut.reset();
#endif

}


void doMatching()
{
  for (int  sec=o2::constants::math::NSectors;sec--;) {
    doMatching(sec);
  }
}


///*
void doMatching(int sec)
{
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
//*/
/*
void doMatching(int sec)
{
  printf("\n\nMatching sector %d\n",sec);
  
  auto &cacheITS = mITSSectIndexCache[sec]; // array of cached ITS track indices for this sector
  auto &cacheTPC = mTPCSectIndexCache[sec]; // array of cached ITS track indices for this sector
  auto &tbinStartTPC = mTPCTimeBinStart[sec]; // array of 1st TPC track with timeMax in ITS ROFrame
  int nTracksITS = cacheITS.size();
  int nTracksTPC = cacheTPC.size();
  int nTBinsCached = tbinStartTPC.size();
  if (!nTBinsCached) {
    printf("No TPC tracks in sector %d\n",sec);
    return;
  }

  for (auto iits=0;iits<nTracksITS;iits++) {
    auto &trefITS = mITSWork[ cacheITS[iits] ];

    // no point in testing TPC tracks whose timeMax < its_track_timeMin
    // we need to check TPC tracks starting from this bin
    int tbin = static_cast<int>(tpcTimeBin2ITSROFrame(trefITS.timeMin));
    if (tbin<0) tbin = 0;
    int itpc0 = tbinStartTPC[ tbin < nTBinsCached ? tbin : nTBinsCached-1 ]; //1st track index to consider
    printf("start with [%d/%d]->[%d/%d]\n",tbin,tbinStartTPC.size(), itpc0, cacheTPC.size());

    for (int itpc=itpc0;itpc<nTracksTPC; itpc++) {
      auto & trefTPC = mTPCWork[ cacheTPC[itpc] ];
      //
      if (trefITS.timeMax < trefTPC.timeMax - mNTPCBinsFullDrift) {
	// the spread between TPC tracks timeMin and timeMax cannot exceed mNTPCBinsFullDrift, hence
	// if ITS.timeMax is less than trefTPC.timeMax - mNTPCBinsFullDrift, then ITS.timeMax < timeMin
	// of all following TPC tracks (with larger timeMax)
	break;
      }
      if (trefITS.timeMax < trefTPC.timeMin) {
        continue; // no overlap, but TPC tracks sorted in timeMax, there might be still matches
      }
      //
      float chi2 = -1;
      int rejFlag = compareITSTPCTracks(trefITS,trefTPC,chi2);
      if (rejFlag) continue;
    }
    
  }
  
}
*/

//______________________________________________
bool registerMatchRecord(const TrackLocITS& tITS,const TrackLocTPC& tTPC, float& chi2)
{
  ///< record matching candidate, making sure that number of candidates per TPC track, sorted
  ///< in matching chi2 does not exceed allowed number
  const int OverrideExisting = MatchRecord::Dummy-100;
  int nextID = mMatchRecordID[tTPC.trOrigID]; // best matchRecord entry in the mMatchRecords
  if (nextID == MatchRecord::Dummy) { // no matches yet, just add new record 
    mMatchRecordID[tTPC.trOrigID] = mMatchRecords.size(); // register new record as top one
    mMatchRecords.emplace_back(tITS.trOrigID, tITS.eventID, chi2); // create new record with empty reference on next match
    return true;
  }
  int count=0, topID=MatchRecord::Dummy;
  do {
    auto& nextMatchRec = mMatchRecords[nextID];
    count++;
    if (chi2<nextMatchRec.chi2) { // need to insert new record before nextMatchRec?
      if (count<mMaxMatchCandidates) {
	break; // will insert in front of nextID
      }
      else { // max number of candidates reached, will overwrite the last one
	nextMatchRec.chi2 = chi2;
	nextMatchRec.trOrigID = tITS.trOrigID;
	nextMatchRec.eventID = tITS.eventID;
	nextID = OverrideExisting; // to flag overriding existing candidate
	break;
      }
    }
    topID = nextID; // register better parent
    nextID = nextMatchRec.nextRecID;
  } while (nextID!=MatchRecord::Dummy);

  // if count == mMaxMatchCandidates, the max number of candidates was already reached, and the
  // new candidated was either discarded (if its chi2 is worst one) or has overwritten worst
  // existing candidate. Otherwise, we need to add new entry
  if (count<mMaxMatchCandidates) {
    if (topID==MatchRecord::Dummy) { // the new match one is top candidate
      mMatchRecordID[tTPC.trOrigID] = mMatchRecords.size(); // register new record as top one
    }
    else { // there are better candidates
      mMatchRecords[topID].nextRecID =  mMatchRecords.size(); // register to his parent
    }
    // nextID==-1 will mean that the while loop run over all candidates->the new one is the worst (goes to the end)
    mMatchRecords.emplace_back(tITS.trOrigID,tITS.eventID, chi2, nextID); // create new record with empty reference on next match
    return true; 
  }
  return nextID==OverrideExisting; // unless nextID was assigned OverrideExisting, new candidate was discarded
}

//______________________________________________
int compareITSTPCTracks(const TrackLocITS& tITS,const TrackLocTPC& tTPC, float& chi2)
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
void fillITSTPCmatchTree(int itsID, int tpcID, int rejFlag, float chi2)
{
  ///< fill debug tree for ITS TPC tracks matching check
  timer.Stop();
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

  timer.Start(kFALSE);
}
#endif

TChain* initInputITS(std::string &inputTracks)
{
  // ITS input    
  TChain *ch = new TChain("o2sim");
  ch->AddFile(inputTracks.data());
  ch->SetBranchAddress("ITSTrack",&mITSTracksArrayInp);
  
  // is there MC info available ?
  if (ch->GetBranch("ITSTrackMCTruth")) {
    ch->SetBranchAddress("ITSTrackMCTruth",&mITSTrkLabels);
  }
  mCurrITSTreeEntry = -1;

  return ch;
}

TChain* initInputTPC(std::string &inputTracks)
{
  // TPC input
  TChain* ch = new TChain("events");
  ch->AddFile(inputTracks.data());
  ch->SetBranchAddress("Tracks",&mTPCTracksArrayInp);

  // is there MC info available ?
  if (ch->GetBranch("TracksMCTruth")) {
    ch->SetBranchAddress("TracksMCTruth",&mTPCTrkLabels);
  }
  mCurrTPCTreeEntry = -1;
  return ch;
}

bool prepareTPCData()
{
  // load next chunk of TPC data and prepare for matching
  mMatchRecordID.clear();
  mMatchRecords.clear();
  
  timer.Stop();
  if (!loadTPCData()) {
    timer.Start(kFALSE);
    return false;
  }
  timer.Start(kFALSE);

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
	      [](int a, int b) {	       
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
bool prepareITSData()
{
  // load next chunk of ITS data and prepare for matching
  timer.Stop();
  if (!loadITSData()) {
    timer.Start(kFALSE);
    return false;
  }
  timer.Start(kFALSE);
  
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
	      [](int a, int b) {
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

bool loadITSData()
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

bool loadTPCData()
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


int initSimGeomAndField(std::string geomFileName,std::string grpFileName,
			std::string geomName,std::string grpName)
{
  LOG(INFO)<<"Loading geometry from "<<geomFileName.data()<<FairLogger::endl;
  TFile flGeom(geomFileName.data());
  if ( flGeom.IsZombie() ) return -1;
  if ( !flGeom.Get(geomName.data()) ) return -2;

  //
  LOG(INFO)<<"Loading field from GRP of "<<grpFileName.data()<<FairLogger::endl;
  TFile flGRP(grpFileName.data());
  if ( flGRP.IsZombie() ) return -10;
  auto grp = static_cast<GRP*>(flGRP.GetObjectChecked(grpName.data(),
						      GRP::Class()));
  if (!grp) return -12;
  grp->print();
  initFieldFromGRP(grp);
  return 0;
}

void initFieldFromGRP(const GRP* grp)
{
  if ( TGeoGlobalMagField::Instance()->IsLocked() ) {
    if (TGeoGlobalMagField::Instance()->GetField()->TestBit(o2field::kOverrideGRP)) {
      printf("ExpertMode!!! GRP information will be ignored\n");
      printf("ExpertMode!!! Running with the externally locked B field\n");
      return;
    }
    else {
      printf("Destroying existing B field instance\n");
      delete TGeoGlobalMagField::Instance();
    }
  }
  auto fld = o2field::createFieldMap(grp->getL3Current(),grp->getDipoleCurrent());
  TGeoGlobalMagField::Instance()->SetField( fld );
  TGeoGlobalMagField::Instance()->Lock();
  printf("Running with the B field constructed out of GRP\n");
  printf("Access field via TGeoGlobalMagField::Instance()->Field(xyz,bxyz) or via\n");
  printf("auto o2field = static_cast<o2::field::MagneticField*>( TGeoGlobalMagField::Instance()->GetField() )\n");
  
}

//______________________________________________
float getPredictedChi2NoZ(const o2::track::TrackParCov& tr1, const o2::track::TrackParCov& tr2) //RS make it const
{
  /// get chi2 between 2 tracks, neglecting Z parameter.
  /// 2 tracks must be defined at the same parameters X,alpha (check is currently commented)

  //  if (std::abs(tr1.getAlpha() - tr2.getAlpha()) > FLT_EPSILON) {
  //    LOG(ERROR) << "The reference Alpha of the tracks differ: "
  //	       << tr1.getAlpha() << " : " << tr2.getAlpha() << FairLogger::endl;
  //    return 2. * HugeF;
  //  }
  //  if (std::abs(tr1.getX() - tr2.getX()) > FLT_EPSILON) {
  //    LOG(ERROR) << "The reference X of the tracks differ: "
  //	       << tr1.getX() << " : " << tr2.getX() << FairLogger::endl;
  //    return 2. * HugeF;
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
void addTrackCloneForNeighbourSector(const TrackLocITS& src, int sector)
{
  // add clone of the src ITS track cashe, propagate it to ref.X in requested sector
  // and register its index in the sector cache. Used for ITS tracks which are so close
  // to their setctor edge that their matching should be checked also in the neighbouring sector
  
  mITSWork.push_back(src); // clone the track defined in given sector
  auto &trc = mITSWork.back().track;
  if ( trc.rotate(o2::utils::Sector2Angle(sector)) &&
       Propagator::Instance()->PropagateToXBxByBz(trc,mXRef)) { //TODO: use faster prop here, no 3d field, materials
    mITSSectIndexCache[sector].push_back( mITSWork.size()-1 ); // register track CLONE
    if (mMCTruthON) {
      mITSLblWork.emplace_back( mITSTrkLabels->getLabels(src.trOrigID)[0] );
    }
  }
  else {
    mITSWork.pop_back(); // rotation / propagation failed
  }
}

//______________________________________________
bool propagateToRefX(o2::track::TrackParCov &trc)
{
  // propagate track to matching reference X, making sure its assigned alpha
  // is consistent with TPC sector
  bool refReached = false;
  refReached = mXRef<10.; // RS: tmp, to cover mXRef~0
  while ( Propagator::Instance()->PropagateToXBxByBz(trc,mXRef) ) {
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

void printCandidates() // temporary
{
  for (int itr=0;itr<mTPCTracksArrayInp->size();itr++) {
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
