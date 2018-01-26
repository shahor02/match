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
#include "DetectorsBase/Utils.h"
#include "DetectorsBase/Constants.h"
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

TStopwatch timer;


using o2field = o2::field::MagneticField;
using o2fieldFast = o2::field::MagFieldFast;
using GRP = o2::parameters::GRPObject;

using MatrixDSym4 = ROOT::Math::SMatrix<double, 4, 4, ROOT::Math::MatRepSym<double, 4>>;
using MatrixD4 = ROOT::Math::SMatrix<double, 4, 4, ROOT::Math::MatRepStd<double, 4>>;

using namespace o2::utils;
using namespace o2::Base;
using namespace o2::Base::Track;
using namespace o2::TPC;
using namespace o2::ITS;


struct TrackLocTPC {
  TrackParCov track;
  int trOrigID = -1; 
  float timeMin = 0.f;
  float timeMax = 0.f;
  TrackLocTPC(const TrackParCov& src, int id) : track(src),trOrigID(id) {}
  ClassDefNV(TrackLocTPC,1);
};

struct TrackLocITS {
  TrackParCov track;
  int trOrigID = -1;
  int roFrame = -1;
  float timeMin = 0.f;
  float timeMax = 0.f;
  TrackLocITS(const TrackParCov& src, int id) : track(src),trOrigID(id) {}
  ClassDefNV(TrackLocITS,1);
};

int mCurrTPCTreeEntry=-1; ///< current TPC tree entry loaded to memory
int mCurrITSTreeEntry=-1; ///< current ITS tree entry loaded to memory

const float xTPCInnerRef = 80.0; ///< reference radius at which TPC provides the tracks 
const float xRef = 70.0;
const float yMaxAtXRef = xRef*std::tan(o2::Base::Constants::kSectorSpanRad*0.5); ///< max Y in the sector at reference X
//const float xRef = 70.0;
const float rMaxITS = 39.3; ///<RS??
const float zMaxITS = 85.; ///<with margin RS??

float mCrudeAbsDiff[Track::kNParams] = {2., 2., 0.2, 0.2, 4.}; ///<tolerance on abs. different of ITS/TPC params
float mCrudeNSigma[Track::kNParams] = {49.,49.,49.,49.,49.}; ///<tolerance on per-component ITS/TPC params NSigma
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

constexpr float Tan70 = 2.74747771e+00; // tg(70 degree): std::tan(70.*o2::Base::Constants::kPI/180.);
constexpr float Cos70I2 = 1.+Tan70*Tan70; // 1/cos^2(70) = 1 + tan^2(70)
constexpr float MaxSnp = 0.85;
constexpr float MaxTgp = 1.61357f; // max tg corresponting to MaxSnp = MaxSnp/std::sqrt((1.-MaxSnp)*(1.+MaxSnp));
//constexpr float Max;

TChain *chainITS=nullptr, *chainTPC=nullptr;

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
std::unique_ptr<TreeStreamRedirector> mDBGOut;
UInt_t mDBGFlags = 0;

enum DebugFlagTypes : UInt_t {
  MatchTreeAll      = 0x1<<1,    ///< produce matching candidates tree for all candidates
  MatchTreeAccOnly  = 0x1<<2     ///< fill the matching candidates tree only once the cut is passed
};

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

///< check if partucular flags are set
bool isDebugFlag(UInt_t flags) { return mDBGFlags & flags; }

///< set or unset debug stream flag
void setDebugFlag(UInt_t flag, bool on=true);

///< get debug trees flags
UInt_t getDebugFlags() { return mDBGFlags; }

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
float getPredictedChi2NoZ(const TrackParCov& tr1, const TrackParCov& tr2); //const; //RS make it const
bool propagateToRefX(o2::Base::Track::TrackParCov &trc);
void addTrackCloneForNeighbourSector(const TrackLocITS& src, int sector);
void fillITSTPCmatchTree(const TrackLocITS& tITS,const TrackLocTPC& tTPC, int rejFlag, float chi2=-1.);

const o2fieldFast* field = nullptr;



///>>>------ these are input arrays which should not be modified by the matching code, but by IO only
std::vector<o2::ITS::CookedTrack> *mITSTracksArrayInp = nullptr;
std::vector<o2::TPC::TrackTPC> *mTPCTracksArrayInp = nullptr;  ///<input tracks

o2::dataformats::MCTruthContainer<o2::MCCompLabel> *mITSTrkLabels = nullptr; ///< ITS Track MC labels
o2::dataformats::MCTruthContainer<o2::MCCompLabel> *mTPCTrkLabels = nullptr; ///< TPC Track MC labels
/// <<<-----


using PointF = Point3D<float>;

std::vector<TrackLocTPC> mTPCWork; ///<TPC track params prepared for matching
std::vector<TrackLocITS> mITSWork; ///<ITS track params prepared for matching
std::array<std::vector<int>,o2::Base::Constants::kNSectors> mTPCSectIndexCache;
std::array<std::vector<int>,o2::Base::Constants::kNSectors> mITSSectIndexCache;

 ///<indices of 1st entries with time-bin above the value
std::array<std::vector<int>,o2::Base::Constants::kNSectors> mTPCTimeBinStart;
///<indices of 1st entries of ITS tracks with givem ROframe
std::array<std::vector<int>,o2::Base::Constants::kNSectors> mITSTimeBinStart; 

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
  setDebugFlag(MatchTreeAll,true);
  //setDebugFlag(MatchTreeAccOnly,true);
  
  // Setup timer
  timer.Start();

  // debug streamer
  if (mDBGFlags) {
    mDBGOut = std::make_unique<TreeStreamRedirector>("dbg_match.root","recreate");
  }
  
  // load geometry and field
  initSimGeomAndField(inputGeom,inputGRP);
  o2field *fieldSlow = (o2field*)TGeoGlobalMagField::Instance()->GetField();
  fieldSlow->AllowFastField(true);
  field = fieldSlow->getFastField();

  // ================= DIFFEREN INITS HERE =======================
  
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
  mSectEdgeMargin2 = mCrudeAbsDiff[Track::kY]*mCrudeAbsDiff[Track::kY];
  
  chainITS = initInputITS(inputTracksITS);
  chainTPC = initInputTPC(inputTracksTPC);
  mMCTruthON = (mITSTrkLabels && mTPCTrkLabels);
  
  while(prepareTPCData()) {
    while(prepareITSData()) {    
      doMatching();
    }
  }

  timer.Stop();
  timer.Print();
  
  mDBGOut.reset();
}


void doMatching()
{
  for (int  sec=o2::Base::Constants::kNSectors;sec--;) {
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
    if (itsROBin>=tbinStartITS.size()) { // time of TPC track exceeds the max time of ITS in the cache
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
      
      if ( mDBGOut && ((rejFlag==Accept && isDebugFlag(MatchTreeAccOnly)) || isDebugFlag(MatchTreeAll)) ) {
	fillITSTPCmatchTree(trefITS,trefTPC,rejFlag, chi2);
      }
      
      if (rejFlag == RejectOnTgl) {
	// ITS tracks in each ROFrame are ordered in Tgl, hence if this check failed on Tgl check
	// (i.e. tgl_its>tgl_tpc+tolerance), tnem all other ITS tracks in this ROFrame will also have tgl too large.
	// Jump on the 1st ITS track of the next ROFrame
	int nextROF = trefITS.roFrame+1;
	if ( nextROF >=  tbinStartITS.size()) { // no more ITS ROFrames in cache
	  break;
	}
	printf("JUMP from %d to %d\n",iits, tbinStartITS[nextROF]-1);
	iits = tbinStartITS[nextROF]-1;  // next track to be checked 
	continue;
      }
      if (rejFlag!=Accept) continue;
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
  diff = trackITS.getParam(Track::kTgl)-trackTPC.getParam(Track::kTgl);
  if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiff[Track::kTgl], RejectOnTgl)) ) {
    return rejFlag;
  }
  diff *= diff/(trackITS.getDiagError2(Track::kTgl)+trackTPC.getDiagError2(Track::kTgl));
  if ( (rejFlag=roughCheckDif(diff,mCrudeNSigma[Track::kTgl], RejectOnTgl+NSigmaShift)) ) {
    return rejFlag;
  }

  diff = trackITS.getParam(Track::kY)-trackTPC.getParam(Track::kY);
  if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiff[Track::kY], RejectOnY)) ) {
    return rejFlag;
  }
  diff *= diff/(trackITS.getDiagError2(Track::kY)+trackTPC.getDiagError2(Track::kY));
  if ( (rejFlag=roughCheckDif(diff,mCrudeNSigma[Track::kY], RejectOnY+NSigmaShift)) ) {
    return rejFlag;
  }
  
  if (mCompareTracksDZ) { // in continuous mode we usually don't use DZ
    diff = trackITS.getParam(Track::kZ)-trackTPC.getParam(Track::kZ);
    if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiff[Track::kZ],RejectOnZ)) ) {
      return rejFlag;
    }
    diff *= diff/(trackITS.getDiagError2(Track::kZ)+trackTPC.getDiagError2(Track::kZ));
    if ( (rejFlag=roughCheckDif(diff,mCrudeNSigma[Track::kZ], RejectOnZ+NSigmaShift)) ) {
      return rejFlag;
    }    
  }

  diff = trackITS.getParam(Track::kSnp)-trackTPC.getParam(Track::kSnp);
  if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiff[Track::kSnp], RejectOnSnp)) ) {
    return rejFlag;
  }
  diff *= diff/(trackITS.getDiagError2(Track::kSnp)+trackTPC.getDiagError2(Track::kSnp));
  if ( (rejFlag=roughCheckDif(diff,mCrudeNSigma[Track::kSnp], RejectOnSnp+NSigmaShift)) ) {
    return rejFlag;
  }
  
  diff = trackITS.getParam(Track::kQ2Pt)-trackTPC.getParam(Track::kQ2Pt);
  if ( (rejFlag=roughCheckDif(diff,mCrudeAbsDiff[Track::kQ2Pt], RejectOnQ2Pt)) ) {
    return rejFlag;
  }
  diff *= diff/(trackITS.getDiagError2(Track::kQ2Pt)+trackTPC.getDiagError2(Track::kQ2Pt));
  if ( (rejFlag=roughCheckDif(diff,mCrudeNSigma[Track::kQ2Pt], RejectOnQ2Pt+NSigmaShift)) ) {
    return rejFlag;
  }

  // calculate mutual chi2 excluding Z in continuos mode
  chi2 = getPredictedChi2NoZ(tITS.track,tTPC.track);
  if (chi2>mCutChi2TPCITS) return RejectOnChi2;

  return Accept;
}

//_________________________________________________________
void fillITSTPCmatchTree(const TrackLocITS& tITS,const TrackLocTPC& tTPC, int rejFlag, float chi2)
{
  ///< fill debug tree for ITS TPC matching check
  timer.Stop();
  
  o2::MCCompLabel lblITS,lblTPC;
  if (mMCTruthON) {
    lblITS = mITSTrkLabels->getLabels(tITS.trOrigID)[0];
    lblTPC = mTPCTrkLabels->getLabels(tTPC.trOrigID)[0];      
  }
  if (chi2<0.) { // need to recalculate
    chi2 = getPredictedChi2NoZ(tITS.track,tTPC.track);
  }
  auto trackITSnc = tITS; // we need non-const versions for streamer
  auto trackTPCnc = tTPC;
  const auto& trTPCorig = (*mTPCTracksArrayInp)[tTPC.trOrigID];
  const auto& trITSorig = (*mITSTracksArrayInp)[tITS.trOrigID];
  auto chi2ITS = trITSorig.getChi2(), chi2TPC = trTPCorig.getChi2();
  auto nclITS = trITSorig.getNumberOfClusters(), nclTPC=trTPCorig.getNClusterReferences();
  (*mDBGOut)<<"match"<<"chi2Match="<<chi2<<"its="<<trackITSnc.track<<"itsROF="<<trackITSnc.roFrame
	    <<"itsTMin="<<trackITSnc.timeMin<<"itsTMax="<<trackITSnc.timeMax<<"chi2ITS="<<chi2ITS<<"nclITS="<<nclITS
	    <<"tpc="<<trackTPCnc.track<<"tpcTMin="<<trackTPCnc.timeMin<<"tpcTMax="<<trackTPCnc.timeMax
	    <<"chi2TPC="<<chi2TPC<<"nclTPC="<<nclTPC;
  if (mMCTruthON) {
    (*mDBGOut)<<"match"<<"itsLbl="<<lblITS<<"tpcLbl="<<lblTPC;
  }
  (*mDBGOut)<<"match"<<"rejFlag="<<rejFlag<<"\n";

  timer.Start(kFALSE);

}


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
  timer.Stop();
  if (!loadTPCData()) {
    timer.Start(kFALSE);
    return false;
  }
  timer.Start(kFALSE);
  
  // copy the track params, propagate to reference X and build sector tables
  int ntr = mTPCTracksArrayInp->size();
  mTPCWork.clear();
  mTPCWork.reserve(ntr);
  for (int sec=o2::Base::Constants::kNSectors;sec--;) {
    mTPCSectIndexCache[sec].clear();
    mTPCTimeBinStart[sec].clear();
  }
  
  for (int it=0;it<ntr;it++) {
    
    o2::TPC::TrackTPC& trcOrig = (*mTPCTracksArrayInp)[it];

    // make sure the track was propagated to inner TPC radius at the ref. radius
    if (trcOrig.getX()>xTPCInnerRef+0.1) continue; // failed propagation to inner TPC radius, cannot be matched

    mTPCWork.emplace_back(static_cast<TrackParCov&>(trcOrig),it); // working copy of track param
    auto & trc = mTPCWork.back();
    
    // propagate to matching Xref
    if (!propagateToRefX(trc.track)) {
      mTPCWork.pop_back(); // discard track whose propagation to xRef failed
      continue;
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
    mTPCSectIndexCache[Utils::Angle2Sector( trc.track.getAlpha() )].push_back( mTPCWork.size()-1 ); 
  }

  // sort tracks in each sector according to their timeMax
  for (int sec=o2::Base::Constants::kNSectors; sec--;) {
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
  for (int sec=o2::Base::Constants::kNSectors;sec--;) {
    mITSSectIndexCache[sec].clear();
  }

  for (int it=0;it<ntr;it++) {
    auto& trcOrig = (*mITSTracksArrayInp)[it];

    if (trcOrig.getParamOut().getX()<1.) {
      continue; // backward refit failed
    }
    mITSWork.emplace_back(static_cast<TrackParCov&>(trcOrig.getParamOut()),it); // working copy of outer track param
    auto & trc = mITSWork.back();

    // TODO: why I did this?
    if ( !trc.track.rotate( Utils::Angle2Alpha(trc.track.getPhiPos()) ) ) {
      mITSWork.pop_back(); // discard failed track
      continue;
    } 
    // make sure the track is at the ref. radius
    if (!propagateToRefX(trc.track)) {
      mITSWork.pop_back(); // discard failed track
      continue; // add to cache only those ITS tracks which reached ref.X and have reasonable snp
    }
    trc.timeMin = itsROFrame2TPCTimeBin(trcOrig.getROFrame());
    trc.timeMax = trc.timeMin + mITSROFrame2TPCBin;
    trc.roFrame = trcOrig.getROFrame();

    // cache work track index
    int sector = Utils::Angle2Sector( trc.track.getAlpha() );
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
    float dy2Up = (yMaxAtXRef-trc.track.getY())/(tgp + Tan70);
    if ( (dy2Up*dy2Up*Cos70I2)<mSectEdgeMargin2) { // need to check this track for matching in sector up
      addTrackCloneForNeighbourSector(trc, sector<(o2::Base::Constants::kNSectors-1) ? sector+1 : 0);
    }
    // sector down
    float dy2Dn = (yMaxAtXRef+trc.track.getY())/(tgp - Tan70);
    if ( (dy2Dn*dy2Dn*Cos70I2)<mSectEdgeMargin2) { // need to check this track for matching in sector down
      addTrackCloneForNeighbourSector(trc, sector>1 ? sector-1 : o2::Base::Constants::kNSectors-1); 
    }
  }
  
  // sort tracks in each sector according to their time, then tgl
  for (int sec=o2::Base::Constants::kNSectors; sec--;) {
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
  while(++mCurrITSTreeEntry < chainITS->GetEntries()) {
    chainITS->GetEntry(mCurrITSTreeEntry);
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
  while(++mCurrTPCTreeEntry < chainTPC->GetEntries()) {
    chainTPC->GetEntry(mCurrTPCTreeEntry);
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
float getPredictedChi2NoZ(const TrackParCov& tr1, const TrackParCov& tr2) //RS make it const
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
    return 2. * HugeF;
  }
  double chi2diag = 0., chi2ndiag = 0., diff[Track::kNParams-1] = {
    tr1.getParam(Track::kY)    - tr2.getParam(Track::kY),
    tr1.getParam(Track::kSnp)  - tr2.getParam(Track::kSnp),
    tr1.getParam(Track::kTgl)  - tr2.getParam(Track::kTgl),
    tr1.getParam(Track::kQ2Pt) - tr2.getParam(Track::kQ2Pt)
  };
  for (int i = Track::kNParams-1; i--;) {
    chi2diag += diff[i] * diff[i] * covMat(i, i);
    for (int j = i; j--;) {
      chi2ndiag += diff[i] * diff[j] * covMat(i, j);
    }
  }
  return chi2diag + 2. * chi2ndiag;
}

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

//______________________________________________
void addTrackCloneForNeighbourSector(const TrackLocITS& src, int sector)
{
  // add clone of the src ITS track cashe, propagate it to ref.X in requested sector
  // and register its index in the sector cache. Used for ITS tracks which are so close
  // to their setctor edge that their matching should be checked also in the neighbouring sector
  
  mITSWork.push_back(src); // clone the track defined in given sector
  auto &trc = mITSWork.back().track;
  if ( trc.rotate(Utils::Sector2Angle(sector)) &&
       Propagator::Instance()->PropagateToXBxByBz(trc,xRef)) { //TODO: use faster prop here, no 3d field, materials
    mITSSectIndexCache[sector].push_back( mITSWork.size()-1 ); // register track CLONE
  }
  else {
    mITSWork.pop_back(); // rotation / propagation failed
  }
}

//______________________________________________
bool propagateToRefX(o2::Base::Track::TrackParCov &trc)
{
  // propagate track to matching reference X, making sure its assigned alpha
  // is consistent with TPC sector
  bool refReached = false;
  refReached = xRef<10.; // RS: tmp, to cover xRef~0
  while ( Propagator::Instance()->PropagateToXBxByBz(trc,xRef) ) {
    if (refReached) break; // RS: tmp
    // make sure the track is indeed within the sector defined by alpha
    if ( fabs(trc.getY()) < xRef*tan(o2::Base::Constants::kSectorSpanRad/2) ) {
      refReached = true;
      break; // ok, within
    }
    auto alphaNew = Utils::Angle2Alpha(trc.getPhiPos());
    if ( !trc.rotate( alphaNew ) != 0 ) {
      break; // failed (RS: check effect on matching tracks to neighbouring sector)
    }
  }
  return refReached && std::abs(trc.getSnp())<MaxSnp;
  
}
