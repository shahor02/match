#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TChain.h>
#include <TGeoGlobalMagField.h>
#include <string>
#include <FairLogger.h>

#include "Field/MagneticField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

#include "MatchTPCITS.h"
#endif

void testMatch(std::string path = "./"
	       ,std::string outputfile="o2match_itstpc.root"
	       ,std::string inputTracksITS="o2track_its.root"
	       ,std::string inputTracksTPC="tracksFromNative.root"
	       ,std::string inputClustersITS="o2clus.root"
	       ,std::string inputGeom="O2geometry.root"
	       ,std::string inputGRP="o2sim_grp.root")
{

  o2::globaltracking::MatchTPCITS matching;

  if (path.back()!='/') {
    path += '/';
  }

  //>>>---------- attach input data --------------->>>
  TChain itsTracks("o2sim");
  itsTracks.AddFile((path+inputTracksITS).data());
  matching.setInputChainITSTracks(&itsTracks);

  TChain tpcTracks("events");
  tpcTracks.AddFile((path+inputTracksTPC).data());
  matching.setInputChainTPCTracks(&tpcTracks);

  TChain itsClusters("o2sim");
  itsClusters.AddFile((path+inputClustersITS).data());
  matching.setInputChainITSClusters(&itsClusters);

  //<<<---------- attach input data ---------------<<<
  
#ifdef _ALLOW_DEBUG_TREES_
  matching.setDebugTreeFileName(path+matching.getDebugTreeFileName());
  // dump accepted pairs only
  matching.setDebugFlag(o2::globaltracking::MatchTPCITS::MatchTreeAccOnly);
  // dump all checked pairs
  // matching.setDebugFlag(o2::globaltracking::MatchTPCITS::MatchTreeAll);
  // dump winner matches
  matching.setDebugFlag(o2::globaltracking::MatchTPCITS::WinnerMatchesTree);
#endif
  
  //-------- init geometry and field --------//
  o2::Base::GeometryManager::loadGeometry(path+inputGeom,"FAIRGeom");
  o2::Base::Propagator::initFieldFromGRP(path+inputGRP);
  
  //-------------------- settings -----------//
  matching.setITSROFrameLengthMUS(5.0f);      // ITS ROFrame duration in \mus
  matching.setCutMatchingChi2(100.);
  std::array<float,o2::track::kNParams> cutsAbs = {2.f, 2.f, 0.2f, 0.2f, 4.f};
  std::array<float,o2::track::kNParams> cutsNSig2 = {49.f,49.f,49.f,49.f,49.f};
  matching.setCrudeAbsDiffCut(cutsAbs);
  matching.setCrudeNSigma2Cut(cutsNSig2);
  
  matching.init();

  matching.run();
  
}
