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

void testMatch(std::string outputfile="o2match_itstpc.root"
	       ,std::string inputTracksITS="o2track_its.root"
	       ,std::string inputTracksTPC="tracksFromNative.root"
	       ,std::string inputGeom="O2geometry.root"
	       ,std::string inputGRP="o2sim_grp.root")
{

  o2::globaltracking::MatchTPCITS matching;
  
  TChain its("o2sim"),tpc("events");
  its.AddFile(inputTracksITS.data());
  tpc.AddFile(inputTracksTPC.data());

  matching.setInputChainITS(&its);
  matching.setInputChainTPC(&tpc);

  //-------- init geometry and field --------//
  o2::Base::GeometryManager::loadGeometry(inputGeom,"FAIRGeom");
  o2::Base::Propagator::initFieldFromGRP(inputGRP);
  
  //-------------------- settings -----------//
  matching.setITSROFrame(10.0f); // ITS ROFrame duration


  // require debug tree output
  matching.setDebugFlag(o2::globaltracking::MatchTPCITS::MatchTreeAll,true);
  
  matching.init();

  matching.run();
  
}
