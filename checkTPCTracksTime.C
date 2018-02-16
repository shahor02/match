#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGeoGlobalMagField.h>
#include <string>
#include <FairLogger.h>
#include "FairEventHeader.h"

#include "Field/MagneticField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "TPCBase/ParameterElectronics.h"
#include "TPCBase/ParameterDetector.h"
#include "TPCBase/ParameterGas.h"

#endif

TH2F *histoA = nullptr,*histoC = nullptr;
TCanvas* canv = nullptr;

void checkTPCTracksTime(std::string path = "./"
			,std::string inputTracksTPC="tracksFromNative.root"
			,std::string inputDigits="o2dig.root"
			,std::string inputGeom="O2geometry.root"
			,std::string inputGRP="o2sim_grp.root")
{

  if (path.back()!='/') {
    path += '/';
  }

  histoA = new TH2F("Aside","A: #Delta T;z_{tr,x=0}, cm;T_{true}-T_{est}, #mus",100,-5,5,100,-2,2);
  histoC = new TH2F("Cside","C: #Delta T;z_{tr,x=0}, cm;T_{true}-T_{est}, #mus",100,-5,5,100,-2,2);
  
  //-------- init geometry and field --------//
  o2::Base::GeometryManager::loadGeometry(path+inputGeom,"FAIRGeom");
  o2::Base::Propagator::initFieldFromGRP(path+inputGRP);
  o2::field::MagneticField* slowField = 
    static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
  double origin[3]={0.};
  float bz = slowField->getBz(origin);

  //>>>---------- attach input data --------------->>>
  std::vector<o2::TPC::TrackTPC> *mTPCTracksArrayInp = nullptr;
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> *mTPCTrkLabels = nullptr;
  FairEventHeader* digiHeader = nullptr;
    
  TChain tpcTracks("events");
  tpcTracks.AddFile((path+inputTracksTPC).data());
  tpcTracks.SetBranchAddress("Tracks",&mTPCTracksArrayInp);
  tpcTracks.SetBranchAddress("TracksMCTruth",&mTPCTrkLabels);

  // digitization event record
  TChain digits("o2sim");
  digits.AddFile((path+inputDigits).data());
  TBranch* brHead = digits.GetBranch("EventHeader.");
  brHead->SetAddress(&digiHeader);
  //<<<---------- attach input data ---------------<<<

  const auto & gasParam = o2::TPC::ParameterGas::defaultInstance();
  const auto & elParam = o2::TPC::ParameterElectronics::defaultInstance();
  const auto & detParam = o2::TPC::ParameterDetector::defaultInstance();
  float mTPCTBinMUS = elParam.getZBinWidth();
  float mTPCVDrift0 = gasParam.getVdrift();
  float mTPCZMax = detParam.getTPClength();
  float mTPCBin2Z = mTPCTBinMUS*mTPCVDrift0;

  
  int lastHeadIDRed = -1;
  float evTime = 0;
  for (int itpc=0;itpc<tpcTracks.GetEntries();itpc++) {
    tpcTracks.GetEntry(itpc);
    for (int itr=0;itr<(*mTPCTracksArrayInp).size();itr++) {
      const auto& track = (*mTPCTracksArrayInp)[itr];
      const auto& lbl = mTPCTrkLabels->getLabels(itr)[0];
      if (lbl.getTrackID()<0 || track.getPt()<1.) continue; // ignore fake or low pt tracks
      // test if the track was indeed constrained to Z=0
      float z = track.getZAt(0,bz);
      if (std::abs(z)>5) continue; // failed propagation
      
      auto evID = lbl.getEventID();
      if (lastHeadIDRed!=evID) {  // get event time in \mus
	brHead->GetEntry( lastHeadIDRed = evID );
	evTime = digiHeader->GetEventTime()/1e3;
	printf("Loading header %d, time = %f\n",evID,evTime);
      }
      // track time in \mus under assumption it comes from Z=0
      float trackTime = (track.getTime0() - mTPCZMax/mTPCBin2Z)*mTPCTBinMUS;
      // fill histos for each side separately
      if (track.hasCSideClustersOnly()) {
	histoC->Fill(z, evTime - trackTime);
      }
      else {
	histoA->Fill(z, evTime - trackTime);
      }
    }
  }

  canv = new TCanvas("time","tine",1000,600);
  canv->Divide(3,2);
  canv->cd(1);
  histoA->Draw();
  canv->cd(2);
  auto hax = histoA->ProjectionX();
  hax->SetTitle("A: Z @ X = 0");
  float meanAX = hax->GetBinCenter(hax->GetMaximumBin());
  hax->Fit("gaus","","", meanAX-0.5,meanAX+0.5);
  canv->cd(3);
  auto hay = histoA->ProjectionY();
  float meanAY = hay->GetBinCenter(hay->GetMaximumBin());
  hay->SetTitle("A: #Delta T");
  hay->Fit("gaus","","",meanAY-0.5,meanAY+0.5);

  canv->cd(4);
  histoC->Draw();
  canv->cd(5);
  auto hcx = histoA->ProjectionX();
  float meanCX = hcx->GetBinCenter(hcx->GetMaximumBin());
  hcx->SetTitle("C: Z @ X = 0");
  hcx->Fit("gaus","","",meanCX-0.5,meanCX+0.5);
  canv->cd(6);
  auto hcy = histoA->ProjectionY();
  float meanCY = hcy->GetBinCenter(hcy->GetMaximumBin());
  hcy->SetTitle("A: #Delta T");
  hcy->Fit("gaus","","",meanCY-0.5,meanCY+0.5);

  
}
