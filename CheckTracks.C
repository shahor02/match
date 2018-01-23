/// \file CheckTracks.C
/// \brief Simple macro to check ITSU tracks

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <array>

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>

#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "ITSMFTReconstruction/Cluster.h"
#include "ITSReconstruction/CookedTrack.h"

#include "CommonUtils/TreeStream.h"
#include "CommonUtils/TreeStreamRedirector.h"

#endif

using namespace o2::utils;

void CheckTracks(std::string paramfile="o2sim.root",
		 std::string inputfile="o2track_its.root") {
  using namespace o2::ITS;

  TreeStreamRedirector mDBGOut("CheckTracks.root","recreate");
  

  char filename[100];

  // MC tracks
  TFile *file0 = TFile::Open(paramfile.data());
  TTree *mcTree=(TTree*)gFile->Get("o2sim");
  std::vector<o2::MCTrack>* mcArr = nullptr;
  mcTree->SetBranchAddress("MCTrack",&mcArr);

  // Reconstructed tracks
  TFile *file1 = TFile::Open(inputfile.data());
  TTree *recTree=(TTree*)gFile->Get("o2sim");
  std::vector<CookedTrack> *recArr=nullptr;
  recTree->SetBranchAddress("ITSTrack",&recArr);
  // Track MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> *trkLabArr=nullptr;
  recTree->SetBranchAddress("ITSTrackMCTruth",&trkLabArr);

  
  Int_t nevMC = mcTree->GetEntries();
  Int_t nevRec = recTree->GetEntries();
  
  for (Int_t n=0; n<nevMC; n++) {
    std::cout<<"Event "<<n<<'/'<<nevMC<<std::endl;
    Int_t nGen=0, nGoo=0;
    mcTree->GetEvent(n);
    Int_t nmc=mcArr->size();

    while(nmc--) {
      const auto& mcTrack = (*mcArr)[nmc];
      Int_t mID = mcTrack.getMotherTrackId();
      if (mID >= 0) continue; // Select primary particles 
      Int_t pdg = mcTrack.GetPdgCode();
      if (TMath::Abs(pdg) != 211) continue;  // Select pions
      
      nGen++; // Generated tracks for the efficiency calculation 
      
      Float_t mcPx = mcTrack.GetStartVertexMomentumX();
      Float_t mcPy = mcTrack.GetStartVertexMomentumY();
      Float_t mcPz = mcTrack.GetStartVertexMomentumZ();
      Float_t mcPt = mcTrack.GetPt();
      Float_t mcPhi= TMath::ATan2(mcPy,mcPx);
      Float_t mcLam= TMath::ATan2(mcPz,mcPt);
      Float_t recPhi=-1.; 
      Float_t recLam=-1.; 
      Float_t recPt=-1.; 
      Float_t ip[2]{0.,0.};
      o2::MCCompLabel label;
      int roFrame = -1;
      for (int ievRec=n;ievRec<nevRec;ievRec++) { // corresponding rec event may start fron n
	recTree->GetEvent(ievRec);
	Int_t nrec=recArr->size();	
	for (Int_t i=0; i<nrec; i++) {
	  const CookedTrack &recTrack = (*recArr)[i];
	  auto mclab = (trkLabArr->getLabels(i))[0];
	  if (mclab.getEventID()!=n || std::abs(mclab.getTrackID())!=nmc) continue;
	  label = mclab;
	  std::array<float,3> p;
	  recTrack.getPxPyPzGlo(p);
	  recPt = recTrack.getPt();
	  recPhi = TMath::ATan2(p[1],p[0]);
	  recLam = TMath::ATan2(p[2],recPt);
	  Float_t vx=0., vy=0., vz=0.;  // Assumed primary vertex
	  Float_t bz=5.;                // Assumed magnetic field 
	  recTrack.getImpactParams(vx, vy, vz, bz, ip);
	  roFrame = recTrack.getROFrame();
	  if (label.getTrackID()==nmc) nGoo++; // Good found tracks for the efficiency calculation
	  break;
	}
	if (!label.isEmpty()) break; // was assigned
      }
      mDBGOut<<"itscheck"
	     <<"mcPhi="<<mcPhi<<"mcLam="<<mcLam<<"mcPt="<<mcPt<<"recPhi="<<recPhi<<"recLam="<<recLam<<"recPt="<<recPt
	     <<"ipD="<<ip[0]<<"ipZ="<<ip[1]<<"roFrame="<<roFrame<<"label="<<label<<"\n";
    }
    Float_t eff = (nGen > 0) ? nGoo/Float_t(nGen) : -1.;
    std::cout<<"Good found tracks: "<<nGoo<<" gen tracks "<<nGen<<",  efficiency: "<<eff<<std::endl;
  }


  // "recPt>0" means "found tracks only"  
  // "label>0" means "found good tracks only"
  //
  //  TFile fltree("CheckTracks.root");
  //  TTree* tree = (TTree*)fltree.Get("itscheck");
  TTree* tree = &(mDBGOut<<"itscheck").getTree();
  new TCanvas; tree->Draw("ipD","recPt>0 && label.getTrackID()>=0");
  new TCanvas; tree->Draw("mcLam-recLam","recPt>0 && label.getTrackID()>=0");
  new TCanvas; tree->Draw("mcPt-recPt","recPt>0 && label.getTrackID()>=0");

  mDBGOut.Close();

  
}
