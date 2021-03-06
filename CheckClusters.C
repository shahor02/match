/// \file CheckDigits.C
/// \brief Simple macro to check ITSU clusters

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TString.h>

#include "ITSMFTSimulation/Hit.h"
#include "MathUtils/Utils.h"
#include "MathUtils/Cartesian3D.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTReconstruction/Cluster.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#endif

void CheckClusters() {
  using namespace o2;
  using namespace o2::ITS;

  using o2::ITSMFT::Hit;
  using o2::ITSMFT::Cluster;

  TFile *f=TFile::Open("CheckClusters.root","recreate");
  TNtuple *nt=new TNtuple("ntc","cluster ntuple","x:y:z:xh:yh:zh:dx:dz:lab:rof:ev:hlx:hlz:clx:clz");

  char filename[100];

  // Geometry
  sprintf(filename, "o2sim_par.root");  TFile *file = TFile::Open(filename);
  gFile->Get("FairGeoParSet");
  
  auto gman =  o2::ITS::GeometryTGeo::Instance();
  gman->fillMatrixCache( utils::bit2Mask(TransformType::T2L, TransformType::T2GRot, TransformType::L2G) ); // request cached transforms
  
  // Hits
  sprintf(filename, "o2sim.root");
  TFile *file0 = TFile::Open(filename);
  TTree *hitTree=(TTree*)gFile->Get("o2sim");
  std::vector<Hit> *hitArray = nullptr;
  hitTree->SetBranchAddress("ITSHit", &hitArray);

  
  // Clusters
  sprintf(filename, "o2clus.root");
  TFile *file1 = TFile::Open(filename);
  TTree *clusTree=(TTree*)gFile->Get("o2sim");
  std::vector<Cluster> *clusArr=nullptr;
  //clusTree->SetBranchAddress("ITSCluster",&clusArr); // Why this does not work ???
  auto *branch = clusTree->GetBranch("ITSCluster");
  if (!branch) {
    std::cout<<"No clusters !"<<std::endl;
    return;
  }
  branch->SetAddress(&clusArr);
  // Cluster MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel> *clusLabArr=nullptr;
  clusTree->SetBranchAddress("ITSClusterMCTruth",&clusLabArr);

  Int_t nevCl = clusTree->GetEntries(); // clusters in cont. readout may be grouped as few events per entry
  Int_t nevH = hitTree->GetEntries(); // hits are stored as one event per entry
  int ievC=0,ievH=0;
  int lastReadHitEv = -1;
  for (ievC=0;ievC<nevCl;ievC++) {
    clusTree->GetEvent(ievC);
    Int_t nc = clusArr->size();
    printf("processing cluster event %d\n",ievC);

    while(nc--) {
      // cluster is in tracking coordinates always
      Cluster &c=(*clusArr)[nc];
      Int_t chipID = c.getSensorID();
      const auto locC = c.getXYZLoc(*gman); // convert from tracking to local frame
      const auto gloC = c.getXYZGloRot(*gman); // convert from tracking to global frame
      auto lab = (clusLabArr->getLabels(nc))[0];

      float dx=0,dz=0;
      int trID = lab.getTrackID();
      int ievH = lab.getEventID();
      
      Point3D<float> locH,locHsta,gloH;
      if (trID>=0) { // is this cluster from hit or noise ?  
	//printf("Ev:%d Cl%d/%d RO:%d Chip:%d -> H %d / %d lastReadHitEv: %d\n",ievC,nc,clusArr->size(),c.getROFrame(),chipID, trID, ievH, lastReadHitEv);
	Hit* p = nullptr;
	if (lastReadHitEv!=ievH) {
	  hitTree->GetEvent(ievH);
	  lastReadHitEv = ievH;
	}
	for (auto& ptmp : *hitArray) {
	  if (ptmp.GetDetectorID() != chipID) continue; 
	  if (ptmp.GetTrackID() != trID) continue;
	  p = &ptmp;
	  break;
	}
	if (!p) {
	  printf("did not find hit (scanned HitEvs %d %d) for cluster of tr%d on chip %d\n",ievH,nevH,trID,chipID);
	  locH.SetXYZ(0.f,0.f,0.f);
	}
	else {
	  // mean local position of the hit
	  locH    = gman->getMatrixL2G(chipID)^( p->GetPos() );  // inverse conversion from global to local
	  locHsta = gman->getMatrixL2G(chipID)^( p->GetPosStart() );
	  locH.SetXYZ( 0.5*(locH.X()+locHsta.X()),0.5*(locH.Y()+locHsta.Y()),0.5*(locH.Z()+locHsta.Z()) );
	  gloH =  p->GetPos();
	  //std::cout << "chip "<< p->GetDetectorID() << "  PposGlo " << p->GetPos() << std::endl;
	  //std::cout << "chip "<< c->getSensorID() << "  PposLoc " << locH << std::endl;
	  dx = locH.X()-locC.X();
	  dz = locH.Z()-locC.Z();
	}
      }
      nt->Fill(gloC.X(),gloC.Y(),gloC.Z(),
	       gloH.X(),gloH.Y(),gloH.Z(),
	       dx, dz, trID, c.getROFrame(), ievC,
	       locH.X(),locH.Z(), locC.X(),locC.Z());
    }
  }
  new TCanvas; nt->Draw("y:x");
  new TCanvas; nt->Draw("dx:dz");
  f->cd();
  nt->Write();
  f->Close();
}
