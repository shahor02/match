#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>

#include "Rtypes.h"
#include "TString.h"
#include "TStopwatch.h"
#include "TGeoManager.h"

#include "FairLogger.h"
#include "FairRunAna.h"
#include "FairFileSource.h"
#include "FairSystemInfo.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"

#include "TPCReconstruction/ClustererTask.h"
#include "ITSReconstruction/ClustererTask.h"

#endif

unsigned threadsTPC = 0;

void run_clus(std::string outputfile="o2clus.root"
	      ,std::string inputfile="o2dig.root"	     
	      ,std::string paramfile="o2sim_par.root"
	      )
{

  // Initialize logger
  FairLogger *logger = FairLogger::GetLogger();
  logger->SetLogVerbosityLevel("LOW");
  logger->SetLogScreenLevel("INFO");
  
  TStopwatch timer;
  
  // Setup FairRoot analysis manager
  FairRunAna * run = new FairRunAna();
  FairFileSource *fFileSource = new FairFileSource(inputfile.data());
  run->SetSource(fFileSource);
  run->SetOutputFile(outputfile.data());
  
  // Setup Runtime DB
  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  FairParRootFileIo* parInput1 = new FairParRootFileIo();
  parInput1->open(paramfile.data());
  rtdb->setFirstInput(parInput1);

  //============================= create clusterizers =============================>>>

  o2::ITS::ClustererTask *clus = new o2::ITS::ClustererTask();
  run->AddTask(clus);

  o2::TPC::ClustererTask *clustTPC = new o2::TPC::ClustererTask();
  clustTPC->setContinuousReadout(true);
  clustTPC->setClustererEnable(o2::TPC::ClustererTask::ClustererType::Box,false);
  clustTPC->setClustererEnable(o2::TPC::ClustererTask::ClustererType::HW,true);
  
  run->AddTask(clustTPC);

  //-----------
  run->Init();
  //===============================================================================<<<
  //  clustTPC->getHwClusterer()->setProcessingType(o2::TPC::HwClusterer::Processing::Sequential);
  clustTPC->getHwClusterer()->setNumThreads(threadsTPC);
  
  // Start simulation
  timer.Start();
  run->Run();

  std::cout << std::endl << std::endl;

  // Extract the maximal used memory an add is as Dart measurement
  // This line is filtered by CTest and the value send to CDash
  FairSystemInfo sysInfo;
  Float_t maxMemory=sysInfo.GetMaxMemory();
  std::cout << R"(<DartMeasurement name="MaxMemory" type="numeric/double">)";
  std::cout << maxMemory;
  std::cout << "</DartMeasurement>" << std::endl;
  
  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  
  Float_t cpuUsage=ctime/rtime;
  std::cout << R"(<DartMeasurement name="CpuLoad" type="numeric/double">)";
  std::cout << cpuUsage;
  std::cout << "</DartMeasurement>" << std::endl;

  std::cout << std::endl << std::endl;
  std::cout << "Output file is "    << outputfile.data() << std::endl;
  //std::cout << "Parameter file is " << parFile << std::endl;
  std::cout << "Real time " << rtime << " s, CPU time " << ctime
	    << "s" << std::endl << std::endl;
  std::cout << "Macro finished succesfully." << std::endl;
  
  
}
