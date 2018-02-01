{
  printf("Executing %s/rootlogon.C\n",gSystem->pwd());

  printf("Executing ~/rootlogon.C\n");
  gROOT->ProcessLine(".x ~/rootlogon.C");
  
  //  gROOT->Macro("$O2_ROOT/share/macro/load_all_libs.C");
  gSystem->Load("libSimulationDataFormat.so");

  /*
  TString so = gSystem->GetMakeSharedLib();
  so.ReplaceAll("-O2","-O0");
  gSystem->SetMakeSharedLib(so.Data());
  printf("Set MakeSharedLib: %s\n", gSystem->GetMakeSharedLib());
  */
}
