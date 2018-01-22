{
  printf("Executing %s/rootlogon.C\n",gSystem->pwd());
  //  gROOT->Macro("$O2_ROOT/share/macro/load_all_libs.C");
  gSystem->Load("libSimulationDataFormat.so");

  TString so = gSystem->GetMakeSharedLib();
  so.ReplaceAll("-O2","-O0");
  gSystem->SetMakeSharedLib(so.Data());
  printf("Set MakeSharedLib: %s\n", gSystem->GetMakeSharedLib());
}
