void proc(const char* dir="./")
{
  gROOT->ProcessLine(".x setAClicDBG.C");
  gROOT->ProcessLine(".L MatchTPCITS.cxx++g");
  gROOT->ProcessLine(Form(".x testMatch.C++g(\"%s\")",dir));
}
