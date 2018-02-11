void proc(const char* dir="./", bool dbg = false)
{
  if (dbg) gROOT->ProcessLine(".x setAClicDBG.C");
  if (dbg) {
    gROOT->ProcessLine(".L MatchTPCITS.cxx++g");
    gROOT->ProcessLine(Form(".x testMatch.C++g(\"%s\")",dir));
  }
  else {
    gROOT->ProcessLine(".L MatchTPCITS.cxx++");
    gROOT->ProcessLine(Form(".x testMatch.C++(\"%s\")",dir));
  }
}
