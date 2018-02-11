void proc(const char* dir="./", float itsOffs=0, bool dbg = false)
{
  if (dbg) gROOT->ProcessLine(".x setAClicDBG.C");
  if (dbg) {
    gROOT->ProcessLine(".L MatchTPCITS.cxx++g");
    gROOT->ProcessLine(Form(".x testMatch.C++g(\"%s\",%f)",dir,itsOffs));
  }
  else {
    gROOT->ProcessLine(".L MatchTPCITS.cxx++");
    gROOT->ProcessLine(Form(".x testMatch.C++(\"%s\",%f)",dir,itsOffs));
  }
}
