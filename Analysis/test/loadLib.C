{
  gSystem->Load("$CMSSW_BASE/src/CMSDataTools/AnalysisTree/test/loadLib.C");

  gSystem->AddIncludePath("-I$CMSSW_BASE/src/STOP_1Lepton/Analysis/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/STOP_1Lepton/Analysis/test/");
  gSystem->Load("libSTOP_1LeptonAnalysis.so");
}
