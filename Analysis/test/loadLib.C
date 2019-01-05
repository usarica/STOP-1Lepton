{
  TString LIBCMS4CORE="CMS3_CORE.so";

  gSystem->Load("$CMSSW_BASE/src/CMSDataTools/AnalysisTree/test/loadLib.C");

  // Top tagger
  gSystem->Load("libTopTaggerDataFormats.so");
  gSystem->Load("libTopTaggerCfgParser.so");
  gSystem->Load("libTopTaggerTopTagger.so");

  gSystem->AddIncludePath("-I$CMSSW_BASE/src/STOP_1Lepton/Analysis/interface/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/STOP_1Lepton/Analysis/test/");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/cmstas/CORE/");
  gSystem->Load("$CMSSW_BASE/src/cmstas/CORE/" + LIBCMS4CORE);
  gSystem->Load("libSTOP_1LeptonAnalysis.so");
}
