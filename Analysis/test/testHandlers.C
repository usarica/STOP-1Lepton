#include "EventAnalyzer.h"


void testHandlers(){
  TDirectory* curdir = gDirectory;

  //std::string stropts = "indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./ outfile=WZZ_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=100 ismc=true";
  //std::string stropts = "indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./ outfile=WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=1000 ismc=true";
  //std::string stropts = "indir=/hadoop/cms/store/group/snt/run2_mc2017 outdir=./ outfile=GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8_0.root sample=/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM year=2017 maxevents=1000 ismc=true inputfiles=merged_ntuple_1.root";
  std::string stropts = "indir=/hadoop/cms/store/group/snt/run2_data2018 outdir=./ outfile=Run2018D-PromptReco-v2_MET.root sample=/MET/Run2018D-PromptReco-v2/MINIAOD inputfiles=merged_ntuple_357.root year=2018 maxevents=1000 ismc=false";
  //std::string stropts = "indir=/hadoop/cms/store/group/snt/run2_data2018 outdir=./ outfile=Run2018C-17Sep2018-v1_JetHT.root sample=/JetHT/Run2018C-17Sep2018-v1/MINIAOD inputfiles=merged_ntuple_19.root year=2018 maxevents=1000 ismc=false";
  //std::string stropts = "indir=/hadoop/cms/store/group/snt/run2_data2018 outdir=./ outfile=Run2018C-17Sep2018-v1_MET.root sample=/MET/Run2018C-17Sep2018-v1/MINIAOD inputfiles=merged_ntuple_19.root year=2018 maxevents=1000 ismc=false";
  //std::string stropts = "indir=/hadoop/cms/store/group/snt/run2_data2017 outdir=./ outfile=Run2017F-31Mar2018-v1_DoubleEG.root sample=/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD year=2017 maxevents=1000 ismc=false inputfiles=merged_ntuple_1.root";
  //std::string stropts = "indir=/hadoop/cms/store/group/snt/run2_data2017 outdir=./ outfile=Run2017F-09May2018-v1_DoubleEG.root sample=/DoubleEG/Run2017F-09May2018-v1/MINIAOD year=2017 maxevents=1000 ismc=false inputfiles=merged_ntuple_1.root";

  FrameworkOptionParser opts(stropts);

  TFile* foutput = TFile::Open((opts.outputDir()+opts.outputFilename()).c_str(), "recreate");
  BaseTree* outtree = new BaseTree("test"); // The tree to record into the ROOT file
  curdir->cd();

  FrameworkSet theSet(opts, CMS4_EVENTS_TREE_NAME);

  WeightsHandler wgtHandler;
  //wgtHandler.setVerbosity(TVar::DEBUG);
  for (auto* tree:theSet.getFrameworkTreeList()) wgtHandler.bookBranches(tree);

  GenInfoHandler genInfoHandler;
  //genInfoHandler.setVerbosity(TVar::DEBUG);
  for (auto* tree:theSet.getFrameworkTreeList()) genInfoHandler.bookBranches(tree);

  EventFilterHandler eventFilter;
  //eventFilter.setVerbosity(TVar::DEBUG_VERBOSE);
  for (auto* tree:theSet.getFrameworkTreeList()) eventFilter.bookBranches(tree);

  PUScaleFactorHandler puSFHandler;
  VertexPUHandler vertexPUHandler;
  //vertexPUHandler.setVerbosity(TVar::DEBUG_VERBOSE);
  for (auto* tree:theSet.getFrameworkTreeList()) vertexPUHandler.bookBranches(tree);

  PFCandHandler pfcandHandler;
  //pfcandHandler.setVerbosity(TVar::DEBUG);
  for (auto* tree:theSet.getFrameworkTreeList()) pfcandHandler.bookBranches(tree);

  MuonScaleFactorHandler muonSFHandler;
  MuonHandler muonHandler;
  //muonHandler.setVerbosity(TVar::DEBUG);
  for (auto* tree:theSet.getFrameworkTreeList()) muonHandler.bookBranches(tree);

  ElectronScaleFactorHandler electronSFHandler;
  ElectronHandler electronHandler;
  //electronHandler.setVerbosity(TVar::DEBUG);
  for (auto* tree:theSet.getFrameworkTreeList()) electronHandler.bookBranches(tree);

  PhotonScaleFactorHandler photonSFHandler;
  PhotonHandler photonHandler;
  //photonHandler.setVerbosity(TVar::DEBUG);
  for (auto* tree:theSet.getFrameworkTreeList()) photonHandler.bookBranches(tree);

  JetMETHandler jetHandler;
  //jetHandler.setVerbosity(TVar::DEBUG);
  METCorrectionHandler metCorrector;
  BtagScaleFactorHandler btagSFHandler_MC_noFS(BtagHelpers::kDeepCSV_Medium, false);
  BtagScaleFactorHandler btagSFHandler_MC_FS(BtagHelpers::kDeepCSV_Medium, true);
  jetHandler.registerBtagSFHandlers(&btagSFHandler_MC_noFS, &btagSFHandler_MC_FS);
  JECScaleFactorHandler jecSFHandler_ak4(JECJERHelpers::kAK4);
  JECScaleFactorHandler jecSFHandler_ak8(JECJERHelpers::kAK8);
  jetHandler.registerJECSFHandlers(&jecSFHandler_ak4, &jecSFHandler_ak8);
  JERScaleFactorHandler jerSFHandler_ak4(JECJERHelpers::kAK4);
  JERScaleFactorHandler jerSFHandler_ak8(JECJERHelpers::kAK8);
  jetHandler.registerJERSFHandlers(&jerSFHandler_ak4, &jerSFHandler_ak8);
  for (auto* tree:theSet.getFrameworkTreeList()) jetHandler.bookBranches(tree);

  IsoTrackHandler isotrkHandler;
  //isotrkHandler.setVerbosity(TVar::DEBUG);
  for (auto* tree:theSet.getFrameworkTreeList()) isotrkHandler.bookBranches(tree);

  TauHandler tauHandler;
  //tauHandler.setVerbosity(TVar::DEBUG);
  for (auto* tree:theSet.getFrameworkTreeList()) tauHandler.bookBranches(tree);

  EventAnalyzer analyzer(&theSet);
  if (opts.isData()) analyzer.setSampleIdStorageOption(FrameworkTreeLooperBase::kStoreByRunAndEventNumber);
  // Set maximum events to process
  analyzer.setMaximumEvents(opts.maxEventsToProcess());
  // Control what is recorded
  //analyzer.setRecordIsoTracksFlag(false);
  //analyzer.setRecordTausFlag(false);
  //analyzer.setWriteFailingObjects(false);
  //analyzer.setWriteSelectionVariables(false);
  // Ivy handlers
  analyzer.addExternalIvyObject("WeightsHandler", &wgtHandler);
  analyzer.addExternalIvyObject("GenInfoHandler", &genInfoHandler);
  analyzer.addExternalIvyObject("EventFilterHandler", &eventFilter);
  analyzer.addExternalIvyObject("VertexPUHandler", &vertexPUHandler);
  analyzer.addExternalIvyObject("PFCandHandler", &pfcandHandler);
  analyzer.addExternalIvyObject("MuonHandler", &muonHandler);
  analyzer.addExternalIvyObject("ElectronHandler", &electronHandler);
  analyzer.addExternalIvyObject("PhotonHandler", &photonHandler);
  analyzer.addExternalIvyObject("JetMETHandler", &jetHandler);
  analyzer.addExternalIvyObject("IsoTrackHandler", &isotrkHandler);
  analyzer.addExternalIvyObject("TauHandler", &tauHandler);
  // SF handlers
  analyzer.addExternalScaleFactorHandler("MuonSFHandler", &muonSFHandler);
  analyzer.addExternalScaleFactorHandler("ElectronSFHandler", &electronSFHandler);
  analyzer.addExternalScaleFactorHandler("PhotonSFHandler", &photonSFHandler);
  analyzer.addExternalScaleFactorHandler("PileUpSFHandler", &puSFHandler);
  // Register JetMET SF handlers so that they can be updated
  analyzer.addExternalScaleFactorHandler("BTagSFHandler_MC_noFS", &btagSFHandler_MC_noFS);
  analyzer.addExternalScaleFactorHandler("BTagSFHandler_MC_FS", &btagSFHandler_MC_FS);
  analyzer.addExternalScaleFactorHandler("JECSFHandler_ak4", &jecSFHandler_ak4);
  analyzer.addExternalScaleFactorHandler("JECSFHandler_ak8", &jecSFHandler_ak8);
  analyzer.addExternalScaleFactorHandler("JERSFHandler_ak4", &jerSFHandler_ak4);
  analyzer.addExternalScaleFactorHandler("JERSFHandler_ak8", &jerSFHandler_ak8);
  analyzer.addExternalScaleFactorHandler("METCorrectionHandler", &metCorrector);
  // Output tree setup
  analyzer.setExternalProductTree(outtree);

  analyzer.loop(true, false, true);
  MELAout << "There are " << outtree->getNEvents() << " products" << endl;
  outtree->writeToFile(foutput);
  delete outtree;

  foutput->Close();
}
