#include "EventAnalyzer.h"


void doFullProduction(int whichSample=-1){
  TDirectory* curdir = gDirectory;
  std::vector<std::string> const optlist={
    /*0*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root sample=/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM year=2018 maxevents=500000 recordeveryn=50000 ismc=true",
    /*1*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root sample=/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM year=2018 maxevents=500000 recordeveryn=50000 ismc=true",
    /*2*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8.root sample=/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=500000 recordeveryn=50000 ismc=true",
    /*3*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=ZZTo4L_TuneCP5_13TeV_powheg_pythia8.root sample=/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM year=2018 maxevents=500000 recordeveryn=50000 ismc=true",
    /*4*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=-1 recordeveryn=50000 ismc=true",
    /*5*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=WZZ_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=-1 recordeveryn=50000 ismc=true",
    /*6*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=TT_DiLept_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/TT_DiLept_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=500000 recordeveryn=50000 ismc=true",
    /*7*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=TTPlus1Jet_DiLept_TuneCP5_13TeV-amcatnloFXFX-pythia8.root sample=/TTPlus1Jet_DiLept_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=500000 recordeveryn=50000 ismc=true",
    /*8*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root sample=/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=500000 recordeveryn=50000 ismc=true",
    /*9*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=500000 recordeveryn=50000 ismc=true"
  };
  if (whichSample>=(int) optlist.size()) return;
  for (std::string const& stropts:optlist){
    if (whichSample>=0 && stropts!=optlist.at(whichSample)) continue;

    FrameworkOptionParser opts(stropts);

    TFile* foutput = TFile::Open((opts.outputDir()+opts.outputFilename()).c_str(), "recreate");
    BaseTree* outtree = new BaseTree("test"); // The tree to record into the ROOT file
    curdir->cd();

    FrameworkSet theSet(opts, CMS4_EVENTS_TREE_NAME);

    WeightsHandler wgtHandler;
    wgtHandler.set2016SchemeFlag(true);
    for (auto* tree:theSet.getFrameworkTreeList()) wgtHandler.bookBranches(tree);

    GenInfoHandler genInfoHandler;
    //genInfoHandler.setParticleInfoFlag(false);
    for (auto* tree:theSet.getFrameworkTreeList()) genInfoHandler.bookBranches(tree);

    EventFilterHandler eventFilter;
    //eventFilter.setVerbosity(TVar::DEBUG_VERBOSE);
    for (auto* tree:theSet.getFrameworkTreeList()) eventFilter.bookBranches(tree);

    PFCandHandler pfcandHandler;
    for (auto* tree:theSet.getFrameworkTreeList()) pfcandHandler.bookBranches(tree);

    MuonScaleFactorHandler muonSFHandler;
    MuonHandler muonHandler;
    for (auto* tree:theSet.getFrameworkTreeList()) muonHandler.bookBranches(tree);

    ElectronScaleFactorHandler electronSFHandler;
    ElectronHandler electronHandler;
    for (auto* tree:theSet.getFrameworkTreeList()) electronHandler.bookBranches(tree);

    //PhotonScaleFactorHandler photonSFHandler;
    PhotonHandler photonHandler;
    for (auto* tree:theSet.getFrameworkTreeList()) photonHandler.bookBranches(tree);

    JetMETHandler jetHandler; // Needed for gen. jets
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

    EventAnalyzer analyzer(&theSet);

    // Set maximum events to process
    analyzer.setMaximumEvents(opts.maxEventsToProcess());
    analyzer.setRecordEveryNEvents(opts.recordEveryNEvents());
    analyzer.setWriteSelectionVariables(false);
    // Ivy handlers
    analyzer.addExternalIvyObject("WeightsHandler", &wgtHandler);
    analyzer.addExternalIvyObject("GenInfoHandler", &genInfoHandler);
    analyzer.addExternalIvyObject("EventFilterHandler", &eventFilter);
    analyzer.addExternalIvyObject("PFCandHandler", &pfcandHandler);
    analyzer.addExternalIvyObject("MuonHandler", &muonHandler);
    analyzer.addExternalIvyObject("ElectronHandler", &electronHandler);
    analyzer.addExternalIvyObject("PhotonHandler", &photonHandler);
    analyzer.addExternalIvyObject("JetMETHandler", &jetHandler);
    // SF handlers
    analyzer.addExternalScaleFactorHandler("MuonSFHandler", &muonSFHandler);
    analyzer.addExternalScaleFactorHandler("ElectronSFHandler", &electronSFHandler);
    //analyzer.addExternalScaleFactorHandler("PhotonSFHandler", &photonSFHandler);
    // Register JetMET SF handlers so that they can be updated
    analyzer.addExternalScaleFactorHandler("BTagSFHandler_MC_noFS", &btagSFHandler_MC_noFS);
    analyzer.addExternalScaleFactorHandler("BTagSFHandler_MC_FS", &btagSFHandler_MC_FS);
    analyzer.addExternalScaleFactorHandler("JECSFHandler_ak4", &jecSFHandler_ak4);
    analyzer.addExternalScaleFactorHandler("JECSFHandler_ak8", &jecSFHandler_ak8);
    analyzer.addExternalScaleFactorHandler("JERSFHandler_ak4", &jerSFHandler_ak4);
    analyzer.addExternalScaleFactorHandler("JERSFHandler_ak8", &jerSFHandler_ak8);
    // Output tree setup
    analyzer.setExternalProductTree(outtree);

    analyzer.loop(true, false, true);
    MELAout << "There are " << outtree->getNEvents() << " products" << endl;
    outtree->writeToFile(foutput);
    delete outtree;

    foutput->Close();
    curdir->cd();
  }
}
