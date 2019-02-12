#include "EventAnalyzer.h"


void doMTDProduction(int whichSample=-1){
  TDirectory* curdir = gDirectory;
  std::vector<std::string> const optlist={
    /*0*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/MTD/2018/ outfile=QCD_Pt-15to20_EMEnriched_TuneCP5_13TeV_pythia8.root sample=/QCD_Pt-15to20_EMEnriched_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=500000 recordeveryn=50000 ismc=true",
    /*1*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/MTD/2018/ outfile=QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8.root sample=/QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM year=2018 maxevents=500000 recordeveryn=50000 ismc=true",
    /*2*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/MTD/2018/ outfile=QCD_Pt-170to300_EMEnriched_TuneCP5_13TeV_pythia8.root sample=/QCD_Pt-170to300_EMEnriched_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM year=2018 maxevents=500000 recordeveryn=50000 ismc=true",
    /*3*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/MTD/2018/ outfile=QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8.root sample=/QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM year=2018 maxevents=500000 recordeveryn=50000 ismc=true"
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
    //analyzer.setWeightsFlag(false);
    //analyzer.setGenInfoFlag(false);
    //analyzer.setGenParticlesFlag(false);
    analyzer.setEventFilterFlag(false);
    analyzer.setPFCandsFlag(false);
    analyzer.setMuonsFlag(false);
    analyzer.setElectronsFlag(false);
    analyzer.setPhotonsFlag(false);
    //analyzer.setJetMETFlag(false);
    analyzer.setWriteSelectionVariables(false);
    // Ivy handlers
    analyzer.addExternalIvyObject("WeightsHandler", &wgtHandler);
    analyzer.addExternalIvyObject("GenInfoHandler", &genInfoHandler);
    analyzer.addExternalIvyObject("JetMETHandler", &jetHandler);
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
