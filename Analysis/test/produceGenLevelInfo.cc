#include "EventAnalyzer.h"


void produceGenLevelInfo(){
  TDirectory* curdir = gDirectory;
  std::vector<std::string> optlist={
    //"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=WZZ_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=-1 ismc=true",
    //"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root sample=/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM year=2018 maxevents=500000 ismc=true",
    //"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root sample=/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM year=2018 maxevents=500000 ismc=true",
    "indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=ZZTo4L_TuneCP5_13TeV_powheg_pythia8.root sample=/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM year=2018 maxevents=500000 ismc=true"
  };
  for (std::string const& stropts:optlist){
    FrameworkOptionParser opts(stropts);

    TFile* foutput = TFile::Open((opts.outputDir()+opts.outputFilename()).c_str(), "recreate");
    BaseTree* outtree = new BaseTree("test"); // The tree to record into the ROOT file
    curdir->cd();

    FrameworkSet theSet(opts, CMS4_EVENTS_TREE_NAME);

    WeightsHandler wgtHandler;
    wgtHandler.set2016SchemeFlag(true);
    for (auto* tree:theSet.getFrameworkTreeList()) wgtHandler.bookBranches(tree);

    GenInfoHandler genInfoHandler;
    for (auto* tree:theSet.getFrameworkTreeList()) genInfoHandler.bookBranches(tree);

    JetMETHandler jetHandler; // Needed for gen. jets
    jetHandler.setAK4JetsFlag(false);
    jetHandler.setAK8JetsFlag(false);
    jetHandler.setMETFlag(false);
    jetHandler.setTopsFlag(false);
    for (auto* tree:theSet.getFrameworkTreeList()) jetHandler.bookBranches(tree);

    EventAnalyzer analyzer(&theSet);
    analyzer.setEventFilterFlag(false);
    analyzer.setElectronsFlag(false);
    analyzer.setMuonsFlag(false);
    //analyzer.setJetMETFlag(false);

    // Set maximum events to process
    analyzer.setMaximumEvents(opts.maxEventsToProcess());
    // Ivy handlers
    analyzer.addExternalIvyObject("WeightsHandler", &wgtHandler);
    analyzer.addExternalIvyObject("GenInfoHandler", &genInfoHandler);
    analyzer.addExternalIvyObject("JetMETHandler", &jetHandler);
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
