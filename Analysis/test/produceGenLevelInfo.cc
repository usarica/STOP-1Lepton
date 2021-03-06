#include "EventAnalyzer.h"


void produceGenLevelInfo(bool doReco=false, int whichSample=-1){
  TDirectory* curdir = gDirectory;
  std::vector<std::string> const optlist={
    /*0*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.root sample=/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM year=2018 maxevents=500000 ismc=true",
    /*1*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root sample=/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM year=2018 maxevents=500000 ismc=true",
    /*2*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8.root sample=/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=500000 ismc=true",
    /*3*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=ZZTo4L_TuneCP5_13TeV_powheg_pythia8.root sample=/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM year=2018 maxevents=500000 ismc=true",
    /*4*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=-1 ismc=true",
    /*5*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=WZZ_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=-1 ismc=true",
    /*6*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=TT_DiLept_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/TT_DiLept_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=500000 ismc=true",
    /*7*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=TTPlus1Jet_DiLept_TuneCP5_13TeV-amcatnloFXFX-pythia8.root sample=/TTPlus1Jet_DiLept_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=500000 ismc=true",
    /*8*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root sample=/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=500000 ismc=true",
    /*9*/"indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./output/test_GenWeights/2018/ outfile=TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=500000 ismc=true"
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
    for (auto* tree:theSet.getFrameworkTreeList()) genInfoHandler.bookBranches(tree);

    MuonHandler muonHandler;
    if (doReco){ for (auto* tree:theSet.getFrameworkTreeList()) muonHandler.bookBranches(tree); }

    ElectronHandler electronHandler;
    if (doReco){ for (auto* tree:theSet.getFrameworkTreeList()) electronHandler.bookBranches(tree); }

    PhotonHandler photonHandler;
    if (doReco){ for (auto* tree:theSet.getFrameworkTreeList()) photonHandler.bookBranches(tree); }

    JetMETHandler jetHandler; // Needed for gen. jets
    if (!doReco){
      jetHandler.setAK4JetsFlag(false);
      jetHandler.setAK8JetsFlag(false);
      jetHandler.setMETFlag(false);
      jetHandler.setTopsFlag(false);
    }
    for (auto* tree:theSet.getFrameworkTreeList()) jetHandler.bookBranches(tree);

    EventAnalyzer analyzer(&theSet);
    analyzer.setEventFilterFlag(false);
    if (!doReco){
      analyzer.setMuonsFlag(false);
      analyzer.setElectronsFlag(false);
      analyzer.setPhotonsFlag(false);
    }

    // Set maximum events to process
    analyzer.setMaximumEvents(opts.maxEventsToProcess());
    // Ivy handlers
    analyzer.addExternalIvyObject("WeightsHandler", &wgtHandler);
    analyzer.addExternalIvyObject("GenInfoHandler", &genInfoHandler);
    if (doReco){
      analyzer.addExternalIvyObject("MuonHandler", &muonHandler);
      analyzer.addExternalIvyObject("ElectronHandler", &electronHandler);
      analyzer.addExternalIvyObject("PhotonHandler", &photonHandler);
    }
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
