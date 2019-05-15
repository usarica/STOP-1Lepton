#include "EventAnalyzer.h"
#include "common_includes.h"
#include "FileTransferHelpers.h"


using namespace FileTransferHelpers;


void doGenProduction(std::string stropts){
  TDirectory* curdir = gDirectory;

  FrameworkOptionParser opts(stropts);

  gSystem->Exec(Form("mkdir -p %s", opts.outputDir().c_str()));
  TFile* foutput = TFile::Open((opts.outputDir()+opts.outputFilename()).c_str(), "recreate");
  BaseTree* outtree = new BaseTree("SelectedTree"); // The tree to record into the ROOT file
  curdir->cd();

  FrameworkSet theSet(opts, CMS4_EVENTS_TREE_NAME);

  WeightsHandler wgtHandler;
  wgtHandler.set2016SchemeFlag(true);
  for (auto* tree:theSet.getFrameworkTreeList()) wgtHandler.bookBranches(tree);

  GenInfoHandler genInfoHandler;
  for (auto* tree:theSet.getFrameworkTreeList()) genInfoHandler.bookBranches(tree);

  JetMETHandler jetHandler; // Needed for gen. jets
  jetHandler.setGenJetsFlag(true);
  jetHandler.setAK4JetsFlag(false);
  jetHandler.setAK8JetsFlag(false);
  jetHandler.setMETFlag(false);
  jetHandler.setTopsFlag(false);
  for (auto* tree:theSet.getFrameworkTreeList()) jetHandler.bookBranches(tree);

  EventAnalyzer analyzer(&theSet);
  if (opts.isData()) analyzer.setSampleIdStorageOption(FrameworkTreeLooperBase::kStoreByRunAndEventNumber);

  // Set maximum events to process
  analyzer.setMaximumEvents(opts.maxEventsToProcess());
  analyzer.setRecordEveryNEvents(opts.recordEveryNEvents());
  analyzer.setJetMETFlag(true);
  analyzer.setEventFilterFlag(false);
  analyzer.setVertexPUInfoFlag(false);
  analyzer.setPFCandsFlag(false);
  analyzer.setMuonsFlag(false);
  analyzer.setElectronsFlag(false);
  analyzer.setPhotonsFlag(false);
  analyzer.setCorrectedMETFlag(false);
  analyzer.setIsoTracksFlag(false);
  analyzer.setTausFlag(false);
  analyzer.setRecordIsoTracksFlag(false);
  analyzer.setRecordTausFlag(false);
  analyzer.setWriteFailingObjects(false);
  analyzer.setWriteSelectionVariables(false);
  if (opts.recordEveryNEvents()){
    BaseTree::setRobustSaveWrite(true);
    outtree->setAutoSave(0);
    outtree->setAutoFlush(0);
  }
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

  MELAout << "doGenProduction: File generation was successful! Initiating copy..." << endl;

  if (opts.condorOutputDir()!="") InitiateCondorFileTransfer(opts.outputDir().c_str(), opts.outputFilename().c_str(), opts.condorSite().c_str(), opts.condorOutputDir().c_str(), "", 5);
}
