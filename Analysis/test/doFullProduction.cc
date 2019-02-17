#include "EventAnalyzer.h"
#include "common_includes.h"
#include "FileTransferHelpers.h"


#define _DOPFCANDS_ false

using namespace FileTransferHelpers;


void doFullProduction(std::string stropts){
  TDirectory* curdir = gDirectory;

  FrameworkOptionParser opts(stropts);

  gSystem->Exec(Form("mkdir -p %s", opts.outputDir().c_str()));
  TFile* foutput = TFile::Open((opts.outputDir()+opts.outputFilename()).c_str(), "recreate");
  BaseTree* outtree = new BaseTree("test"); // The tree to record into the ROOT file
  outtree->setAutoSave(0);
  curdir->cd();

  FrameworkSet theSet(opts, CMS4_EVENTS_TREE_NAME);

  WeightsHandler wgtHandler;
  wgtHandler.set2016SchemeFlag(true);
  for (auto* tree:theSet.getFrameworkTreeList()) wgtHandler.bookBranches(tree);

  GenInfoHandler genInfoHandler;
  for (auto* tree:theSet.getFrameworkTreeList()) genInfoHandler.bookBranches(tree);

  EventFilterHandler eventFilter;
  for (auto* tree:theSet.getFrameworkTreeList()) eventFilter.bookBranches(tree);

  PUScaleFactorHandler puSFHandler;
  VertexPUHandler vertexPUHandler;
  for (auto* tree:theSet.getFrameworkTreeList()) vertexPUHandler.bookBranches(tree);

  PFCandHandler pfcandHandler;
  if (_DOPFCANDS_){ for (auto* tree:theSet.getFrameworkTreeList()) pfcandHandler.bookBranches(tree); }

  MuonScaleFactorHandler muonSFHandler;
  MuonHandler muonHandler;
  for (auto* tree:theSet.getFrameworkTreeList()) muonHandler.bookBranches(tree);

  ElectronScaleFactorHandler electronSFHandler;
  ElectronHandler electronHandler;
  for (auto* tree:theSet.getFrameworkTreeList()) electronHandler.bookBranches(tree);

  PhotonScaleFactorHandler photonSFHandler;
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
  analyzer.setPFCandsFlag(_DOPFCANDS_);
  // Ivy handlers
  analyzer.addExternalIvyObject("WeightsHandler", &wgtHandler);
  analyzer.addExternalIvyObject("GenInfoHandler", &genInfoHandler);
  analyzer.addExternalIvyObject("EventFilterHandler", &eventFilter);
  analyzer.addExternalIvyObject("VertexPUHandler", &vertexPUHandler);
  if (_DOPFCANDS_) analyzer.addExternalIvyObject("PFCandHandler", &pfcandHandler);
  analyzer.addExternalIvyObject("MuonHandler", &muonHandler);
  analyzer.addExternalIvyObject("ElectronHandler", &electronHandler);
  analyzer.addExternalIvyObject("PhotonHandler", &photonHandler);
  analyzer.addExternalIvyObject("JetMETHandler", &jetHandler);
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
  // Output tree setup
  analyzer.setExternalProductTree(outtree);

  analyzer.loop(true, false, true);
  MELAout << "There are " << outtree->getNEvents() << " products" << endl;
  outtree->writeToFile(foutput);
  delete outtree;

  foutput->Close();
  curdir->cd();

  if (opts.condorOutputDir()!="") InitiateCondorFileTransfer(opts.outputDir().c_str(), opts.outputFilename().c_str(), opts.condorSite().c_str(), opts.condorOutputDir().c_str());
}
