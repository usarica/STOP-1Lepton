#include "common_includes.h"


class EventAnalyzer : public FrameworkTreeLooperBase{
protected:

  bool runEvent(FrameworkTree* tree, float const& externalWgt, SimpleEntry& product);

public:
  EventAnalyzer();
  EventAnalyzer(FrameworkTree* inTree);
  EventAnalyzer(std::vector<FrameworkTree*> const& inTreeList);
  EventAnalyzer(FrameworkSet const* inTreeSet);

};

EventAnalyzer::EventAnalyzer() : FrameworkTreeLooperBase() {}
EventAnalyzer::EventAnalyzer(FrameworkTree* inTree) : FrameworkTreeLooperBase(inTree) {}
EventAnalyzer::EventAnalyzer(std::vector<FrameworkTree*> const& inTreeList) : FrameworkTreeLooperBase(inTreeList) {}
EventAnalyzer::EventAnalyzer(FrameworkSet const* inTreeSet) : FrameworkTreeLooperBase(inTreeSet) {}

bool EventAnalyzer::runEvent(FrameworkTree* tree, float const& externalWgt, SimpleEntry& product){
  constexpr bool doWeights=false;
  constexpr bool doElectrons=true;

  bool validProducts = (tree!=nullptr);
  if (!validProducts) return validProducts;

  WeightsHandler* wgtHandler=nullptr;
  ElectronHandler* eleHandler=nullptr;
  for (auto it_ivy=this->externalIvyObjects.begin(); it_ivy!=this->externalIvyObjects.end(); it_ivy++){
    if (doWeights){ WeightsHandler* wgt_ivy = dynamic_cast<WeightsHandler*>(it_ivy->second); if (!wgtHandler && wgt_ivy){ wgtHandler = wgt_ivy; } }
    if (doElectrons){ ElectronHandler* ele_ivy = dynamic_cast<ElectronHandler*>(it_ivy->second); if (!eleHandler && ele_ivy){ eleHandler = ele_ivy; } }
  }

  /*********************/
  /**  WEIGHTS BLOCK  **/
  /*********************/
  float wgt = externalWgt;
  validProducts &= (!doWeights || wgtHandler!=nullptr || tree->isData());
  if (!validProducts){
    MELAerr << "EventAnalyzer::runEvent: Weight handle is invalid (Tree: " << tree->sampleIdentifier << ")." << endl;
    return validProducts;
  }
  if (wgtHandler){
    validProducts &= wgtHandler->constructWeights();
    if (!validProducts){
      MELAerr << "EventAnalyzer::runEvent: Weight product could not be constructed (Tree: " << tree->sampleIdentifier << ")." << endl;
      return validProducts;
    }
    validProducts &= wgtHandler->recordWeights(product, externalWgt);
    product.getNamedVal<float>(WeightVariables::getWeightName(WeightVariables::wCentral), wgt);
  }
  else product.setNamedVal<float>(WeightVariables::getWeightName(WeightVariables::wCentral), wgt);
  if (std::isnan(wgt) || std::isinf(wgt) || wgt==0.f){
    if (wgt!=0.f){
      MELAerr << "EventAnalyzer::runEvent: Invalid weight " << wgt << " is being discarded (Tree " << tree->sampleIdentifier << ")." << endl;
      exit(1);
    }
    validProducts=false;
  }
  if (!validProducts) return validProducts;


  /***********************/
  /**  ELECTRONS BLOCK  **/
  /***********************/
  validProducts &= (!doElectrons || eleHandler!=nullptr);
  if (!validProducts){
    MELAerr << "EventAnalyzer::runEvent: Electron handle is invalid (Tree: " << tree->sampleIdentifier << ")." << endl;
    return validProducts;
  }
  if (eleHandler){
    validProducts &= eleHandler->constructElectrons();
    if (!validProducts){
      MELAerr << "EventAnalyzer::runEvent: Electrons could not be constructed (Tree: " << tree->sampleIdentifier << ")." << endl;
      tree->print();
      exit(1);
      return validProducts;
    }
    std::vector<ElectronObject*> const& electrons = eleHandler->getProducts();
    std::vector<bool> conv_vtx_flag;
    std::vector<int> expectedMissingInnerHits;
    std::vector<float> energySC;
    std::vector<float> etaSC;
    std::vector<float> etaSeedSC;
    std::vector<float> rho;
    std::vector<float> sigmaIEtaIEta_full5x5;
    std::vector<float> dEtaIn;
    std::vector<float> dPhiIn;
    std::vector<float> hOverE;
    std::vector<float> ecalEnergy;
    std::vector<float> eOverPIn;
    std::vector<float> dxyPV;
    std::vector<float> dzPV;
    std::vector<float> miniIso_ch;
    std::vector<float> miniIso_nh;
    std::vector<float> miniIso_em;
    std::vector<bool> isVeto;
    std::vector<bool> isLoose;
    std::vector<bool> isMedium;
    for (ElectronObject const* electron:electrons){
      if (!electron) continue;
      ElectronVariables const& extras = electron->extras;

      conv_vtx_flag.push_back(extras.conv_vtx_flag);
      expectedMissingInnerHits.push_back(extras.expectedMissingInnerHits);
      energySC.push_back(extras.energySC);
      etaSC.push_back(extras.etaSC);
      etaSeedSC.push_back(extras.etaSeedSC);
      rho.push_back(extras.rho);
      sigmaIEtaIEta_full5x5.push_back(extras.sigmaIEtaIEta_full5x5);
      dEtaIn.push_back(extras.dEtaIn);
      dPhiIn.push_back(extras.dPhiIn);
      hOverE.push_back(extras.hOverE);
      ecalEnergy.push_back(extras.ecalEnergy);
      eOverPIn.push_back(extras.eOverPIn);
      dxyPV.push_back(extras.dxyPV);
      dzPV.push_back(extras.dzPV);
      miniIso_ch.push_back(extras.miniIso_ch);
      miniIso_nh.push_back(extras.miniIso_nh);
      miniIso_em.push_back(extras.miniIso_em);

      isVeto.push_back(ElectronSelectionHelpers::testVetoSelection(*electron));
      isLoose.push_back(ElectronSelectionHelpers::testLooseSelection(*electron));
      isMedium.push_back(ElectronSelectionHelpers::testMediumSelection(*electron));
    }
    product.setNamedVal<std::vector<bool>>("electrons_conv_vtx_flag", conv_vtx_flag);
    product.setNamedVal<std::vector<int>>("electrons_expectedMissingInnerHits", expectedMissingInnerHits);
    product.setNamedVal<std::vector<float>>("electrons_energySC", energySC);
    product.setNamedVal<std::vector<float>>("electrons_etaSC", etaSC);
    product.setNamedVal<std::vector<float>>("electrons_etaSeedSC", etaSeedSC);
    product.setNamedVal<std::vector<float>>("electrons_rho", rho);
    product.setNamedVal<std::vector<float>>("electrons_sigmaIEtaIEta_full5x5", sigmaIEtaIEta_full5x5);
    product.setNamedVal<std::vector<float>>("electrons_dEtaIn", dEtaIn);
    product.setNamedVal<std::vector<float>>("electrons_dPhiIn", dPhiIn);
    product.setNamedVal<std::vector<float>>("electrons_hOverE", hOverE);
    product.setNamedVal<std::vector<float>>("electrons_ecalEnergy", ecalEnergy);
    product.setNamedVal<std::vector<float>>("electrons_eOverPIn", eOverPIn);
    product.setNamedVal<std::vector<float>>("electrons_dxyPV", dxyPV);
    product.setNamedVal<std::vector<float>>("electrons_dzPV", dzPV);
    product.setNamedVal<std::vector<float>>("electrons_miniIso_ch", miniIso_ch);
    product.setNamedVal<std::vector<float>>("electrons_miniIso_nh", miniIso_nh);
    product.setNamedVal<std::vector<float>>("electrons_miniIso_em", miniIso_em);

    product.setNamedVal<std::vector<bool>>("electrons_isVeto", isVeto);
    product.setNamedVal<std::vector<bool>>("electrons_isLoose", isLoose);
    product.setNamedVal<std::vector<bool>>("electrons_isMedium", isMedium);
  }

  return validProducts;
}


void testWeights(){
  TDirectory* curdir = gDirectory;

  std::string stropts = "indir=/hadoop/cms/store/group/snt/run2_mc2018 outdir=./ outfile=WZZ_TuneCP5_13TeV-amcatnlo-pythia8.root sample=/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM year=2018 maxevents=100 ismc=true";
  FrameworkOptionParser opts(stropts);

  BaseTree outtree("test"); // The tree to record into the ROOT file
  TFile* foutput = TFile::Open((opts.outputDir()+opts.outputFilename()).c_str(), "recreate");
  curdir->cd();

  FrameworkSet theSet(opts, CMS4_EVENTS_TREE_NAME);

  WeightsHandler wgtHandler;
  wgtHandler.setVerbosity(TVar::DEBUG);
  for (auto* tree:theSet.getFrameworkTreeList()) wgtHandler.bookBranches(tree);

  ElectronHandler eleHandler;
  eleHandler.setVerbosity(TVar::DEBUG);
  for (auto* tree:theSet.getFrameworkTreeList()) eleHandler.bookBranches(tree);

  EventAnalyzer analyzer(&theSet);
  analyzer.setMaximumEvents(opts.maxEventsToProcess());
  analyzer.addExternalIvyObject("WeightsHandler", &wgtHandler);
  analyzer.addExternalIvyObject("ElectronHandler", &eleHandler);
  analyzer.setExternalProductTree(&outtree);

  analyzer.loop(true, false, true);
  MELAout << "There are " << outtree.getNEvents() << " products" << endl;
  outtree.writeToFile(foutput);

  foutput->Close();
}
