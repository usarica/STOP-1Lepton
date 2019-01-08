#include <cassert>
#include "ParticleObjectHelpers.h"
#include "JetMETHandler.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "TFTopTaggerHelpers.h"
#include "FrameworkVariables.hh"
#include "FrameworkTree.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


JetMETHandler::JetMETHandler() :
  IvyBase(),
  metobj(nullptr),
  registeredElectrons(nullptr),
  registeredMuons(nullptr)
{
  // ak4 jet variables
  this->addConsumed<float>(_ak4jets_rho_);

  this->addConsumed<std::vector<int>*>(_ak4jets_npfcands_);
  this->addConsumed<std::vector<int>*>(_ak4jets_chargedHadronMultiplicity_);
  this->addConsumed<std::vector<int>*>(_ak4jets_neutralHadronMultiplicity_);
  this->addConsumed<std::vector<int>*>(_ak4jets_photonMultiplicity_);
  this->addConsumed<std::vector<int>*>(_ak4jets_electronMultiplicity_);
  this->addConsumed<std::vector<int>*>(_ak4jets_muonMultiplicity_);
  this->addConsumed<std::vector<int>*>(_ak4jets_chargedMultiplicity_);
  this->addConsumed<std::vector<int>*>(_ak4jets_neutralMultiplicity_);
  this->addConsumed<std::vector<int>*>(_ak4jets_totalMultiplicity_);

  this->addConsumed<std::vector<float>*>(_ak4jets_area_);
  this->addConsumed<std::vector<float>*>(_ak4jets_undoJEC_);
  this->addConsumed<std::vector<float>*>(_ak4jets_chargedHadronE_);
  this->addConsumed<std::vector<float>*>(_ak4jets_chargedEmE_);
  this->addConsumed<std::vector<float>*>(_ak4jets_neutralHadronE_);
  this->addConsumed<std::vector<float>*>(_ak4jets_neutralEmE_);
  this->addConsumed<std::vector<float>*>(_ak4jets_hfHadronE_);
  this->addConsumed<std::vector<float>*>(_ak4jets_hfEmE_);
  this->addConsumed<std::vector<float>*>(_ak4jets_photonE_);
  this->addConsumed<std::vector<float>*>(_ak4jets_electronE_);
  this->addConsumed<std::vector<float>*>(_ak4jets_muonE_);

  this->addConsumed<std::vector<float>*>(_ak4jets_pfCombinedInclusiveSecondaryVertexV2BJetTag_);
  this->addConsumed<std::vector<float>*>(_ak4jets_ptDistribution_);
  this->addConsumed<std::vector<float>*>(_ak4jets_axis1_);
  this->addConsumed<std::vector<float>*>(_ak4jets_axis2_);

  this->addConsumed<std::vector<TString>*>(_ak4jets_bDiscriminatorNames);

  this->addConsumed<std::vector<std::vector<float>>*>(_ak4jets_bDiscriminators);

  this->addConsumed<std::vector<CMSLorentzVector>*>(_ak4jets_momentum_);


  // ak8 jet variables
  this->addConsumed<float>(_ak8jets_rho_);

  this->addConsumed<std::vector<int>*>(_ak8jets_parton_flavor_);

  this->addConsumed<std::vector<float>*>(_ak8jets_area_);
  this->addConsumed<std::vector<float>*>(_ak8jets_undoJEC_);
  this->addConsumed<std::vector<float>*>(_ak8jets_tau1_);
  this->addConsumed<std::vector<float>*>(_ak8jets_tau2_);
  this->addConsumed<std::vector<float>*>(_ak8jets_tau3_);
  this->addConsumed<std::vector<float>*>(_ak8jets_deepdisc_qcd_);
  this->addConsumed<std::vector<float>*>(_ak8jets_deepdisc_top_);
  this->addConsumed<std::vector<float>*>(_ak8jets_deepdisc_w_);
  this->addConsumed<std::vector<float>*>(_ak8jets_deepdisc_z_);
  this->addConsumed<std::vector<float>*>(_ak8jets_deepdisc_zbb_);
  this->addConsumed<std::vector<float>*>(_ak8jets_deepdisc_hbb_);
  this->addConsumed<std::vector<float>*>(_ak8jets_deepdisc_h4q_);

  this->addConsumed<std::vector<CMSLorentzVector>*>(_ak8jets_momentum_);


  // MET variables
  this->addConsumed<float>(_pfmet_);
  this->addConsumed<float>(_pfmetPhi_);
}

void JetMETHandler::clear(){
  for (auto& obj:ak4jets) delete obj;
  ak4jets.clear();
  for (auto& obj:ak8jets) delete obj;
  ak8jets.clear();
  for (auto& obj:tftops) delete obj;
  tftops.clear();
  delete metobj; metobj=nullptr;
}


bool JetMETHandler::constructAK4Jets(){
  float rho = 0;

  std::vector<int>* npfcands = nullptr;
  std::vector<int>* chargedHadronMultiplicity = nullptr;
  std::vector<int>* neutralHadronMultiplicity = nullptr;
  std::vector<int>* photonMultiplicity = nullptr;
  std::vector<int>* electronMultiplicity = nullptr;
  std::vector<int>* muonMultiplicity = nullptr;
  std::vector<int>* chargedMultiplicity = nullptr;
  std::vector<int>* neutralMultiplicity = nullptr;
  std::vector<int>* totalMultiplicity = nullptr;

  std::vector<float>* area = nullptr;
  std::vector<float>* undoJEC = nullptr;
  std::vector<float>* chargedHadronE = nullptr;
  std::vector<float>* chargedEmE = nullptr;
  std::vector<float>* neutralHadronE = nullptr;
  std::vector<float>* neutralEmE = nullptr;
  std::vector<float>* hfHadronE = nullptr;
  std::vector<float>* hfEmE = nullptr;
  std::vector<float>* photonE = nullptr;
  std::vector<float>* electronE = nullptr;
  std::vector<float>* muonE = nullptr;

  std::vector<float>* pfCombinedInclusiveSecondaryVertexV2BJetTag = nullptr;
  std::vector<float>* ptDistribution = nullptr;
  std::vector<float>* axis1 = nullptr;
  std::vector<float>* axis2 = nullptr;

  std::vector<TString>* bDiscriminatorNames = nullptr;

  std::vector<std::vector<float>>* bDiscriminators = nullptr;

  std::vector<CMSLorentzVector>* momentum = nullptr;

  // Beyond this point starts checks and selection
  bool allVariablesPresent = (
    this->getConsumedValue(_ak4jets_rho_, rho)

    && this->getConsumedValue(_ak4jets_npfcands_, npfcands)
    && this->getConsumedValue(_ak4jets_chargedHadronMultiplicity_, chargedHadronMultiplicity)
    && this->getConsumedValue(_ak4jets_neutralHadronMultiplicity_, neutralHadronMultiplicity)
    && this->getConsumedValue(_ak4jets_photonMultiplicity_, photonMultiplicity)
    && this->getConsumedValue(_ak4jets_electronMultiplicity_, electronMultiplicity)
    && this->getConsumedValue(_ak4jets_muonMultiplicity_, muonMultiplicity)
    && this->getConsumedValue(_ak4jets_chargedMultiplicity_, chargedMultiplicity)
    && this->getConsumedValue(_ak4jets_neutralMultiplicity_, neutralMultiplicity)
    && this->getConsumedValue(_ak4jets_totalMultiplicity_, totalMultiplicity)

    && this->getConsumedValue(_ak4jets_area_, area)
    && this->getConsumedValue(_ak4jets_undoJEC_, undoJEC)
    && this->getConsumedValue(_ak4jets_chargedHadronE_, chargedHadronE)
    && this->getConsumedValue(_ak4jets_chargedEmE_, chargedEmE)
    && this->getConsumedValue(_ak4jets_neutralHadronE_, neutralHadronE)
    && this->getConsumedValue(_ak4jets_neutralEmE_, neutralEmE)
    && this->getConsumedValue(_ak4jets_hfHadronE_, hfHadronE)
    && this->getConsumedValue(_ak4jets_hfEmE_, hfEmE)
    && this->getConsumedValue(_ak4jets_photonE_, photonE)
    && this->getConsumedValue(_ak4jets_electronE_, electronE)
    && this->getConsumedValue(_ak4jets_muonE_, muonE)

    && this->getConsumedValue(_ak4jets_pfCombinedInclusiveSecondaryVertexV2BJetTag_, pfCombinedInclusiveSecondaryVertexV2BJetTag)
    && this->getConsumedValue(_ak4jets_ptDistribution_, ptDistribution)
    && this->getConsumedValue(_ak4jets_axis1_, axis1)
    && this->getConsumedValue(_ak4jets_axis2_, axis2)

    && this->getConsumedValue(_ak4jets_bDiscriminatorNames, bDiscriminatorNames)

    && this->getConsumedValue(_ak4jets_bDiscriminators, bDiscriminators)

    && this->getConsumedValue(_ak4jets_momentum_, momentum)
    );

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "JetMETHandler::constructAK4Jets: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK4Jets: All variables are set up!" << endl;
  if (!momentum){
    if (this->verbosity>=TVar::ERROR) MELAerr << "JetMETHandler::constructAK4Jets: Jets could not be linked! Pointer to " << _ak4jets_momentum_ << " is null!" << endl;
    return false;
  }

  if (momentum->empty()) return true; // Construction is successful, it is just that no muons exist.

  unsigned int nProducts = momentum->size();
  ak4jets.reserve(nProducts);
  for (unsigned int ip=0; ip<nProducts; ip++){
    if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK4Jets: Attempting jet " << ip << "..." << endl;

    ak4jets.push_back(new AK4JetObject(0, momentum->at(ip)));
    AK4JetObject*& obj = ak4jets.back();

    obj->extras.rho = rho;

    obj->extras.npfcands = npfcands->at(ip);
    obj->extras.chargedHadronMultiplicity = chargedHadronMultiplicity->at(ip);
    obj->extras.neutralHadronMultiplicity = neutralHadronMultiplicity->at(ip);
    obj->extras.photonMultiplicity = photonMultiplicity->at(ip);
    obj->extras.electronMultiplicity = electronMultiplicity->at(ip);
    obj->extras.muonMultiplicity = muonMultiplicity->at(ip);
    obj->extras.chargedMultiplicity = chargedMultiplicity->at(ip);
    obj->extras.neutralMultiplicity = neutralMultiplicity->at(ip);
    obj->extras.totalMultiplicity = totalMultiplicity->at(ip);

    obj->extras.area = area->at(ip);
    obj->extras.undoJEC = undoJEC->at(ip);
    obj->extras.chargedHadronE = chargedHadronE->at(ip);
    obj->extras.chargedEmE = chargedEmE->at(ip);
    obj->extras.neutralHadronE = neutralHadronE->at(ip);
    obj->extras.neutralEmE = neutralEmE->at(ip);
    obj->extras.hfHadronE = hfHadronE->at(ip);
    obj->extras.hfEmE = hfEmE->at(ip);
    obj->extras.photonE = photonE->at(ip);
    obj->extras.electronE = electronE->at(ip);
    obj->extras.muonE = muonE->at(ip);

    static const TString strDeepFlavorPrefix = JetMETHandler::getAK4JetDeepFlavorPrefix(*bDiscriminatorNames);

    obj->extras.deepCSVb = JetMETHandler::getBtagValueFromLists(*bDiscriminatorNames, *bDiscriminators, ip, strDeepFlavorPrefix+"JetTags:probb");
    obj->extras.deepCSVc = JetMETHandler::getBtagValueFromLists(*bDiscriminatorNames, *bDiscriminators, ip, strDeepFlavorPrefix+"JetTags:probc");
    obj->extras.deepCSVl = JetMETHandler::getBtagValueFromLists(*bDiscriminatorNames, *bDiscriminators, ip, strDeepFlavorPrefix+"JetTags:probudsg");
    obj->extras.deepCSVbb = JetMETHandler::getBtagValueFromLists(*bDiscriminatorNames, *bDiscriminators, ip, strDeepFlavorPrefix+"JetTags:probbb");
    obj->extras.deepCSVcc = -9000;
    obj->extras.pfCombinedInclusiveSecondaryVertexV2BJetTag = pfCombinedInclusiveSecondaryVertexV2BJetTag->at(ip);
    obj->extras.ptDistribution = ptDistribution->at(ip);
    obj->extras.axis1 = axis1->at(ip);
    obj->extras.axis2 = axis2->at(ip);

    // DO NOT SET THE SELECTION BITS YET!

    if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(ak4jets);

  return true;
}

bool JetMETHandler::constructAK8Jets(){
  float rho = 0;

  std::vector<int>* parton_flavor = nullptr;

  std::vector<float>* area = nullptr;
  std::vector<float>* undoJEC = nullptr;
  std::vector<float>* tau1 = nullptr;
  std::vector<float>* tau2 = nullptr;
  std::vector<float>* tau3 = nullptr;
  std::vector<float>* deepdisc_qcd = nullptr;
  std::vector<float>* deepdisc_top = nullptr;
  std::vector<float>* deepdisc_w = nullptr;
  std::vector<float>* deepdisc_z = nullptr;
  std::vector<float>* deepdisc_zbb = nullptr;
  std::vector<float>* deepdisc_hbb = nullptr;
  std::vector<float>* deepdisc_h4q = nullptr;

  std::vector<CMSLorentzVector>* momentum = nullptr;

  // Beyond this point starts checks and selection
  bool allVariablesPresent = (
    this->getConsumedValue(_ak8jets_rho_, rho)

    && this->getConsumedValue(_ak8jets_parton_flavor_, parton_flavor)

    && this->getConsumedValue(_ak8jets_area_, area)
    && this->getConsumedValue(_ak8jets_undoJEC_, undoJEC)
    && this->getConsumedValue(_ak8jets_tau1_, tau1)
    && this->getConsumedValue(_ak8jets_tau2_, tau2)
    && this->getConsumedValue(_ak8jets_tau3_, tau3)
    && this->getConsumedValue(_ak8jets_deepdisc_qcd_, deepdisc_qcd)
    && this->getConsumedValue(_ak8jets_deepdisc_top_, deepdisc_top)
    && this->getConsumedValue(_ak8jets_deepdisc_w_, deepdisc_w)
    && this->getConsumedValue(_ak8jets_deepdisc_z_, deepdisc_z)
    && this->getConsumedValue(_ak8jets_deepdisc_zbb_, deepdisc_zbb)
    && this->getConsumedValue(_ak8jets_deepdisc_hbb_, deepdisc_hbb)
    && this->getConsumedValue(_ak8jets_deepdisc_h4q_, deepdisc_h4q)

    && this->getConsumedValue(_ak8jets_momentum_, momentum)
    );

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "JetMETHandler::constructAK8Jets: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK8Jets: All variables are set up!" << endl;
  if (!momentum){
    if (this->verbosity>=TVar::ERROR) MELAerr << "JetMETHandler::constructAK8Jets: Jets could not be linked! Pointer to " << _ak8jets_momentum_ << " is null!" << endl;
    return false;
  }

  if (momentum->empty()) return true; // Construction is successful, it is just that no muons exist.

  unsigned int nProducts = momentum->size();
  ak8jets.reserve(nProducts);
  for (unsigned int ip=0; ip<nProducts; ip++){
    if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK8Jets: Attempting jet " << ip << "..." << endl;

    ak8jets.push_back(new AK8JetObject(0, momentum->at(ip)));
    AK8JetObject*& obj = ak8jets.back();

    obj->extras.rho = rho;

    obj->extras.parton_flavor = parton_flavor->at(ip);

    obj->extras.area = area->at(ip);
    obj->extras.undoJEC = undoJEC->at(ip);
    obj->extras.tau1 = tau1->at(ip);
    obj->extras.tau2 = tau2->at(ip);
    obj->extras.tau3 = tau3->at(ip);
    obj->extras.deepdisc_qcd = deepdisc_qcd->at(ip);
    obj->extras.deepdisc_top = deepdisc_top->at(ip);
    obj->extras.deepdisc_w = deepdisc_w->at(ip);
    obj->extras.deepdisc_z = deepdisc_z->at(ip);
    obj->extras.deepdisc_zbb = deepdisc_zbb->at(ip);
    obj->extras.deepdisc_hbb = deepdisc_hbb->at(ip);
    obj->extras.deepdisc_h4q = deepdisc_h4q->at(ip);

    // DO NOT SET THE SELECTION BITS YET!

    if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(ak8jets);

  return true;
}

bool JetMETHandler::constructMET(){
  float met = 0;
  float metPhi = 0;

  // Beyond this point starts checks and selection
  bool allVariablesPresent = (
    this->getConsumedValue(_pfmet_, met)
    && this->getConsumedValue(_pfmetPhi_, metPhi)
    );

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "JetMETHandler::constructMET: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructMET: All variables are set up!" << endl;

  metobj = new METObject;
  metobj->extras.met = met;
  metobj->extras.phi = metPhi;

  if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

  return true;
}

bool JetMETHandler::constructTFTops(){
  tftops = TFTopTaggerHelpers::getTopsFromResolvedJets(ak4jets);
  return true;
}

bool JetMETHandler::constructJetMET(){
  clear();
  if (!currentTree) return false;
  return (constructAK4Jets() && constructAK8Jets() && constructMET() && constructTFTops() && applyJetCleaning() && applySelections());
}
bool JetMETHandler::applyJetCleaning(){
  if (registeredMuons){

    registeredMuons=nullptr; // De-register muons now
  }
  if (registeredElectrons){

    registeredElectrons=nullptr; // De-register electrons now
  }
  return true;
}
bool JetMETHandler::applySelections(){
  for (auto* obj:ak4jets) AK4JetSelectionHelpers::setSelectionBits(*obj);
  for (auto* obj:ak8jets) AK8JetSelectionHelpers::setSelectionBits(*obj);
  return true;
}

TString JetMETHandler::getAK4JetDeepFlavorPrefix(std::vector<TString> const& bDiscriminatorNames){
  static TString deepCSV_prefix = "NULL";
  if (deepCSV_prefix == "NULL"){
    for (TString const& discName : bDiscriminatorNames){
      if (discName.Contains("pfDeepCSV")){ // 2017 convention
        deepCSV_prefix = "pfDeepCSV";
        break;
      }
      else if (discName.Contains("deepFlavour")){ // 2016 convention
        deepCSV_prefix = "deepFlavour";
        break;
      }
    }
    if (deepCSV_prefix == "NULL"){
      MELAerr << "JetMETHandler::getAK4JetDeepFlavorTag: Cannot find the DeepFlavor/DeepCSV discriminator name!" << endl;
      exit(1);
    }
  }
  return deepCSV_prefix;
}
float JetMETHandler::getBtagValueFromLists(std::vector<TString> const& bDiscriminatorNames, std::vector<std::vector<float>> const& btagvals, size_t ijet, TString btagname){
  size_t index;
  std::vector<TString>::const_iterator begin_it = bDiscriminatorNames.begin();
  std::vector<TString>::const_iterator end_it = bDiscriminatorNames.end();
  std::vector<TString>::const_iterator found_it = std::find(begin_it, end_it, btagname);
  if (found_it != end_it) index = found_it - begin_it;
  else{
    MELAerr << "JetMETHandler::getBtagValueFromLists: Cannot find b-discriminator " << btagname << endl;
    MELAerr << "\t- Available tag names: " << bDiscriminatorNames << endl;
    exit(1);
  }
  return btagvals.at(ijet).at(index);
}

void JetMETHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(tree);
  if (!fwktree) return;

  // ak4 jet variables
  fwktree->bookEDMBranch<float>(_ak4jets_rho_, 0);

  fwktree->bookEDMBranch<std::vector<int>*>(_ak4jets_npfcands_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_ak4jets_chargedHadronMultiplicity_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_ak4jets_neutralHadronMultiplicity_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_ak4jets_photonMultiplicity_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_ak4jets_electronMultiplicity_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_ak4jets_muonMultiplicity_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_ak4jets_chargedMultiplicity_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_ak4jets_neutralMultiplicity_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_ak4jets_totalMultiplicity_, nullptr);

  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_area_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_undoJEC_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_chargedHadronE_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_chargedEmE_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_neutralHadronE_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_neutralEmE_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_hfHadronE_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_hfEmE_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_photonE_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_electronE_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_muonE_, nullptr);

  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_pfCombinedInclusiveSecondaryVertexV2BJetTag_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_ptDistribution_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_axis1_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak4jets_axis2_, nullptr);

  fwktree->bookEDMBranch<std::vector<TString>*>(_ak4jets_bDiscriminatorNames, nullptr);

  fwktree->bookEDMBranch<std::vector<std::vector<float>>*>(_ak4jets_bDiscriminators, nullptr);

  fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_ak4jets_momentum_, nullptr);


  // ak8 jet variables
  if (std::string(_ak4jets_rho_)!=std::string(_ak8jets_rho_)) fwktree->bookEDMBranch<float>(_ak8jets_rho_, 0);

  fwktree->bookEDMBranch<std::vector<int>*>(_ak8jets_parton_flavor_, nullptr);

  fwktree->bookEDMBranch<std::vector<float>*>(_ak8jets_area_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak8jets_undoJEC_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak8jets_tau1_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak8jets_tau2_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak8jets_tau3_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak8jets_deepdisc_qcd_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak8jets_deepdisc_top_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak8jets_deepdisc_w_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak8jets_deepdisc_z_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak8jets_deepdisc_zbb_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak8jets_deepdisc_hbb_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_ak8jets_deepdisc_h4q_, nullptr);

  fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_ak8jets_momentum_, nullptr);


  // MET variables
  fwktree->bookEDMBranch<float>(_pfmet_, 0);
  fwktree->bookEDMBranch<float>(_pfmetPhi_, 0);
}
