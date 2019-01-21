#include <cassert>
#include "Samples.h"
#include "ParticleObjectHelpers.h"
#include "JetMETHandler.h"
#include "ElectronSelectionHelpers.h"
#include "MuonSelectionHelpers.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "JECJERHelpers.h"
#include "TFTopTaggerHelpers.h"
#include "FrameworkVariables.hh"
#include "FrameworkTree.h"
#include "MELAStreamHelpers.hh"
#include "TRandom3.h"
#include <cmstas/CORE/Tools/jetcorr/JetCorrectionUncertainty.h>
#include <cmstas/CORE/Tools/jetcorr/JetCorrectorParameters.h>


using namespace std;
using namespace MELAStreamHelpers;


JetMETHandler::JetMETHandler() :
  IvyBase(),

  doGenJets(true),
  doAK4Jets(true),
  doAK8Jets(true),
  doMET(true),
  doTops(true),

  metobj(nullptr),
  registeredBtagSFHandler(nullptr),
  registeredBtagSFHandler_FastSim(nullptr),
  registeredJECSFHandler_ak4jets(nullptr),
  registeredJECSFHandler_ak8jets(nullptr),
  registeredJERSFHandler_ak4jets(nullptr),
  registeredJERSFHandler_ak8jets(nullptr),
  registeredElectrons(nullptr),
  registeredMuons(nullptr)
{}

void JetMETHandler::clear(){
  // Clear everything in reverse order of creation
  delete metobj; metobj=nullptr;

  for (auto& obj:tftops) delete obj;
  tftops.clear();

  for (auto& obj:ak8jets) delete obj;
  ak8jets.clear();

  for (auto& obj:ak4jets) delete obj;
  ak4jets.clear();

  for (auto& obj:genjets) delete obj;
  genjets.clear();

  // Do not clear registered objects here
}


bool JetMETHandler::constructGenJets(){
  if (!doGenJets) return true;

  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree) return false;
  if (!(fwktree->isMC() && fwktree->branchExists(_genjets_momentum_))) return true;

  std::vector<CMSLorentzVector>* genjets_momentum = nullptr;
  if (!this->getConsumedValue(_genjets_momentum_, genjets_momentum)) return false;
  if (!genjets_momentum) return false;

  // Yes, this is the only thing to do...
  for (auto const& mom:(*genjets_momentum)) genjets.push_back(new GenJetObject(mom));

  return true;
}

bool JetMETHandler::constructAK4Jets(){
  if (!doAK4Jets) return true;

  float rho = 0;

  std::vector<int>* npfcands = nullptr;
  std::vector<int>* parton_flavor = nullptr;
  std::vector<int>* hadron_flavor = nullptr;
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
    && this->getConsumedValue(_ak4jets_parton_flavor_, parton_flavor)
    && this->getConsumedValue(_ak4jets_hadron_flavor_, hadron_flavor)
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

  if (momentum->empty()) return true; // Construction is successful, it is just that no jets exist.

  unsigned int nProducts = momentum->size();
  ak4jets.reserve(nProducts);
  for (unsigned int ip=0; ip<nProducts; ip++){
    if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK4Jets: Attempting jet " << ip << "..." << endl;

    ak4jets.push_back(new AK4JetObject(0, momentum->at(ip)));
    AK4JetObject*& obj = ak4jets.back();

    obj->extras.rho = rho;

    obj->extras.npfcands = npfcands->at(ip);
    obj->extras.parton_flavor = parton_flavor->at(ip);
    obj->extras.hadron_flavor = hadron_flavor->at(ip);
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

  return true;
}

bool JetMETHandler::constructAK8Jets(){
  if (!doAK8Jets) return true;

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

  if (momentum->empty()) return true; // Construction is successful, it is just that no jets exist.

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

  return true;
}

bool JetMETHandler::constructMET(){
  if (!doMET) return true;

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
  if (!doTops) return true;

  if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructTFTops: Attempting to construct tops from ak4jets!" << endl;
  tftops = TFTopTaggerHelpers::getTopsFromResolvedJets(ak4jets);
  if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success! Constructed " << tftops.size() << " tops." << endl;
  return true;
}

bool JetMETHandler::constructJetMET(){
  clear();
  if (!currentTree) return false;
  bool res = (constructGenJets() && constructAK4Jets() && constructAK8Jets() && constructMET());

  // Apply selections after sorting
  res &= applyJEC(); // Re-apply JEC before calculating b-tag SFs
  res &= applyBtagSFs(); // Also handles b-tagging itself
  res &= matchGenJets(); // Match before JER
  res &= applyJER(); // Apply JER AFTER calculating b-tag SFs
  res &= applyJetCleaning();

  // Sort particles after jet cleaning is done
  ParticleObjectHelpers::sortByGreaterPt(genjets);
  ParticleObjectHelpers::sortByGreaterPt(ak4jets);
  ParticleObjectHelpers::sortByGreaterPt(ak8jets);

  res &= applySelections(); // Apply selections
  res &= constructTFTops(); // Construct top candidates at the latest stage

  // Sort the tops as well
  ParticleObjectHelpers::sortByGreaterPt(tftops);

  return res;
}

bool JetMETHandler::applyJetCleaning(){
  std::vector<AK4JetObject*> ak4jets_new; ak4jets_new.reserve(ak4jets.size());
  for (auto*& jet:ak4jets){
    bool doSkip=false;
    if (registeredMuons){
      for (auto const* part:*(registeredMuons)){
        if (!part->testSelection(MuonSelectionHelpers::kVetoIDReco) || !part->testSelection(MuonSelectionHelpers::kSkimPtEta)) continue;
        if (reco::deltaR(jet->getFinalMomentum(), part->momentum)<jet->ConeRadiusConstant){ doSkip=true; break; }
      }
    }
    if (registeredElectrons){
      for (auto const* part:*(registeredElectrons)){
        if (!part->testSelection(ElectronSelectionHelpers::kVetoIDReco) || !part->testSelection(ElectronSelectionHelpers::kSkimPtEta)) continue;
        if (reco::deltaR(jet->getFinalMomentum(), part->momentum)<jet->ConeRadiusConstant){ doSkip=true; break; }
      }
    }
    if (!doSkip) ak4jets_new.push_back(jet);
    else delete jet;
  }
  ak4jets = ak4jets_new;

  std::vector<AK8JetObject*> ak8jets_new; ak8jets_new.reserve(ak8jets.size());
  for (auto*& jet:ak8jets){
    bool doSkip=false;
    if (registeredMuons){
      for (auto const* part:*(registeredMuons)){
        if (!part->testSelection(MuonSelectionHelpers::kVetoIDReco) || !part->testSelection(MuonSelectionHelpers::kSkimPtEta)) continue;
        if (reco::deltaR(jet->getFinalMomentum(), part->momentum)<jet->ConeRadiusConstant){ doSkip=true; break; }
      }
    }
    if (registeredElectrons){
      for (auto const* part:*(registeredElectrons)){
        if (!part->testSelection(ElectronSelectionHelpers::kVetoIDReco) || !part->testSelection(ElectronSelectionHelpers::kSkimPtEta)) continue;
        if (reco::deltaR(jet->getFinalMomentum(), part->momentum)<jet->ConeRadiusConstant){ doSkip=true; break; }
      }
    }
    if (!doSkip) ak8jets_new.push_back(jet);
    else delete jet;
  }
  ak8jets = ak8jets_new;

  registeredMuons=nullptr; // De-register muons now
  registeredElectrons=nullptr; // De-register electrons now
  return true;
}
bool JetMETHandler::applyJEC(){
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree) return false;

  bool const doCorrectAK4Jets = true;
  if (registeredJECSFHandler_ak4jets && doCorrectAK4Jets){
    for (AK4JetObject* jet:ak4jets) registeredJECSFHandler_ak4jets->applyJEC(jet, fwktree->isMC(), fwktree->isFastSim());
  }

  bool const doCorrectAK8Jets = (SampleHelpers::theDataYear == 2016); // FIXME: Need the prescription for 2018 as well
  if (registeredJECSFHandler_ak8jets && doCorrectAK8Jets){
    for (AK8JetObject* jet:ak8jets) registeredJECSFHandler_ak8jets->applyJEC(jet, fwktree->isMC(), fwktree->isFastSim());
  }

  return true;
}
bool JetMETHandler::matchGenJets(){
  if (genjets.empty()) return true;
  // Determine the pairing absolutely

  // AK4 jets
  if (!ak4jets.empty()){
    std::vector<unsigned int> remaining_recojets; if (ak4jets.size()>0) remaining_recojets.reserve(ak4jets.size());
    for (unsigned int ijet=0; ijet<ak4jets.size(); ijet++) remaining_recojets.push_back(ijet);
    std::vector<unsigned int> remaining_genjets; remaining_genjets.reserve(genjets.size());
    for (unsigned int ijet=0; ijet<genjets.size(); ijet++) remaining_genjets.push_back(ijet);
    while (!remaining_recojets.empty() && !remaining_genjets.empty()){
      int chosenRecoJet=-1;
      int chosenGenJet=-1;
      float minDeltaR=-1;
      for (unsigned int const& rjet:remaining_recojets){
        CMSLorentzVector pReco = ak4jets.at(rjet)->getFinalMomentum();
        for (unsigned int const& gjet:remaining_genjets){
          CMSLorentzVector const& pGen = genjets.at(gjet)->momentum;
          float deltaR = reco::deltaR(pGen, pReco);
          if (minDeltaR==-1. || deltaR<minDeltaR){
            minDeltaR=deltaR;
            chosenRecoJet=rjet;
            chosenGenJet=gjet;
          }
        }
      }

      if (chosenRecoJet>=0 && chosenGenJet>=0) ak4jets.at(chosenRecoJet)->associatedGenJet = genjets.at(chosenGenJet);
      for (auto it=remaining_recojets.begin(); it!=remaining_recojets.end(); it++){ if ((int) *it == chosenRecoJet){ remaining_recojets.erase(it); break; } }
      for (auto it=remaining_genjets.begin(); it!=remaining_genjets.end(); it++){ if ((int) *it == chosenGenJet){ remaining_genjets.erase(it); break; } }
    }
  }

  // AK8 jets
  if (!ak8jets.empty()){
    std::vector<unsigned int> remaining_recojets; if (ak8jets.size()>0) remaining_recojets.reserve(ak8jets.size());
    for (unsigned int ijet=0; ijet<ak8jets.size(); ijet++) remaining_recojets.push_back(ijet);
    std::vector<unsigned int> remaining_genjets; remaining_genjets.reserve(genjets.size());
    for (unsigned int ijet=0; ijet<genjets.size(); ijet++) remaining_genjets.push_back(ijet);
    while (!remaining_recojets.empty() && !remaining_genjets.empty()){
      int chosenRecoJet=-1;
      int chosenGenJet=-1;
      float minDeltaR=-1;
      for (unsigned int const& rjet:remaining_recojets){
        CMSLorentzVector pReco = ak8jets.at(rjet)->getFinalMomentum();
        for (unsigned int const& gjet:remaining_genjets){
          CMSLorentzVector const& pGen = genjets.at(gjet)->momentum;
          float deltaR = reco::deltaR(pGen, pReco);
          if (minDeltaR==-1. || deltaR<minDeltaR){
            minDeltaR=deltaR;
            chosenRecoJet=rjet;
            chosenGenJet=gjet;
          }
        }
      }

      if (chosenRecoJet>=0 && chosenGenJet>=0) ak8jets.at(chosenRecoJet)->associatedGenJet = genjets.at(chosenGenJet);
      for (auto it=remaining_recojets.begin(); it!=remaining_recojets.end(); it++){ if ((int) *it == chosenRecoJet){ remaining_recojets.erase(it); break; } }
      for (auto it=remaining_genjets.begin(); it!=remaining_genjets.end(); it++){ if ((int) *it == chosenGenJet){ remaining_genjets.erase(it); break; } }
    }
  }

  return true;
}
bool JetMETHandler::applyJER(){
  if (!doAK4Jets && !doAK8Jets) return true;

  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree) return false;
  if (!(fwktree->isMC() && fwktree->branchExists(_genjets_momentum_))) return true;
  std::vector<CMSLorentzVector>* genjets_momentum = nullptr;
  if (!this->getConsumedValue(_genjets_momentum_, genjets_momentum)) return false;
  if (!genjets_momentum) return false;

  if (registeredJERSFHandler_ak4jets && !ak4jets.empty()){
    float rho=0;
    for (auto*& jet:ak4jets){
      AK4JetVariables const& extras = jet->extras;
      rho = extras.rho;
      break;
    }
    registeredJERSFHandler_ak4jets->smear(ak4jets, rho);
  }
  if (registeredJERSFHandler_ak8jets && !ak8jets.empty()){
    float rho = 0;
    for (auto*& jet:ak8jets){
      AK8JetVariables const& extras = jet->extras;
      rho = extras.rho;
      break;
    }
    registeredJERSFHandler_ak8jets->smear(ak8jets, rho);
  }

  return true;
}
bool JetMETHandler::applySelections(){
  for (auto* obj:ak4jets){ /*obj->resetSelectionBits();*/ AK4JetSelectionHelpers::setSelectionBits(*obj); }
  for (auto* obj:ak8jets){ /*obj->resetSelectionBits();*/ AK8JetSelectionHelpers::setSelectionBits(*obj); }
  return true;
}
bool JetMETHandler::applyBtagSFs(){
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree) return false;

  BtagScaleFactorHandler* theBTagSFHandler = nullptr;
  if (fwktree->isMC() && fwktree->isFastSim()) theBTagSFHandler = registeredBtagSFHandler_FastSim;
  else if (fwktree->isMC()) theBTagSFHandler = registeredBtagSFHandler;
  if (!theBTagSFHandler) return true;

  float bTaggerThreshold;
  if (theBTagSFHandler) bTaggerThreshold = theBTagSFHandler->WPval;
  else bTaggerThreshold = BtagHelpers::getBtagWP(AK4JetSelectionHelpers::AK4Jets_BTagWPType);

  for (auto*& jet:ak4jets){
    AK4JetVariables const& extras = jet->extras;

    float bTagger = extras.deepCSVb + extras.deepCSVbb;
    bool isBtagged = (bTagger > bTaggerThreshold);

    bool isBtaggedWithSF   = isBtagged;
    bool isBtaggedWithSFUp = isBtagged;
    bool isBtaggedWithSFDn = isBtagged;
    if (!theBTagSFHandler){
      int const& flav = extras.hadron_flavor;
      CMSLorentzVector correctedMomentum = jet->getFinalMomentum();
      float jpt = correctedMomentum.Pt();
      float jphi = correctedMomentum.Phi();
      float jeta = correctedMomentum.Eta();
      TRandom3 rand;
      rand.SetSeed(std::abs(static_cast<int>(sin(jphi)*100000)));
      float R = rand.Uniform();
      float SF   = theBTagSFHandler->getSF(0, flav, jpt, jeta);
      float SFUp = theBTagSFHandler->getSF(1, flav, jpt, jeta);
      float SFDn = theBTagSFHandler->getSF(-1, flav, jpt, jeta);
      float bTagMCEff = theBTagSFHandler->getEff(flav, jpt, jeta);
      if (SF  <=1.f && isBtagged && R<1.f-SF) isBtaggedWithSF   = false;
      if (SFUp<=1.f && isBtagged && R<1.f-SFUp) isBtaggedWithSFUp = false;
      if (SFDn<=1.f && isBtagged && R<1.f-SFDn) isBtaggedWithSFDn = false;
      if (SF  >1.f && !isBtagged && R<(1.f-SF)/(1.f-1.f/bTagMCEff)) isBtaggedWithSF   = true;
      if (SFUp>1.f && !isBtagged && R<(1.f-SFUp)/(1.f-1.f/bTagMCEff)) isBtaggedWithSFUp = true;
      if (SFDn>1.f && !isBtagged && R<(1.f-SFDn)/(1.f-1.f/bTagMCEff)) isBtaggedWithSFDn = true;
    }
    if (isBtaggedWithSF) jet->setSelectionBit(AK4JetSelectionHelpers::kIsBTagged);
    if (isBtaggedWithSFUp) jet->setSelectionBit(AK4JetSelectionHelpers::kIsBTagged_SFUp);
    if (isBtaggedWithSFDn) jet->setSelectionBit(AK4JetSelectionHelpers::kIsBTagged_SFDn);
  }

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
  SampleHelpers::setupUsingOptions(fwktree->getOptions());

  // ak4 jet variables
  if (doAK4Jets){
    this->addConsumed<float>(_ak4jets_rho_);

    this->addConsumed<std::vector<int>*>(_ak4jets_npfcands_);
    this->addConsumed<std::vector<int>*>(_ak4jets_parton_flavor_);
    this->addConsumed<std::vector<int>*>(_ak4jets_hadron_flavor_);
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

    fwktree->bookEDMBranch<float>(_ak4jets_rho_, 0);

    fwktree->bookEDMBranch<std::vector<int>*>(_ak4jets_npfcands_, nullptr);
    fwktree->bookEDMBranch<std::vector<int>*>(_ak4jets_parton_flavor_, nullptr);
    fwktree->bookEDMBranch<std::vector<int>*>(_ak4jets_hadron_flavor_, nullptr);
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
  }

  // ak8 jet variables
  if (doAK8Jets){
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
  }

  // Gen. jet variables
  if (doGenJets){
    this->addConsumed<std::vector<CMSLorentzVector>*>(_genjets_momentum_);
    this->defineConsumedSloppy(_genjets_momentum_); // Define this sloppy for data might not have it
    if (fwktree->isMC()){
      bool hasGenJets = (SampleHelpers::branchExists(fwktree->getSelectedTree(), _genjets_momentum_) || SampleHelpers::aliasExists(fwktree->getSelectedTree(), _genjets_momentum_));
      if (hasGenJets) fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_genjets_momentum_, nullptr);
    }
  }

  // MET variables
  if (doMET){
    this->addConsumed<float>(_pfmet_);
    this->addConsumed<float>(_pfmetPhi_);
    fwktree->bookEDMBranch<float>(_pfmet_, 0);
    fwktree->bookEDMBranch<float>(_pfmetPhi_, 0);
  }
}
