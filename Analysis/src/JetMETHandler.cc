#include <cassert>
#include "Samples.h"
#include "ParticleObjectHelpers.h"
#include "JetMETHandler.h"
#include "MuonSelectionHelpers.h"
#include "ElectronSelectionHelpers.h"
#include "PhotonSelectionHelpers.h"
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


#define AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVES \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, npfcands) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, parton_flavor) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, hadron_flavor) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, chargedHadronMultiplicity) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, neutralHadronMultiplicity) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, photonMultiplicity) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, electronMultiplicity) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, muonMultiplicity) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, chargedMultiplicity) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, neutralMultiplicity) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, totalMultiplicity) \
\
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, area) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, undoJEC) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, chargedHadronE) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, chargedEmE) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, neutralHadronE) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, neutralEmE) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, hfHadronE) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, hfEmE) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, photonE) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, electronE) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, muonE) \
\
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, pfCombinedInclusiveSecondaryVertexV2BJetTag) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, ptDistribution) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, axis1) \
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, axis2) \
\
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<std::vector<float>>, bDiscriminators) \
\
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<CMSLorentzVector>, momentum) \
\
AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<std::vector<CMSLorentzVector>>, mucands_momentum)
// Note that (std::vector<TString>) bDiscriminatorNames is NOT part of this set of directive macros because its size is a fized size of names, not N_ak4jets


#define AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVES \
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, parton_flavor) \
\
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, area) \
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, undoJEC) \
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, tau1) \
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, tau2) \
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, tau3) \
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, deepdisc_qcd) \
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, deepdisc_top) \
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, deepdisc_w) \
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, deepdisc_z) \
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, deepdisc_zbb) \
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, deepdisc_hbb) \
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, deepdisc_h4q) \
\
AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<CMSLorentzVector>, momentum) \


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
  registeredMuons(nullptr),
  registeredElectrons(nullptr),
  registeredPhotons(nullptr)
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
  //if (!fwktree->isMC()) return true;
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

  bool allVariablesPresent = true;

  float rho = 0; allVariablesPresent &= this->getConsumedValue(_ak4jets_rho_, rho);
  std::vector<TString>* bDiscriminatorNames = nullptr; allVariablesPresent &= this->getConsumedValue(_ak4jets_bDiscriminatorNames_, bDiscriminatorNames);

#define AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(TYPE, NAME) \
TYPE::const_iterator itBegin_##NAME, itEnd_##NAME; \
allVariablesPresent &= this->getConsumedCIterators<TYPE>(_ak4jets_##NAME##_, &itBegin_##NAME, &itEnd_##NAME);
  AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "JetMETHandler::constructAK4Jets: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK4Jets: All variables are set up!" << endl;

  if (itBegin_momentum == itEnd_momentum) return true; // Construction is successful, it is just that no jets exist.

  size_t nProducts = (itEnd_momentum - itBegin_momentum);
  ak4jets.reserve(nProducts);
#define AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(TYPE, NAME) auto it_##NAME = itBegin_##NAME;
  AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE
#define AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(TYPE, NAME) it_##NAME++;
  {
    size_t ip=0;
    while (it_momentum != itEnd_momentum){
      if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK4Jets: Attempting jet " << ip << "..." << endl;

      ak4jets.push_back(new AK4JetObject(0, *it_momentum));
      AK4JetObject*& obj = ak4jets.back();

      obj->extras.rho = rho;

      obj->extras.npfcands = *it_npfcands;
      obj->extras.parton_flavor = *it_parton_flavor;
      obj->extras.hadron_flavor = *it_hadron_flavor;
      obj->extras.chargedHadronMultiplicity = *it_chargedHadronMultiplicity;
      obj->extras.neutralHadronMultiplicity = *it_neutralHadronMultiplicity;
      obj->extras.photonMultiplicity = *it_photonMultiplicity;
      obj->extras.electronMultiplicity = *it_electronMultiplicity;
      obj->extras.muonMultiplicity = *it_muonMultiplicity;
      obj->extras.chargedMultiplicity = *it_chargedMultiplicity;
      obj->extras.neutralMultiplicity = *it_neutralMultiplicity;
      obj->extras.totalMultiplicity = *it_totalMultiplicity;

      obj->extras.area = *it_area;
      obj->extras.undoJEC = *it_undoJEC;
      obj->extras.chargedHadronE = *it_chargedHadronE;
      obj->extras.chargedEmE = *it_chargedEmE;
      obj->extras.neutralHadronE = *it_neutralHadronE;
      obj->extras.neutralEmE = *it_neutralEmE;
      obj->extras.hfHadronE = *it_hfHadronE;
      obj->extras.hfEmE = *it_hfEmE;
      obj->extras.photonE = *it_photonE;
      obj->extras.electronE = *it_electronE;
      obj->extras.muonE = *it_muonE;

      obj->extras.momentum_nomus_uncor = obj->momentum*obj->extras.undoJEC;
      for (CMSLorentzVector const& mom:(*it_mucands_momentum)) obj->extras.momentum_nomus_uncor -= mom;

      static const TString strDeepFlavorPrefix = JetMETHandler::getAK4JetDeepFlavorPrefix(*bDiscriminatorNames);

      obj->extras.deepCSVb = JetMETHandler::getBtagValueFromLists(*bDiscriminatorNames, *it_bDiscriminators, strDeepFlavorPrefix+"JetTags:probb");
      obj->extras.deepCSVc = JetMETHandler::getBtagValueFromLists(*bDiscriminatorNames, *it_bDiscriminators, strDeepFlavorPrefix+"JetTags:probc");
      obj->extras.deepCSVl = JetMETHandler::getBtagValueFromLists(*bDiscriminatorNames, *it_bDiscriminators, strDeepFlavorPrefix+"JetTags:probudsg");
      obj->extras.deepCSVbb = JetMETHandler::getBtagValueFromLists(*bDiscriminatorNames, *it_bDiscriminators, strDeepFlavorPrefix+"JetTags:probbb");
      obj->extras.deepCSVcc = -9000;
      obj->extras.pfCombinedInclusiveSecondaryVertexV2BJetTag = *it_pfCombinedInclusiveSecondaryVertexV2BJetTag;
      obj->extras.ptDistribution = *it_ptDistribution;
      obj->extras.axis1 = *it_axis1;
      obj->extras.axis2 = *it_axis2;

      // DO NOT SET THE SELECTION BITS YET!

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
      AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVES
    }
  }
#undef AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE
  // Sorting is skipped intentionally here because JECs and JERs will be re-applied.

  return true;
}

bool JetMETHandler::constructAK8Jets(){
  if (!doAK8Jets) return true;

  bool allVariablesPresent = true;

  float rho = 0; allVariablesPresent &= this->getConsumedValue(_ak8jets_rho_, rho);

#define AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(TYPE, NAME) \
TYPE::const_iterator itBegin_##NAME, itEnd_##NAME; \
allVariablesPresent &= this->getConsumedCIterators<TYPE>(_ak8jets_##NAME##_, &itBegin_##NAME, &itEnd_##NAME);
  AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "JetMETHandler::constructAK8Jets: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK8Jets: All variables are set up!" << endl;

  if (itBegin_momentum == itEnd_momentum) return true; // Construction is successful, it is just that no jets exist.

  size_t nProducts = (itEnd_momentum - itBegin_momentum);
  ak8jets.reserve(nProducts);
#define AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(TYPE, NAME) auto it_##NAME = itBegin_##NAME;
  AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE
#define AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(TYPE, NAME) it_##NAME++;
  {
    size_t ip=0;
    while (it_momentum != itEnd_momentum){
      if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK8Jets: Attempting jet " << ip << "..." << endl;

      ak8jets.push_back(new AK8JetObject(0, *it_momentum));
      AK8JetObject*& obj = ak8jets.back();

      obj->extras.rho = rho;

      obj->extras.parton_flavor = *it_parton_flavor;

      obj->extras.area = *it_area;
      obj->extras.undoJEC = *it_undoJEC;
      obj->extras.tau1 = *it_tau1;
      obj->extras.tau2 = *it_tau2;
      obj->extras.tau3 = *it_tau3;
      obj->extras.deepdisc_qcd = *it_deepdisc_qcd;
      obj->extras.deepdisc_top = *it_deepdisc_top;
      obj->extras.deepdisc_w = *it_deepdisc_w;
      obj->extras.deepdisc_z = *it_deepdisc_z;
      obj->extras.deepdisc_zbb = *it_deepdisc_zbb;
      obj->extras.deepdisc_hbb = *it_deepdisc_hbb;
      obj->extras.deepdisc_h4q = *it_deepdisc_h4q;

      // DO NOT SET THE SELECTION BITS YET!

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
      AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVES
    }
  }
#undef AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE
  // Sorting is skipped intentionally here because JECs and JERs will be re-applied.

  return true;
}

bool JetMETHandler::constructMET(){
  if (!doMET) return true;

  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree) return false;

  float metraw = 0;
  float metrawPhi = 0;
  float metoriginal = 0;
  float metoriginalPhi = 0;

  // Beyond this point starts checks and selection
  bool allVariablesPresent = false;
  // Nowhere should we use "old" MET
  /*if (SampleHelpers::theDataYear == 2017 && SampleHelpers::theDataVersion == SampleHelpers::kCMSSW_9_4_X && fwktree->sampleIdentifier.Contains("09May2018")){
    // For 2017F data samples with this version, the EE noise issue is already fixed, so use the 'old' raw MET
    allVariablesPresent &= (
      this->getConsumedValue(_pfmetraw_old_, metraw)
      && this->getConsumedValue(_pfmetrawPhi_old_, metrawPhi)
      );
  }
  else */if (SampleHelpers::theDataYear == 2016 && SampleHelpers::theDataVersion == SampleHelpers::kCMSSW_8_0_X){
    // For 2016 data, MuEG cleaning needs to be done
    allVariablesPresent &= (
      this->getConsumedValue(_pfmetraw_muegclean_, metraw)
      && this->getConsumedValue(_pfmetrawPhi_muegclean_, metrawPhi)
      );
  }
  else{
    // For 2017 and 2018, this variables corresponds to having all fixes applied.
    // For 2016, this variable does not apply fixes, so the one above should be used instead.
    allVariablesPresent = (
      this->getConsumedValue(_pfmetraw_, metraw)
      && this->getConsumedValue(_pfmetrawPhi_, metrawPhi)
      );
  }
  allVariablesPresent &= (
    this->getConsumedValue(_pfmet_, metoriginal)
    && this->getConsumedValue(_pfmetPhi_, metoriginalPhi)
    );

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "JetMETHandler::constructMET: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructMET: All variables are set up!" << endl;

  metobj = new METObject;
  metobj->extras.met
    = metobj->extras.met_METup = metobj->extras.met_METdn
    = metobj->extras.met_JECup = metobj->extras.met_JECdn
    = metobj->extras.met_JERup = metobj->extras.met_JERdn
    = metobj->extras.met_PUup = metobj->extras.met_PUdn
    = metobj->extras.met_raw = metraw;
  metobj->extras.phi
    = metobj->extras.phi_METup = metobj->extras.phi_METdn
    = metobj->extras.phi_JECup = metobj->extras.phi_JECdn
    = metobj->extras.phi_JERup = metobj->extras.phi_JERdn
    = metobj->extras.phi_PUup = metobj->extras.phi_PUdn
    = metobj->extras.phi_raw = metrawPhi;
  metobj->extras.met_original = metoriginal;
  metobj->extras.phi_original = metoriginalPhi;

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
  res &= applyBtaggingAndSFs(); // Also handles b-tagging itself
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
    if (registeredPhotons){
      for (auto const* part:*(registeredPhotons)){
        if (!part->testSelection(PhotonSelectionHelpers::kLooseIDReco) || !part->testSelection(PhotonSelectionHelpers::kSkimPtEta)) continue;
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
    if (registeredPhotons){
      for (auto const* part:*(registeredPhotons)){
        if (!part->testSelection(PhotonSelectionHelpers::kLooseIDReco) || !part->testSelection(PhotonSelectionHelpers::kSkimPtEta)) continue;
        if (reco::deltaR(jet->getFinalMomentum(), part->momentum)<jet->ConeRadiusConstant){ doSkip=true; break; }
      }
    }
    if (!doSkip) ak8jets_new.push_back(jet);
    else delete jet;
  }
  ak8jets = ak8jets_new;

  registeredMuons=nullptr; // De-register muons now
  registeredElectrons=nullptr; // De-register electrons now
  registeredPhotons=nullptr; // De-register electrons now
  return true;
}
bool JetMETHandler::applyJEC(){
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree) return false;

  bool const use2017EENoiseRecipe = (SampleHelpers::theDataYear == 2017 && SampleHelpers::theDataVersion == SampleHelpers::kCMSSW_9_4_X && !fwktree->sampleIdentifier.Contains("09May2018"));
  //bool const useMuEGCleanedMETRecipe = (SampleHelpers::theDataYear == 2016 && SampleHelpers::theDataVersion == SampleHelpers::kCMSSW_8_0_X);

  constexpr bool doCorrectAK4Jets = true;
  if (registeredJECSFHandler_ak4jets && doCorrectAK4Jets){
    for (AK4JetObject* jet:ak4jets) registeredJECSFHandler_ak4jets->applyJEC(jet, fwktree->isMC(), fwktree->isFastSim());
    if (metobj){ // MET is corrected based on the ak4jets
      CMSLorentzVector v_met(metobj->extras.met_raw*cos(metobj->extras.phi_raw), metobj->extras.met_raw*sin(metobj->extras.phi_raw), 0, metobj->extras.met_raw);
      //CMSLorentzVector unshiftedSum(0, 0, 0, 0);
      CMSLorentzVector shiftedSum_nominal(0, 0, 0, 0);
      CMSLorentzVector shiftedSum_L1only(0, 0, 0, 0);
      CMSLorentzVector shiftedSum_JECup(0, 0, 0, 0);
      CMSLorentzVector shiftedSum_JECdn(0, 0, 0, 0);
      for (AK4JetObject* jet:ak4jets){
        if ((jet->extras.chargedEmE + jet->extras.neutralEmE)/(jet->energy()*jet->extras.undoJEC)>0.9) continue;
        if (fabs(jet->eta())>9.9) continue;

        CMSLorentzVector const& jet_p4_nomus_uncor = jet->extras.momentum_nomus_uncor;
        CMSLorentzVector const jet_p4_nomus_cor = jet_p4_nomus_uncor*jet->extras.JEC_raw_nomus;

        if (use2017EENoiseRecipe && jet_p4_nomus_uncor.pt() < 50. && abs(jet_p4_nomus_uncor.eta()) >= 2.65 && abs(jet_p4_nomus_uncor.eta()) <= 3.139) continue; // Cuts on uncorected momentum
        if (jet_p4_nomus_cor.pt() <= 15.) continue; // Use corected momentum without muon contributions for the rest

        //unshiftedSum += jet_p4_nomus_uncor;
        shiftedSum_nominal += jet_p4_nomus_cor;
        shiftedSum_L1only += jet_p4_nomus_uncor*jet->extras.JEC_L1_raw_nomus;
        shiftedSum_JECup += jet_p4_nomus_cor*(1.f + jet->extras.JEC_raw_unc_nomus);
        shiftedSum_JECdn += jet_p4_nomus_cor*(1.f - jet->extras.JEC_raw_unc_nomus);
      }
      CMSLorentzVector v_met_JEC = v_met - shiftedSum_nominal + shiftedSum_L1only;
      CMSLorentzVector v_met_JECup = v_met - shiftedSum_JECup + shiftedSum_L1only;
      CMSLorentzVector v_met_JECdn = v_met - shiftedSum_JECdn + shiftedSum_L1only;
      metobj->extras.met
        = metobj->extras.met_METup = metobj->extras.met_METdn
        = metobj->extras.met_JERup = metobj->extras.met_JERdn
        = metobj->extras.met_PUup = metobj->extras.met_PUdn
        = v_met_JEC.Pt();
      metobj->extras.phi
        = metobj->extras.phi_JERup = metobj->extras.phi_JERdn
        = metobj->extras.phi_PUup = metobj->extras.phi_PUdn
        = v_met_JEC.Phi();
      metobj->extras.met_JECup = v_met_JECup.Pt();
      metobj->extras.phi_JECup = v_met_JECup.Phi();
      metobj->extras.met_JECdn = v_met_JECdn.Pt();
      metobj->extras.phi_JECdn = v_met_JECdn.Phi();
      if (this->verbosity>=TVar::DEBUG) MELAout
        << "JetMETHandler::applyJEC: Original MET vector: " << v_met
        << " | shifted MET vector: " << v_met_JEC
        << " | shifted MET vector up: " << v_met_JECup
        << " | shifted MET vector dn: " << v_met_JECdn
        << endl;
    }
  }

  //bool const doCorrectAK8Jets = !(SampleHelpers::theDataYear == 2016 && SampleHelpers::theDataVersion == SampleHelpers::kCMSSW_8_0_X);
  //bool const doCorrectAK8Jets = (SampleHelpers::theDataYear != 2016);
  constexpr bool doCorrectAK8Jets = true;
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

  if (registeredJERSFHandler_ak4jets && !ak4jets.empty()){
    float rho=0;
    for (auto*& jet:ak4jets){
      AK4JetVariables const& extras = jet->extras;
      rho = extras.rho;
      break;
    }
    registeredJERSFHandler_ak4jets->smear(ak4jets, rho, fwktree->isMC());
  }
  if (registeredJERSFHandler_ak8jets && !ak8jets.empty()){
    float rho = 0;
    for (auto*& jet:ak8jets){
      AK8JetVariables const& extras = jet->extras;
      rho = extras.rho;
      break;
    }
    registeredJERSFHandler_ak8jets->smear(ak8jets, rho, fwktree->isMC());
  }

  return true;
}
bool JetMETHandler::applySelections(){
  for (auto* obj:ak4jets){
    AK4JetSelectionHelpers::setSelectionBits(*obj);
    if (metobj) AK4JetSelectionHelpers::isBadMuonJet(*obj, *metobj);
  }
  for (auto* obj:ak8jets) AK8JetSelectionHelpers::setSelectionBits(*obj);
  return true;
}
bool JetMETHandler::applyBtaggingAndSFs(){
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree) return false;

  BtagScaleFactorHandler* theBTagSFHandler = nullptr;
  if (fwktree->isMC() && fwktree->isFastSim()) theBTagSFHandler = registeredBtagSFHandler_FastSim;
  else if (fwktree->isMC()) theBTagSFHandler = registeredBtagSFHandler;
  if (!theBTagSFHandler && fwktree->isMC()) return true;

  float bTaggerThreshold;
  if (theBTagSFHandler) bTaggerThreshold = theBTagSFHandler->WPval;
  else bTaggerThreshold = BtagHelpers::getBtagWP(AK4JetSelectionHelpers::AK4Jets_BTagWPType);

  for (auto*& jet:ak4jets){
    AK4JetVariables const& extras = jet->extras;

    float bTagger=0;
    switch (AK4JetSelectionHelpers::AK4Jets_BTagWPType){
    case BtagHelpers::kCSVv2_Loose:
    case BtagHelpers::kCSVv2_Medium:
    case BtagHelpers::kCSVv2_Tight:
      bTagger = extras.pfCombinedInclusiveSecondaryVertexV2BJetTag;
      break;
    case BtagHelpers::kDeepCSV_Loose:
    case BtagHelpers::kDeepCSV_Medium:
    case BtagHelpers::kDeepCSV_Tight:
      bTagger = extras.deepCSVb + extras.deepCSVbb;
      break;
    }

    bool isBtagged = (bTagger > bTaggerThreshold);
    bool isBtaggedWithSF   = isBtagged;
    bool isBtaggedWithSF_JECUp   = isBtagged;
    bool isBtaggedWithSF_JECDn   = isBtagged;
    bool isBtaggedWithSFUp = isBtagged;
    bool isBtaggedWithSFDn = isBtagged;

    if (theBTagSFHandler){
      int const& flav = extras.hadron_flavor;
      CMSLorentzVector correctedMomentum = jet->getFinalMomentum();
      float jpt = correctedMomentum.Pt();
      float jpt_jecup = jet->getCorrectedMomentum(+1).Pt();
      float jpt_jecdn = jet->getCorrectedMomentum(-1).Pt();
      float jphi = correctedMomentum.Phi();
      float jeta = correctedMomentum.Eta();
      TRandom3 rand;
      rand.SetSeed(std::abs(static_cast<int>(sin(jphi)*100000)));
      float R = rand.Uniform();
      float SF   = theBTagSFHandler->getSF(0, flav, jpt, jeta);
      float SF_JECUp   = theBTagSFHandler->getSF(0, flav, jpt_jecup, jeta);
      float SF_JECDn   = theBTagSFHandler->getSF(0, flav, jpt_jecdn, jeta);
      float SFUp = theBTagSFHandler->getSF(1, flav, jpt, jeta);
      float SFDn = theBTagSFHandler->getSF(-1, flav, jpt, jeta);
      float bTagMCEff = theBTagSFHandler->getEff(flav, jpt, jeta);
      if (SF  <=1.f && isBtagged && R<1.f-SF) isBtaggedWithSF   = false;
      if (SF_JECUp<=1.f && isBtagged && R<1.f-SF_JECUp) isBtaggedWithSF_JECUp = false;
      if (SF_JECDn<=1.f && isBtagged && R<1.f-SF_JECDn) isBtaggedWithSF_JECDn = false;
      if (SFUp<=1.f && isBtagged && R<1.f-SFUp) isBtaggedWithSFUp = false;
      if (SFDn<=1.f && isBtagged && R<1.f-SFDn) isBtaggedWithSFDn = false;
      if (SF  >1.f && !isBtagged && R<(1.f-SF)/(1.f-1.f/bTagMCEff)) isBtaggedWithSF   = true;
      if (SFUp>1.f && !isBtagged && R<(1.f-SFUp)/(1.f-1.f/bTagMCEff)) isBtaggedWithSFUp = true;
      if (SFDn>1.f && !isBtagged && R<(1.f-SFDn)/(1.f-1.f/bTagMCEff)) isBtaggedWithSFDn = true;
    }
    if (isBtaggedWithSF) jet->setSelectionBit(AK4JetSelectionHelpers::kIsBTagged);
    if (isBtaggedWithSF_JECUp) jet->setSelectionBit(AK4JetSelectionHelpers::kIsBTagged_JECUp);
    if (isBtaggedWithSF_JECDn) jet->setSelectionBit(AK4JetSelectionHelpers::kIsBTagged_JECDn);
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
  return JetMETHandler::getBtagValueFromLists(bDiscriminatorNames, btagvals.at(ijet), btagname);
}
float JetMETHandler::getBtagValueFromLists(std::vector<TString> const& bDiscriminatorNames, std::vector<float> const& btagvals, TString btagname){
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
  return btagvals.at(index);
}

void JetMETHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(tree);
  if (!fwktree) return;
  SampleHelpers::setupUsingOptions(fwktree->getOptions());

  // ak4 jet variables
  if (doAK4Jets){
    this->addConsumed<float>(_ak4jets_rho_);
    fwktree->bookEDMBranch<float>(_ak4jets_rho_, 0);

    // This variable is special, jsut like rho
    this->addConsumed<std::vector<TString>*>(_ak4jets_bDiscriminatorNames_);
    fwktree->bookEDMBranch<std::vector<TString>*>(_ak4jets_bDiscriminatorNames_, nullptr);

#define AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(TYPE, NAME) this->addConsumed<TYPE*>(_ak4jets_##NAME##_); fwktree->bookEDMBranch<TYPE*>(_ak4jets_##NAME##_, nullptr);
    AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE
  }

  // ak8 jet variables
  if (doAK8Jets){
    this->addConsumed<float>(_ak8jets_rho_);
    if (std::string(_ak4jets_rho_)!=std::string(_ak8jets_rho_)) fwktree->bookEDMBranch<float>(_ak8jets_rho_, 0);

#define AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE(TYPE, NAME) this->addConsumed<TYPE*>(_ak8jets_##NAME##_); fwktree->bookEDMBranch<TYPE*>(_ak8jets_##NAME##_, nullptr);
    AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVE
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
    this->addConsumed<float>(_pfmetraw_);
    this->addConsumed<float>(_pfmetrawPhi_);
    fwktree->bookEDMBranch<float>(_pfmet_, 0);
    fwktree->bookEDMBranch<float>(_pfmetPhi_, 0);
    fwktree->bookEDMBranch<float>(_pfmetraw_, 0);
    fwktree->bookEDMBranch<float>(_pfmetrawPhi_, 0);
    // Nowhere should we use "old" MET
    /*if (SampleHelpers::theDataYear == 2017 && SampleHelpers::theDataVersion == SampleHelpers::kCMSSW_9_4_X && fwktree->sampleIdentifier.Contains("09May2018")){
      this->addConsumed<float>(_pfmetraw_old_);
      this->addConsumed<float>(_pfmetrawPhi_old_);
      this->defineConsumedSloppy(_pfmetraw_old_);
      this->defineConsumedSloppy(_pfmetrawPhi_old_);
      fwktree->bookEDMBranch<float>(_pfmetraw_old_, 0);
      fwktree->bookEDMBranch<float>(_pfmetrawPhi_old_, 0);
    }
    else */if (SampleHelpers::theDataYear == 2016 && SampleHelpers::theDataVersion == SampleHelpers::kCMSSW_8_0_X){
      this->addConsumed<float>(_pfmetraw_muegclean_);
      this->addConsumed<float>(_pfmetrawPhi_muegclean_);
      this->defineConsumedSloppy(_pfmetraw_muegclean_);
      this->defineConsumedSloppy(_pfmetrawPhi_muegclean_);
      fwktree->bookEDMBranch<float>(_pfmetraw_muegclean_, 0);
      fwktree->bookEDMBranch<float>(_pfmetrawPhi_muegclean_, 0);
    }
  }
}


#undef AK4JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef AK8JETS_VECTOR_ITERATOR_HANDLER_DIRECTIVES
