#include <cassert>
#include "TRandom3.h"
#include "FrameworkVariables.hh"
#include "FrameworkSet.h"
#include "EventFilterHandler.h"
#include "GoodEventFilter.h"
#include "ElectronSelectionHelpers.h"
#include "PhotonSelectionHelpers.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


EventFilterHandler::EventFilterHandler() :
  IvyBase(),
  passEventFilterFlag(true)
{
  // HLT variables
  this->addConsumed<TBits>(_hlt_bits_);
  this->addConsumed<std::vector<unsigned int>*>(_hlt_prescales_);
  this->addConsumed<std::vector<unsigned int>*>(_hlt_l1prescales_);
  this->addConsumed<std::vector<TString>*>(_hlt_trigNames_);

  // Vertex variables
  this->addConsumed<std::vector<int>*>(_vertex_isFake_);
  this->addConsumed<std::vector<float>*>(_vertex_Ndof_);
  this->addConsumed<std::vector<CMSLorentzVector>*>(_vertex_positions_);
}

void EventFilterHandler::clear(){
  passEventFilterFlag=true;
  product_HLTpaths.clear();
}

bool EventFilterHandler::constructFilter(){
  clear();
  if (!currentTree){
    if (verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::constructFilter: Current tree is null!" << endl;
    return false;
  }
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree){
    if (verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::constructFilter: Current tree is not derived from a FrameworkTree class!" << endl;
    return false;
  }
  bool allVariablesPresent = true;

  // HLT trigger
  if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << "EventFilterHandler::constructFilter: Checking HLT paths..." << endl;
  bool passAtLeastOneTrigger = false;
  TBits hlt_bits = TBits();
  std::vector<unsigned int>* hlt_prescales = nullptr;
  std::vector<unsigned int>* hlt_l1prescales = nullptr;
  std::vector<TString>* hlt_trigNames = nullptr;
  allVariablesPresent &= (
    this->getConsumedValue(_hlt_bits_, hlt_bits)
    &&
    this->getConsumedValue(_hlt_prescales_, hlt_prescales)
    &&
    this->getConsumedValue(_hlt_l1prescales_, hlt_l1prescales)
    &&
    this->getConsumedValue(_hlt_trigNames_, hlt_trigNames)
    );
  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "EventFilterHandler::constructFilter: Not all trigger variables are consumed properly!" << endl;
    assert(0);
  }
  bool skipTriggerChecks=false;
  if (!hlt_trigNames || !hlt_l1prescales || !hlt_prescales){
    if (this->verbosity>=TVar::DEBUG){
      if (!hlt_trigNames) MELAerr << "EventFilterHandler::constructFilter: Trigger names pointer is null, but this issue should have been caught earlier!" << endl;
      if (!hlt_l1prescales) MELAerr << "EventFilterHandler::constructFilter: L1 prescales pointer is null, but this issue should have been caught earlier!" << endl;
      if (!hlt_prescales) MELAerr << "EventFilterHandler::constructFilter: HLT prescales pointer is null, but this issue should have been caught earlier!" << endl;
    }
    skipTriggerChecks = true;
  }
  std::unordered_map<TString, std::vector<TString>> HLTPathMap = EventFilterHandler::getHLTPaths(fwktree);
  for (auto const& pair:HLTPathMap){
    TString HLTPathType = pair.first;
    bool doHLTPrescale = HLTPathType.Contains("HLTprescale");
    bool doL1Prescale = HLTPathType.Contains("L1prescale");
    std::vector<TString> const& HLTPathsToTest = pair.second;
    if (HLTPathType.Contains(":")){
      TString wish, value;
      HelperFunctions::splitOption(HLTPathType, wish, value, ':');
      // Don't really care about the option, already tested it...
      HLTPathType = wish;
    }

    if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << "EventFilterHandler::constructFilter: Checking HLT path type " << HLTPathType << endl;

    std::unordered_map<TString, bool>::iterator it_product = product_HLTpaths.find(HLTPathType);
    if (it_product == product_HLTpaths.end()){
      product_HLTpaths[HLTPathType] = false;
      it_product = product_HLTpaths.find(HLTPathType);
    }
    for (TString const& hltpath:HLTPathsToTest){
      if (skipTriggerChecks) break;
      if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << "\t- Testing HLT path " << hltpath << endl;
      int trigIndx = -1;
      std::vector<TString>::const_iterator begin_it = hlt_trigNames->cbegin();
      std::vector<TString>::const_iterator end_it = hlt_trigNames->cend();
      std::vector<TString>::const_iterator found_it = end_it;
      for (auto it=begin_it; it!=end_it; it++){
        if (it->Contains(hltpath)){
          found_it=it;
          break;
        }
      }
      if (found_it != end_it) trigIndx = found_it - begin_it;
      else continue;
      bool testBit = hlt_bits.TestBitNumber(trigIndx);
      if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << "\t- HLT path " << *found_it << " = " << testBit;
      if (doHLTPrescale) testBit &= (hlt_prescales->at(trigIndx)>0);
      if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << " (after HLT prescale: " << testBit;
      if (doL1Prescale) testBit &= (hlt_l1prescales->at(trigIndx)>0);
      if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << ", after L1 prescale: " << testBit << ")" << endl;
      passAtLeastOneTrigger |= testBit;
      it_product->second |= testBit;
    }
  }
  passEventFilterFlag &= passAtLeastOneTrigger;
  if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << "\t- Event flag is " << passEventFilterFlag << " after HLT path filter." << endl;
  if (!passEventFilterFlag) return true; // No need to proceed further

  // PV filter
  if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << "EventFilterHandler::constructFilter: Checking valid PV..." << endl;
  std::vector<int>::const_iterator itBegin_vertex_isFake;
  std::vector<float>::const_iterator itBegin_vertex_Ndof;
  std::vector<CMSLorentzVector>::const_iterator itBegin_vertex_positions, itEnd_vertex_positions;
  allVariablesPresent &= (
    this->getConsumedCIterators<std::vector<int>>(_vertex_isFake_, &itBegin_vertex_isFake)
    && this->getConsumedCIterators<std::vector<float>>(_vertex_Ndof_, &itBegin_vertex_Ndof)
    && this->getConsumedCIterators<std::vector<CMSLorentzVector>>(_vertex_positions_, &itBegin_vertex_positions, &itEnd_vertex_positions)
    );
  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "EventFilterHandler::constructFilter: Not all vertex variables are consumed properly!" << endl;
    assert(0);
  }
  { // Test for at least one good vertex
    bool hasGoodVertex = false;
    auto it_isFake = itBegin_vertex_isFake;
    auto it_Ndof = itBegin_vertex_Ndof;
    auto it_position = itBegin_vertex_positions;
    while (it_position!=itEnd_vertex_positions){
      hasGoodVertex |= (
        !(*it_isFake)
        && (*it_Ndof > 4.)
        && (it_position->Rho() <= 2.f)
        && (fabs(it_position->Z()) <= 24.f)
        );
      if (hasGoodVertex) break;
      it_isFake++;
      it_Ndof++;
      it_position++;
    }
    passEventFilterFlag &= hasGoodVertex;
    if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << "\t- Event flag is " << passEventFilterFlag << " after PV filter." << endl;
    if (!passEventFilterFlag) return true; // No need to proceed further
  }

  // MET filters
  if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << "EventFilterHandler::constructFilter: Checking MET filters..." << endl;
  std::vector<TString> metfilters = EventFilterHandler::getMETFilterFlags(fwktree);
  std::vector<bool> metfiltervals(metfilters.size(), true);
  {
    std::vector<TString>::iterator itS=metfilters.begin();
    std::vector<bool>::iterator itV=metfiltervals.begin();
    for (; itS!=metfilters.end(); itS++){
      bool value;
      allVariablesPresent &= this->getConsumedValue((*itS), value);
      *itV = value;
      itV++;
    }
  }
  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "EventFilterHandler::constructFilter: Not all MET variables are consumed properly!" << endl;
    assert(0);
  }
  for (bool const& flag:metfiltervals) passEventFilterFlag &= flag;
  if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << "\t- Event flag is " << passEventFilterFlag << " after MET filters." << endl;
  if (!passEventFilterFlag) return true; // No need to proceed further

  // Test GoodEventFilter for the JSON file
  if (fwktree->isData()){
    if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << "EventFilterHandler::constructFilter: Checking the JSON file..." << endl;
    RunNumber_t run = 0;
    Lumisection_t ls = 0;
    EventNumber_t evt = 0;
    allVariablesPresent &= (
      this->getConsumedValue(_event_RunNumber_, run)
      && this->getConsumedValue(_event_Lumisection_, ls)
      && this->getConsumedValue(_event_EventNumber_, evt)
      );
    if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
      MELAerr << "EventFilterHandler::constructFilter: Not all data variables are consumed properly!" << endl;
      assert(0);
    }
    passEventFilterFlag &= GoodEventFilter::testEvent(run, ls);
    if (verbosity>=TVar::DEBUG_VERBOSE) MELAout << "\t- Event flag is " << passEventFilterFlag << " after the JSON filter." << endl;
    if (!passEventFilterFlag) return true; // No need to proceed further
  }

  return true;
}

bool EventFilterHandler::test2018HEMFilter(
  std::vector<ElectronObject*> const* electrons,
  std::vector<PhotonObject*> const* photons,
  std::vector<AK4JetObject*> const* ak4jets,
  std::vector<AK8JetObject*> const* ak8jets,
  SystematicsHelpers::SystematicVariationTypes syst
) const{
  if (SampleHelpers::theDataYear != 2018) return true;

  // Do not run clear because this is a special filter that does not modify the actual class
  if (!currentTree){
    if (verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::test2018HEMFilter: Current tree is null!" << endl;
    return false;
  }
  FrameworkTree* const fwktree = dynamic_cast<FrameworkTree* const>(currentTree);
  if (!fwktree){
    if (verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::test2018HEMFilter: Current tree is not derived from a FrameworkTree class!" << endl;
    return false;
  }

  bool allVariablesPresent = true;
  RunNumber_t run = 0;
  Lumisection_t ls = 0;
  EventNumber_t evt = 0;
  allVariablesPresent &= (
    this->getConsumedValue(_event_RunNumber_, run)
    && this->getConsumedValue(_event_Lumisection_, ls)
    && this->getConsumedValue(_event_EventNumber_, evt)
    );
  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "EventFilterHandler::test2018HEMFilter: Not all data variables are consumed properly!" << endl;
    assert(0);
  }

  bool checkVeto = false;
  if (fwktree->isMC()){
    TRandom3 rand;
    unsigned long const seed = (run>0 ? run : 0) + (ls>0 ? ls : 0) + (evt>0 ? evt : 0);
    assert(seed>0);
    rand.SetSeed(seed);
    float veto_x = rand.Uniform();
    static const float lumi_affected = SampleHelpers::getIntegratedLuminosity("2018_HEMaffected");
    static const float lumi_total = SampleHelpers::getIntegratedLuminosity("2018");
    checkVeto = (veto_x<lumi_affected/lumi_total);
  }
  else checkVeto = (run>=319077);

  if (!checkVeto) return true;

  using namespace SystematicsHelpers;
  if (syst != sNominal && fwktree->isData()) syst = sNominal;

  bool doVeto = false;
  static const std::pair<float, float> eta_region(-4.7, -1.4);
  static const std::pair<float, float> phi_region(-1.6, -0.8);

  if (!doVeto && electrons){
    for (auto const* part:(*electrons)){
      if (!part->testSelection(ElectronSelectionHelpers::kVetoIDReco)) continue;
      float const eta = part->eta();
      float const phi = part->phi();
      if (eta>=eta_region.first && eta<=eta_region.second && phi>=phi_region.first && phi<=phi_region.second){ doVeto=true; break; }
    }
    if (verbosity>=TVar::DEBUG && doVeto) MELAout << "EventFilterHandler::test2018HEMFilter: Found at least one electron satisfying HEM15/16." << endl;
  }
  if (!doVeto && photons){
    for (auto const* part:(*photons)){
      if (!part->testSelection(PhotonSelectionHelpers::kLooseIDReco)) continue;
      float const eta = part->eta();
      float const phi = part->phi();
      if (eta>=eta_region.first && eta<=eta_region.second && phi>=phi_region.first && phi<=phi_region.second){ doVeto=true; break; }
    }
    if (verbosity>=TVar::DEBUG && doVeto) MELAout << "EventFilterHandler::test2018HEMFilter: Found at least one photon satisfying HEM15/16." << endl;
  }
  // Require a pT>30 GeV cut on jets
  if (!doVeto && ak4jets){
    for (auto const* part:(*ak4jets)){
      bool isSelectedJet=false;
      if (syst == eJECUp) isSelectedJet = (part->testSelection(AK4JetSelectionHelpers::kTightID) && part->getCorrectedMomentum(1).Pt()>=30.f);
      else if (syst == eJECDn) isSelectedJet = (part->testSelection(AK4JetSelectionHelpers::kTightID) && part->getCorrectedMomentum(-1).Pt()>=30.f);
      else if (syst == eJERUp) isSelectedJet = (part->testSelection(AK4JetSelectionHelpers::kTightID) && part->getCorrectedMomentum(2).Pt()>=30.f);
      else if (syst == eJERDn) isSelectedJet = (part->testSelection(AK4JetSelectionHelpers::kTightID) && part->getCorrectedMomentum(-2).Pt()>=30.f);
      else isSelectedJet = (part->testSelection(AK4JetSelectionHelpers::kTightID) && part->getCorrectedMomentum(0).Pt()>=30.f);

      if (!isSelectedJet) continue;

      float const eta = part->eta();
      float const phi = part->phi();
      if (eta>=eta_region.first && eta<=eta_region.second && phi>=phi_region.first && phi<=phi_region.second){ doVeto=true; break; }
    }
    if (verbosity>=TVar::DEBUG && doVeto) MELAout << "EventFilterHandler::test2018HEMFilter: Found at least one AK4 jet satisfying HEM15/16." << endl;
  }
  // Be careful! There is no equivalent of tight ID in ak8 jets, so testing is done on all jets
  if (!doVeto && ak8jets){
    for (auto const* part:(*ak8jets)){
      bool isSelectedJet=false;
      if (syst == eJECUp) isSelectedJet = (part->getCorrectedMomentum(1).Pt()>=30.f);
      else if (syst == eJECDn) isSelectedJet = (part->getCorrectedMomentum(-1).Pt()>=30.f);
      else if (syst == eJERUp) isSelectedJet = (part->getCorrectedMomentum(2).Pt()>=30.f);
      else if (syst == eJERDn) isSelectedJet = (part->getCorrectedMomentum(-2).Pt()>=30.f);
      else isSelectedJet = (part->getCorrectedMomentum(0).Pt()>=30.f);

      if (!isSelectedJet) continue;

      float const eta = part->eta();
      float const phi = part->phi();
      if (eta>=eta_region.first && eta<=eta_region.second && phi>=phi_region.first && phi<=phi_region.second){ doVeto=true; break; }
    }
    if (verbosity>=TVar::DEBUG && doVeto) MELAout << "EventFilterHandler::test2018HEMFilter: Found at least one AK8 jet satisfying HEM15/16." << endl;
  }

  return !doVeto;
}


void EventFilterHandler::bookBranches(BaseTree* tree){
  if (!tree || !tree->isValid()) return;
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(tree);
  if (!fwktree) return;
  SampleHelpers::setupUsingOptions(fwktree->getOptions());

  // HLT variables (non-sloppy)
  fwktree->bookEDMBranch<TBits>(_hlt_bits_, TBits());
  fwktree->bookEDMBranch<std::vector<unsigned int>*>(_hlt_prescales_, nullptr);
  fwktree->bookEDMBranch<std::vector<unsigned int>*>(_hlt_l1prescales_, nullptr);
  fwktree->bookEDMBranch<std::vector<TString>*>(_hlt_trigNames_, nullptr);

  // Vertex variables
  fwktree->bookEDMBranch<std::vector<int>*>(_vertex_isFake_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_vertex_Ndof_, nullptr);
  fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_vertex_positions_, nullptr);

  // MET filters
  std::vector<TString> metfilters = EventFilterHandler::getMETFilterFlags(fwktree);
  for (TString const& filt:metfilters){
    fwktree->bookEDMBranch<bool>(filt, true); // Default should be true to avoid non-existing branches
    this->addConsumed<bool>(filt);
    this->defineConsumedSloppy(filt); // Define as sloppy so that different samples from different years/versions can be processed.
  }
  
  // For data, will also need access to RunNumber, LumiSection and EventNumber
  // This is needed to test for the JSON file
  if (fwktree->isData() || SampleHelpers::theDataYear == 2018){ // 2018 MC also needs these for the HEM filter seed
    fwktree->bookEDMBranch<RunNumber_t>(_event_RunNumber_, 0);
    this->addConsumed<RunNumber_t>(_event_RunNumber_);
    this->defineConsumedSloppy(_event_RunNumber_);

    fwktree->bookEDMBranch<Lumisection_t>(_event_Lumisection_, 0);
    this->addConsumed<Lumisection_t>(_event_Lumisection_);
    this->defineConsumedSloppy(_event_Lumisection_);

    fwktree->bookEDMBranch<EventNumber_t>(_event_EventNumber_, 0);
    this->addConsumed<EventNumber_t>(_event_EventNumber_);
    this->defineConsumedSloppy(_event_EventNumber_);
  }
}
std::vector<TString> EventFilterHandler::getMETFilterFlags(FrameworkTree* fwktree){
  std::vector<TString> metfilters;
  if (!fwktree) return metfilters;

  // MET filters
  // Recommendations are from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
  if (SampleHelpers::theDataYear == 2016){
    metfilters=std::vector<TString>{
      _filt_goodVertices_,
      _filt_hbheNoise_,
      _filt_hbheNoiseIso_,
      _filt_ecalTP_
    };
    if (!fwktree->isFastSim()) metfilters.push_back(_filt_globalSuperTightHalo2016_); // For data or non-FS MC
    if (fwktree->isData()) metfilters.push_back(_filt_eeBadSc_); // Only for data

    if (SampleHelpers::theDataVersion >= SampleHelpers::kCMSSW_9_4_X){ // These MET filters are available in CMSSW_VERSION>=94X
      metfilters.push_back(_filt_BadPFMuon_filter_);
      metfilters.push_back(_filt_BadChargedCandidate_filter_);
    }
    // Else need "Bad PF Muon Filter" and "Bad Charged Hadron Filter" to be calculated on the fly for data, MC and FastSim
    // FIXME: Needs some branches for PFcands and more for muons. Ignore for the moment...
  }
  else if (SampleHelpers::theDataYear == 2017){
    metfilters=std::vector<TString>{
      _filt_goodVertices_,
      _filt_hbheNoise_,
      _filt_hbheNoiseIso_,
      _filt_ecalTP_,
      _filt_BadPFMuon_filter_,
      //_filt_BadChargedCandidate_filter_, // FIXME: To be updated following https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2017_data
      _filt_ecalBadCalib_filterUpdate_
    };
    if (!fwktree->isFastSim()) metfilters.push_back(_filt_globalSuperTightHalo2016_); // For data or non-FS MC
    if (fwktree->isData()) metfilters.push_back(_filt_eeBadSc_); // Only for data
  }
  else if (SampleHelpers::theDataYear == 2018){
    metfilters=std::vector<TString>{
      _filt_goodVertices_,
      _filt_hbheNoise_,
      _filt_hbheNoiseIso_,
      _filt_ecalTP_,
      _filt_BadPFMuon_filter_,
      //_filt_BadChargedCandidate_filter_, // FIXME: To be updated following https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2018_data
      _filt_ecalBadCalib_filterUpdate_
    };
    if (!fwktree->isFastSim()) metfilters.push_back(_filt_globalSuperTightHalo2016_); // For data or non-FS MC
    if (fwktree->isData()) metfilters.push_back(_filt_eeBadSc_); // Only for data
  }

  return metfilters;
}
std::unordered_map<TString, std::vector<TString>> EventFilterHandler::getHLTPaths(FrameworkTree* fwktree){
  std::unordered_map<TString, std::vector<TString>> res;
  if (!fwktree) return res;


  if (SampleHelpers::theDataYear == 2016){
    res["SingleEle"] ={
      "HLT_Ele25_eta2p1_WPTight_Gsf_v", "HLT_Ele27_eta2p1_WPTight_Gsf_v", "HLT_Ele27_WP85_Gsf_v", "HLT_Ele27_eta2p1_WPLoose_Gsf_v"
    };
    res["DoubleEle"] ={
      "HLT_DoubleEle33_CaloIdL_MW_v", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
    res["SingleMu"] ={
      "HLT_IsoMu20_v", "HLT_IsoTkMu20_v", "HLT_IsoMu22_v", "HLT_IsoTkMu22_v", "HLT_IsoMu24_v", "HLT_IsoTkMu24_v"
    };
    res["DoubleMu"] ={
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"
    };
    res["MuEle"] ={
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"
    };
    res["MET_MHT"] ={
      "HLT_PFMET100_PFMHT100_IDTight_v", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v",
      "HLT_PFMET110_PFMHT110_IDTight_v", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v",
      "HLT_PFMET120_PFMHT120_IDTight_v", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
      "HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v", "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v",
      //
      // disabled most of the time
      "HLT_PFMET100_PFMHT100_IDTight_v", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v",
      //
      "HLT_PFMET110_PFMHT110_IDTight_v", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v",
      //
      "HLT_PFMET120_PFMHT120_IDTight_v", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v", "HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v", "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v"
    };
    res["MET"] ={
      "HLT_PFMET170_NoiseCleaned_v",
      "HLT_PFMET170_JetIdCleaned_v",
      "HLT_PFMET170_HBHECleaned_v",
      "HLT_PFMET170_NotCleaned_v"
    };
    // HT triggers for trigger efficiency - from MT2
    res["HT"] ={
      "HLT_PFHT800_v", "HLT_PFHT900_v",
      "HLT_PFHT125_v", "HLT_PFHT200_v",
      "HLT_PFHT300_v", "HLT_PFHT350_v",
      "HLT_PFHT475_v", "HLT_PFHT600_v"
    };
    //as we use those only for trigger efficiency measurements, we actually don't care about the exact prescale ...
    res["AK8Jet"] ={
      "HLT_AK8PFJet360_TrimMass30_v", "HLT_AK8PFJet400_TrimMass30_v",
      "HLT_AK8PFJet500_v", "HLT_AK8PFJet450_v",
      "HLT_AK8PFJet260_v", "HLT_AK8PFJet320_v", "HLT_AK8PFJet400_v"
    };
  }
  else if (SampleHelpers::theDataYear == 2017){
    res["SingleEle"] ={
      "HLT_Ele35_WPTight_Gsf_v", "HLT_Ele32_WPTight_Gsf_v", "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"
    };
    res["DoubleEle"] ={
      "HLT_DoubleEle33_CaloIdL_MW_v", "HLT_DoubleEle25_CaloIdL_MW_v", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"
    };
    res["SingleMu"] ={
      "HLT_IsoMu27_v", "HLT_IsoMu24_v"
    };
    res["DoubleMu"] ={
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"
    };
    res["MuEle"] ={
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"
    };
    res["MET_MHT"] ={
      "HLT_PFMET120_PFMHT120_IDTight_v", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v", "HLT_PFMET100_PFMHT100_IDTight_PFHT60_v", "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",
      //
      "HLT_PFMET110_PFMHT110_IDTight_v", "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v",
      //
      "HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v", "HLT_PFMET120_PFMHT120_IDTight_v",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v",
      "HLT_PFMETTypeOne120_PFMHT120_IDTight_HFCleaned_v", "HLT_PFMETTypeOne120_PFMHT120_IDTight_v",
      //
      "HLT_PFMET130_PFMHT130_IDTight_v",
    };
    res["MET"] ={
      "HLT_PFMET250_HBHECleaned_v", "HLT_PFMET300_HBHECleaned_v", "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v"
    };
    res["HT"] ={
      "HLT_PFHT1050_v", // Unprescaled
      "HLT_PFHT430_v", "HLT_PFHT510_v", "HLT_PFHT590_v", "HLT_PFHT680_v", "HLT_PFHT780_v", "HLT_PFHT890_v"
    };
    res["AK8Jet"] ={
      // Unprescaled
      "HLT_AK8PFJet360_TrimMass30_v", "HLT_AK8PFJet380_TrimMass30_v",
      "HLT_AK8PFJet400_TrimMass30_v", "HLT_AK8PFJet420_TrimMass30_v",
      "HLT_AK8PFJet500_v", "HLT_AK8PFJet550_v",
      // Prescalsed
      "HLT_AK8PFJet260_v", "HLT_AK8PFJet320_v",
      "HLT_AK8PFJet400_v", "HLT_AK8PFJet450_v"
    };
  }
  else if (SampleHelpers::theDataYear == 2018){
    res["SingleEle"] ={
      "HLT_Ele35_WPTight_Gsf_v",
      "HLT_Ele32_WPTight_Gsf_v"
    };
    res["DoubleEle"] ={
      "HLT_DoubleEle33_CaloIdL_MW_v",
      "HLT_DoubleEle25_CaloIdL_MW_v",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"
    };
    res["SingleMu"] ={
      "HLT_IsoMu27_v",
      "HLT_IsoMu24_v"
    };
    res["DoubleMu"] ={
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"
    };
    res["MuEle"] ={
      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"
    };
    res["MET_MHT"] ={
      "HLT_PFMET120_PFMHT120_IDTight_v",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
      "HLT_PFMET100_PFMHT100_IDTight_PFHT60_v",
      "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",
      //
      "HLT_PFMET110_PFMHT110_IDTight_v",
      "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v",
      //
      "HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v",
      "HLT_PFMET120_PFMHT120_IDTight_v",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
      "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v",
      "HLT_PFMETTypeOne120_PFMHT120_IDTight_HFCleaned_v",
      "HLT_PFMETTypeOne120_PFMHT120_IDTight_v",
      //
      "HLT_PFMET130_PFMHT130_IDTight_v"
    };
    res["MET"] ={
      "HLT_PFMET250_HBHECleaned_v",
      "HLT_PFMET300_HBHECleaned_v",
      "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v"
    };
    res["HT"] ={
      // Unprescaled
      "HLT_PFHT1050_v",
      // Prescaled
      "HLT_PFHT430_v", "HLT_PFHT510_v",
      "HLT_PFHT590_v", "HLT_PFHT680_v",
      "HLT_PFHT780_v", "HLT_PFHT890_v"
    };
    res["AK8Jet"] ={
      // Unprescaled
      "HLT_AK8PFJet400_TrimMass30_v", "HLT_AK8PFJet420_TrimMass30_v",
      "HLT_AK8PFJet500_v", "HLT_AK8PFJet550_v",
      // Prescaled
      "HLT_AK8PFJet260_v", "HLT_AK8PFJet320_v",
      "HLT_AK8PFJet400_v", "HLT_AK8PFJet450_v"
    };
  }
  // Additional unprescaled jet trigger in JetHT
  res["CaloJet"] ={ "HLT_CaloJet500_NoJetID_v" };
  // Photon triggers, store individual ones for all 3 years
  res["SinglePhoton:HLTprescale,L1prescale"] ={
    "HLT_Photon50_R9Id90_HE10_IsoM_v",
    "HLT_Photon75_R9Id90_HE10_IsoM_v",
    "HLT_Photon90_R9Id90_HE10_IsoM_v",
    "HLT_Photon120_R9Id90_HE10_IsoM_v",
    "HLT_Photon120_v",
    "HLT_Photon165_R9Id90_HE10_IsoM_v",
    "HLT_Photon165_HE10_v",
    "HLT_Photon175_v"
  };
  // Unprescaled triggers in SinglePhoton
  res["SinglePhoton"] ={
    "HLT_Photon200_v",      // in 2017 & 2018 only
    "HLT_Photon250_NoHE_v", // in 2016 only
    "HLT_Photon300_NoHE_v"
  };

  return res;
}
