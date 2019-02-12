#include "GenericEventAnalyzer.h"
#include "HelperFunctions.h"
#include "WeightsHandler.h"
#include "GenInfoHandler.h"
#include "PFCandHandler.h"
#include "MuonHandler.h"
#include "MuonSelectionHelpers.h"
#include "MuonScaleFactorHandler.h"
#include "ElectronHandler.h"
#include "ElectronSelectionHelpers.h"
#include "ElectronScaleFactorHandler.h"
#include "PhotonHandler.h"
#include "PhotonSelectionHelpers.h"
//#include "PhotonScaleFactorHandler.h"
#include "JetMETHandler.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "TopnessCalculator.h"
#include "EventFilterHandler.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


GenericEventAnalyzer::GenericEventAnalyzer() :
  FrameworkTreeLooperBase(),
  doWeights(true),
  doGenInfo(true),
  doGenParticles(true),
  doEventFilter(true),
  doPFCands(true),
  doMuons(true),
  doElectrons(true),
  doPhotons(true),
  doJetMET(true),
  doWriteSelectionVariables(true)
{}
GenericEventAnalyzer::GenericEventAnalyzer(FrameworkTree* inTree) :
  FrameworkTreeLooperBase(inTree),
  doWeights(true),
  doGenInfo(true),
  doGenParticles(true),
  doEventFilter(true),
  doPFCands(true),
  doMuons(true),
  doElectrons(true),
  doPhotons(true),
  doJetMET(true),
  doWriteSelectionVariables(true)
{}
GenericEventAnalyzer::GenericEventAnalyzer(std::vector<FrameworkTree*> const& inTreeList) :
  FrameworkTreeLooperBase(inTreeList),
  doWeights(true),
  doGenInfo(true),
  doGenParticles(true),
  doEventFilter(true),
  doPFCands(true),
  doMuons(true),
  doElectrons(true),
  doPhotons(true),
  doJetMET(true),
  doWriteSelectionVariables(true)
{}
GenericEventAnalyzer::GenericEventAnalyzer(FrameworkSet const* inTreeSet) :
  FrameworkTreeLooperBase(inTreeSet),
  doWeights(true),
  doGenInfo(true),
  doGenParticles(true),
  doEventFilter(true),
  doPFCands(true),
  doMuons(true),
  doElectrons(true),
  doPhotons(true),
  doJetMET(true),
  doWriteSelectionVariables(true)
{}

bool GenericEventAnalyzer::runEvent(FrameworkTree* tree, float const& externalWgt, SimpleEntry& product){
  bool validProducts = (tree!=nullptr);
  if (!validProducts) return validProducts;

  // Ivy handlers
  WeightsHandler* wgtHandler=nullptr;
  GenInfoHandler* genInfoHandler = nullptr;
  EventFilterHandler* eventFilter = nullptr;
  PFCandHandler* pfcandHandler=nullptr;
  MuonHandler* muonHandler=nullptr;
  ElectronHandler* eleHandler=nullptr;
  PhotonHandler* photonHandler=nullptr;
  JetMETHandler* jetHandler=nullptr;
  for (auto it=this->externalIvyObjects.begin(); it!=this->externalIvyObjects.end(); it++){
    if (doWeights && !tree->isData()){ WeightsHandler* tmp_ivy = dynamic_cast<WeightsHandler*>(it->second); if (!wgtHandler && tmp_ivy){ wgtHandler = tmp_ivy; } }
    if (doGenInfo && !tree->isData()){ GenInfoHandler* tmp_ivy = dynamic_cast<GenInfoHandler*>(it->second); if (!genInfoHandler && tmp_ivy){ genInfoHandler = tmp_ivy; } }
    if (doEventFilter){ EventFilterHandler* tmp_ivy = dynamic_cast<EventFilterHandler*>(it->second); if (!eventFilter && tmp_ivy){ eventFilter = tmp_ivy; } }
    if (doPFCands){ PFCandHandler* tmp_ivy = dynamic_cast<PFCandHandler*>(it->second); if (!pfcandHandler && tmp_ivy){ pfcandHandler = tmp_ivy; } }
    if (doMuons){ MuonHandler* tmp_ivy = dynamic_cast<MuonHandler*>(it->second); if (!muonHandler && tmp_ivy){ muonHandler = tmp_ivy; } }
    if (doElectrons){ ElectronHandler* tmp_ivy = dynamic_cast<ElectronHandler*>(it->second); if (!eleHandler && tmp_ivy){ eleHandler = tmp_ivy; } }
    if (doPhotons){ PhotonHandler* tmp_ivy = dynamic_cast<PhotonHandler*>(it->second); if (!photonHandler && tmp_ivy){ photonHandler = tmp_ivy; } }
    if (doJetMET){ JetMETHandler* tmp_ivy = dynamic_cast<JetMETHandler*>(it->second); if (!jetHandler && tmp_ivy){ jetHandler = tmp_ivy; } }
  }

  // SF handlers
  MuonScaleFactorHandler* muonSFHandler=nullptr;
  ElectronScaleFactorHandler* eleSFHandler=nullptr;
  for (auto it=this->externalScaleFactorHandlers.begin(); it!=this->externalScaleFactorHandlers.end(); it++){
    if (doMuons){ MuonScaleFactorHandler* muon_sfs = dynamic_cast<MuonScaleFactorHandler*>(it->second); if (!muonSFHandler && muon_sfs){ muonSFHandler = muon_sfs; } }
    if (doElectrons){ ElectronScaleFactorHandler* ele_sfs = dynamic_cast<ElectronScaleFactorHandler*>(it->second); if (!eleSFHandler && ele_sfs){ eleSFHandler = ele_sfs; } }
  }


  /*********************/
  /**  WEIGHTS BLOCK  **/
  /*********************/
  float wgt = externalWgt;
  validProducts &= (!doWeights || wgtHandler!=nullptr || tree->isData());
  if (!validProducts){
    MELAerr << "GenericEventAnalyzer::runEvent: Weight handle is invalid (Tree: " << tree->sampleIdentifier << ")." << endl;
    return validProducts;
  }
  if (wgtHandler){
    validProducts &= wgtHandler->constructWeights();
    if (!validProducts){
      MELAerr << "GenericEventAnalyzer::runEvent: Weight product could not be constructed (Tree: " << tree->sampleIdentifier << ")." << endl;
      return validProducts;
    }
    validProducts &= wgtHandler->recordWeights(product, externalWgt);
    product.getNamedVal<float>(WeightVariables::getWeightName(WeightVariables::wCentral), wgt);
  }
  else product.setNamedVal<float>(WeightVariables::getWeightName(WeightVariables::wCentral), wgt);
  if (std::isnan(wgt) || std::isinf(wgt) || wgt==0.f){
    if (wgt!=0.f){
      MELAerr << "GenericEventAnalyzer::runEvent: Invalid weight " << wgt << " is being discarded (Tree " << tree->sampleIdentifier << ")." << endl;
      exit(1);
    }
    validProducts=false;
  }
  if (!validProducts) return validProducts;


  /*********************/
  /**  GENINFO BLOCK  **/
  /*********************/
  validProducts &= (!doGenInfo || genInfoHandler!=nullptr || tree->isData());
  if (!validProducts){
    MELAerr << "GenericEventAnalyzer::runEvent: Weight handle is invalid (Tree: " << tree->sampleIdentifier << ")." << endl;
    return validProducts;
  }
  if (genInfoHandler){
    validProducts &= genInfoHandler->constructGenInfo();
    if (!validProducts){
      MELAerr << "GenericEventAnalyzer::runEvent: Weight product could not be constructed (Tree: " << tree->sampleIdentifier << ")." << endl;
      return validProducts;
    }

    GenEventInfo const* genInfo = genInfoHandler->getGenInfo();
    if (genInfo){
      product.setNamedVal("pid", genInfo->processID);
      product.setNamedVal("qScale", genInfo->qscale);
      product.setNamedVal("alphaS", genInfo->alphaS);
      product.setNamedVal("genMET", genInfo->genMET);
      product.setNamedVal("genMETPhi", genInfo->genMETPhi);
      product.setNamedVal("xsec", genInfo->xsec);
    }
    if (genInfoHandler->getParticleInfoFlag() && doGenParticles){
      auto const& genparts = genInfoHandler->getGenParticles();

      std::vector<bool> genparticles_isPromptFinalState;
      std::vector<bool> genparticles_isPromptDecayed;
      std::vector<bool> genparticles_isDirectPromptTauDecayProductFinalState;
      std::vector<bool> genparticles_isHardProcess;
      std::vector<bool> genparticles_fromHardProcessFinalState;
      std::vector<bool> genparticles_fromHardProcessDecayed;
      std::vector<bool> genparticles_isDirectHardProcessTauDecayProductFinalState;
      std::vector<bool> genparticles_fromHardProcessBeforeFSR;
      std::vector<bool> genparticles_isLastCopy;
      std::vector<bool> genparticles_isLastCopyBeforeFSR;

      std::vector<int> genparticles_id;
      std::vector<int> genparticles_status;

      std::vector<float> genparticles_px;
      std::vector<float> genparticles_py;
      std::vector<float> genparticles_pz;
      std::vector<float> genparticles_mass;

      for (auto const& part:genparts){
        auto const& extras = part->extras;

        // Anything that has status 21 (incoming), 22 (intermediate), 23 (outgoing hard process), 1 (outgoing without treatment except momentum balance) and 2 (only in hard process, like 1 but with subsequent Pythia decay, like in taus)
        // Note that momentum balance among isHardProcess particles might not be exact since status==1 particles have adjustment
        if (
          !(
            extras.isHardProcess
            ||
            extras.status==1
            )
          ) continue;

        genparticles_isPromptFinalState.push_back(extras.isPromptFinalState);
        genparticles_isPromptDecayed.push_back(extras.isPromptDecayed);
        genparticles_isDirectPromptTauDecayProductFinalState.push_back(extras.isDirectPromptTauDecayProductFinalState);
        genparticles_isHardProcess.push_back(extras.isHardProcess);
        genparticles_fromHardProcessFinalState.push_back(extras.fromHardProcessFinalState);
        genparticles_fromHardProcessDecayed.push_back(extras.fromHardProcessDecayed);
        genparticles_isDirectHardProcessTauDecayProductFinalState.push_back(extras.isDirectHardProcessTauDecayProductFinalState);
        genparticles_fromHardProcessBeforeFSR.push_back(extras.fromHardProcessBeforeFSR);
        genparticles_isLastCopy.push_back(extras.isLastCopy);
        genparticles_isLastCopyBeforeFSR.push_back(extras.isLastCopyBeforeFSR);

        genparticles_status.push_back(extras.status);

        genparticles_id.push_back(part->id);
        genparticles_px.push_back(part->x());
        genparticles_py.push_back(part->y());
        genparticles_pz.push_back(part->z());
        genparticles_mass.push_back(part->m());
      }

      product.setNamedVal("genparticles_isPromptFinalState", genparticles_isPromptFinalState);
      product.setNamedVal("genparticles_isPromptDecayed", genparticles_isPromptDecayed);
      product.setNamedVal("genparticles_isDirectPromptTauDecayProductFinalState", genparticles_isDirectPromptTauDecayProductFinalState);
      product.setNamedVal("genparticles_isHardProcess", genparticles_isHardProcess);
      product.setNamedVal("genparticles_fromHardProcessFinalState", genparticles_fromHardProcessFinalState);
      product.setNamedVal("genparticles_fromHardProcessDecayed", genparticles_fromHardProcessDecayed);
      product.setNamedVal("genparticles_isDirectHardProcessTauDecayProductFinalState", genparticles_isDirectHardProcessTauDecayProductFinalState);
      product.setNamedVal("genparticles_fromHardProcessBeforeFSR", genparticles_fromHardProcessBeforeFSR);
      product.setNamedVal("genparticles_isLastCopy", genparticles_isLastCopy);
      product.setNamedVal("genparticles_isLastCopyBeforeFSR", genparticles_isLastCopyBeforeFSR);

      product.setNamedVal("genparticles_id", genparticles_id);
      product.setNamedVal("genparticles_status", genparticles_status);

      product.setNamedVal("genparticles_px", genparticles_px);
      product.setNamedVal("genparticles_py", genparticles_py);
      product.setNamedVal("genparticles_pz", genparticles_pz);
      product.setNamedVal("genparticles_mass", genparticles_mass);
    }
  }


  /********************/
  /**  EVENT FILTER  **/
  /********************/
  validProducts &= (!doEventFilter || eventFilter!=nullptr);
  if (!validProducts){
    MELAerr << "GenericEventAnalyzer::runEvent: Event filter handle is invalid (Tree: " << tree->sampleIdentifier << ")." << endl;
    return validProducts;
  }
  if (eventFilter){
    bool passEventFilters=true;
    validProducts &= eventFilter->constructFilter();
    if (!validProducts){
      MELAerr << "GenericEventAnalyzer::runEvent: Event filter HLT paths are empty. (Tree: " << tree->sampleIdentifier << ")." << endl;
      return validProducts;
    }
    passEventFilters = eventFilter->passEventFilters();
    std::unordered_map<TString, bool> const* passingHLTPaths = &(eventFilter->getHLTProduct());
    for (auto const& pair:(*passingHLTPaths)) product.setNamedVal<bool>((TString("passHLTPath_")+pair.first), pair.second);
    product.setNamedVal<bool>("passEventFilters", passEventFilters);
  }


  /*********************/
  /**  PFCANDS BLOCK  **/
  /*********************/
  std::vector<PFCandObject*> pfcands;
  validProducts &= (!doPFCands || pfcandHandler!=nullptr);
  if (!validProducts){
    MELAerr << "GenericEventAnalyzer::runEvent: PFCand handle is invalid (Tree: " << tree->sampleIdentifier << ")." << endl;
    return validProducts;
  }
  if (pfcandHandler){
    validProducts &= pfcandHandler->constructPFCands();
    if (!validProducts){
      MELAerr << "GenericEventAnalyzer::runEvent: PFCands could not be constructed (Tree: " << tree->sampleIdentifier << ")." << endl;
      tree->print();
      exit(1);
      return validProducts;
    }
    pfcands = pfcandHandler->getProducts();

    std::vector<int> id;
    std::vector<int> charge;

    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> mass;

    std::vector<bool> trackHighPurity;

    std::vector<float> dxy;
    std::vector<float> dz;
    std::vector<float> dxyError;
    std::vector<float> dzError;

    for (PFCandObject const* pfcand:pfcands){
      if (!pfcand) continue;

      PFCandVariables const& extras = pfcand->extras;

      id.push_back(pfcand->id);
      charge.push_back(extras.charge);

      pt.push_back(pfcand->pt());
      eta.push_back(pfcand->eta());
      phi.push_back(pfcand->phi());
      mass.push_back(pfcand->m());

      trackHighPurity.push_back(extras.trackHighPurity);
      dxy.push_back(extras.dxy);
      dz.push_back(extras.dz);
      dxyError.push_back(extras.dxyError);
      dzError.push_back(extras.dzError);
    }
    product.setNamedVal<std::vector<int>>("pfcands_id", id);
    product.setNamedVal<std::vector<int>>("pfcands_charge", charge);

    product.setNamedVal<std::vector<float>>("pfcands_pt", pt);
    product.setNamedVal<std::vector<float>>("pfcands_eta", eta);
    product.setNamedVal<std::vector<float>>("pfcands_phi", phi);
    product.setNamedVal<std::vector<float>>("pfcands_mass", mass);

    product.setNamedVal<std::vector<bool>>("pfcands_trackHighPurity", trackHighPurity);

    product.setNamedVal<std::vector<float>>("pfcands_dxy", dxy);
    product.setNamedVal<std::vector<float>>("pfcands_dz", dz);
    product.setNamedVal<std::vector<float>>("pfcands_dxyError", dxyError);
    product.setNamedVal<std::vector<float>>("pfcands_dzError", dzError);
  }


  /*******************/
  /**  MUONS BLOCK  **/
  /*******************/
  std::vector<MuonObject*> muons;
  validProducts &= (!doMuons || muonHandler!=nullptr);
  if (!validProducts){
    MELAerr << "GenericEventAnalyzer::runEvent: Muon handle is invalid (Tree: " << tree->sampleIdentifier << ")." << endl;
    return validProducts;
  }
  if (muonHandler){
    validProducts &= muonHandler->constructMuons();
    if (!validProducts){
      MELAerr << "GenericEventAnalyzer::runEvent: Muons could not be constructed (Tree: " << tree->sampleIdentifier << ")." << endl;
      tree->print();
      exit(1);
      return validProducts;
    }
    muons = muonHandler->getProducts();

    std::vector<int> id;
    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> mass;

    std::vector<bool> isPFMuon;

    std::vector<long long> POGSelectorBit;

    std::vector<int> type;
    std::vector<int> validHits;
    std::vector<int> lostHits;
    std::vector<int> expectedMissingInnerHits;
    std::vector<int> expectedMissingOuterHits;
    std::vector<int> GlobalFit_Ndof;

    std::vector<float> GlobalFit_Chisq;
    std::vector<float> LocalPos_Chisq;
    std::vector<float> TrkKink;
    std::vector<float> SegComp;
    std::vector<float> dxyPV;
    std::vector<float> dzPV;
    std::vector<float> miniIso_ch;
    std::vector<float> miniIso_nh;
    std::vector<float> miniIso_em;

    std::vector<long long> selectionBits;

    float SF_IdIso=1;
    float SF_Reco=1;
    float SF_Gen=1;
    float SF_IdIso_Up=1;
    float SF_Reco_Up=1;
    float SF_Gen_Up=1;
    float SF_IdIso_Dn=1;
    float SF_Reco_Dn=1;
    float SF_Gen_Dn=1;

    for (MuonObject const* muon:muons){
      if (!muon) continue;
      if (!HelperFunctions::test_bit(muon->selectionBits, MuonSelectionHelpers::kGenPtEta)) continue;

      MuonVariables const& extras = muon->extras;

      id.push_back(muon->id);
      pt.push_back(muon->pt());
      eta.push_back(muon->eta());
      phi.push_back(muon->phi());
      mass.push_back(muon->m());

      isPFMuon.push_back(extras.isPFMuon);

      POGSelectorBit.push_back(extras.POGSelectorBit);

      type.push_back(extras.type);
      validHits.push_back(extras.validHits);
      lostHits.push_back(extras.lostHits);
      expectedMissingInnerHits.push_back(extras.expectedMissingInnerHits);
      expectedMissingOuterHits.push_back(extras.expectedMissingOuterHits);
      GlobalFit_Ndof.push_back(extras.GlobalFit_Ndof);

      GlobalFit_Chisq.push_back(extras.GlobalFit_Chisq);
      LocalPos_Chisq.push_back(extras.LocalPos_Chisq);
      TrkKink.push_back(extras.TrkKink);
      SegComp.push_back(extras.SegComp);
      dxyPV.push_back(extras.dxyPV);
      dzPV.push_back(extras.dzPV);
      miniIso_ch.push_back(extras.miniIso_ch);
      miniIso_nh.push_back(extras.miniIso_nh);
      miniIso_em.push_back(extras.miniIso_em);

      selectionBits.push_back(muon->selectionBits);

      if (muonSFHandler){
        float tmp_SF_IdIso=1;
        float tmp_SF_Reco=1;
        float tmp_SF_Gen=1;
        float tmp_SFerr_IdIso=0;
        float tmp_SFerr_Reco=0;
        float tmp_SFerr_Gen=0;

        muonSFHandler->getIdIsoSFAndError(tmp_SF_IdIso, tmp_SFerr_IdIso, muon, tree->isFastSim());
        muonSFHandler->getRecoSFAndError(tmp_SF_Reco, tmp_SFerr_Reco, muon);
        muonSFHandler->getGenSFAndError(tmp_SF_Gen, tmp_SFerr_Gen, muon, tmp_SF_IdIso, tmp_SFerr_IdIso);

        if (!(isfinite(tmp_SF_IdIso) && isfinite(tmp_SF_Reco) && isfinite(tmp_SF_Gen) && isfinite(tmp_SFerr_IdIso) && isfinite(tmp_SFerr_Reco) && isfinite(tmp_SFerr_Gen))){
          if (verbosity>=TVar::ERROR){
            MELAerr << "GenericEventAnalyzer::runEvent: Some muon scale factors are not finite!" << endl;
            MELAerr << "\t- ID+iso = " << tmp_SF_IdIso << " * (1 +- " << tmp_SFerr_IdIso << ")" << endl;
            MELAerr << "\t- ID+iso = " << tmp_SF_Reco << " * (1 +- " << tmp_SFerr_Reco << ")" << endl;
            MELAerr << "\t- ID+iso = " << tmp_SF_Gen << " * (1 +- " << tmp_SFerr_Gen << ")" << endl;
          }
          exit(1);
        }
        SF_IdIso *= std::min(1.f, std::max(0.f, tmp_SF_IdIso));
        SF_Reco *= std::min(1.f, std::max(0.f, tmp_SF_Reco));
        SF_Gen *= std::min(1.f, std::max(0.f, tmp_SF_Gen));
        SF_IdIso_Up *= std::min(1.f, std::max(0.f, tmp_SF_IdIso * (1.f + tmp_SFerr_IdIso)));
        SF_Reco_Up *= std::min(1.f, std::max(0.f, tmp_SF_Reco * (1.f + tmp_SFerr_Reco)));
        SF_Gen_Up *= std::min(1.f, std::max(0.f, tmp_SF_Gen * (1.f + tmp_SFerr_Gen)));
        SF_IdIso_Dn *= std::min(1.f, std::max(0.f, tmp_SF_IdIso * (1.f - tmp_SFerr_IdIso)));
        SF_Reco_Dn *= std::min(1.f, std::max(0.f, tmp_SF_Reco * (1.f - tmp_SFerr_Reco)));
        SF_Gen_Dn *= std::min(1.f, std::max(0.f, tmp_SF_Gen * (1.f - tmp_SFerr_Gen)));
      }
    }
    product.setNamedVal<std::vector<int>>("muons_id", id);
    product.setNamedVal<std::vector<float>>("muons_pt", pt);
    product.setNamedVal<std::vector<float>>("muons_eta", eta);
    product.setNamedVal<std::vector<float>>("muons_phi", phi);
    product.setNamedVal<std::vector<float>>("muons_mass", mass);
    if (doWriteSelectionVariables){
      product.setNamedVal<std::vector<bool>>("muons_isPFMuon", isPFMuon);

      product.setNamedVal<std::vector<long long>>("muons_POGSelectorBit", POGSelectorBit);

      product.setNamedVal<std::vector<int>>("muons_type", type);
      product.setNamedVal<std::vector<int>>("muons_validHits", validHits);
      product.setNamedVal<std::vector<int>>("muons_lostHits", lostHits);
      product.setNamedVal<std::vector<int>>("muons_expectedMissingInnerHits", expectedMissingInnerHits);
      product.setNamedVal<std::vector<int>>("muons_expectedMissingOuterHits", expectedMissingOuterHits);
      product.setNamedVal<std::vector<int>>("muons_GlobalFit_Ndof", GlobalFit_Ndof);

      product.setNamedVal<std::vector<float>>("muons_GlobalFit_Chisq", GlobalFit_Chisq);
      product.setNamedVal<std::vector<float>>("muons_LocalPos_Chisq", LocalPos_Chisq);
      product.setNamedVal<std::vector<float>>("muons_TrkKink", TrkKink);
      product.setNamedVal<std::vector<float>>("muons_SegComp", SegComp);
      product.setNamedVal<std::vector<float>>("muons_dxyPV", dxyPV);
      product.setNamedVal<std::vector<float>>("muons_dzPV", dzPV);
      product.setNamedVal<std::vector<float>>("muons_miniIso_ch", miniIso_ch);
      product.setNamedVal<std::vector<float>>("muons_miniIso_nh", miniIso_nh);
      product.setNamedVal<std::vector<float>>("muons_miniIso_em", miniIso_em);
    }
    product.setNamedVal<std::vector<long long>>("muons_selectionBits", selectionBits);

    float muonSF = SF_IdIso*SF_Reco;
    float muonSFUp = SF_IdIso_Up*SF_Reco_Up;
    float muonSFDn = SF_IdIso_Dn*SF_Reco_Dn;
    product.setNamedVal<float>("weight_muons", muonSF);
    product.setNamedVal<float>("weight_muons_SFUp", muonSFUp);
    product.setNamedVal<float>("weight_muons_SFDn", muonSFDn);
  }


  /***********************/
  /**  ELECTRONS BLOCK  **/
  /***********************/
  std::vector<ElectronObject*> electrons;
  validProducts &= (!doElectrons || eleHandler!=nullptr);
  if (!validProducts){
    MELAerr << "GenericEventAnalyzer::runEvent: Electron handle is invalid (Tree: " << tree->sampleIdentifier << ")." << endl;
    return validProducts;
  }
  if (eleHandler){
    validProducts &= eleHandler->constructElectrons();
    if (!validProducts){
      MELAerr << "GenericEventAnalyzer::runEvent: Electrons could not be constructed (Tree: " << tree->sampleIdentifier << ")." << endl;
      tree->print();
      exit(1);
      return validProducts;
    }
    electrons = eleHandler->getProducts();

    std::vector<int> id;
    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> mass;

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
    std::vector<long long> selectionBits;

    float SF_IdIso=1;
    float SF_Reco=1;
    float SF_Gen=1;
    float SF_IdIso_Up=1;
    float SF_Reco_Up=1;
    float SF_Gen_Up=1;
    float SF_IdIso_Dn=1;
    float SF_Reco_Dn=1;
    float SF_Gen_Dn=1;

    for (ElectronObject const* electron:electrons){
      if (!electron) continue;
      if (!HelperFunctions::test_bit(electron->selectionBits, ElectronSelectionHelpers::kGenPtEta)) continue;

      ElectronVariables const& extras = electron->extras;

      id.push_back(electron->id);
      pt.push_back(electron->pt());
      eta.push_back(electron->eta());
      phi.push_back(electron->phi());
      mass.push_back(electron->m());

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

      selectionBits.push_back(electron->selectionBits);

      if (eleSFHandler){
        float tmp_SF_IdIso=1;
        float tmp_SF_Reco=1;
        float tmp_SF_Gen=1;
        float tmp_SFerr_IdIso=0;
        float tmp_SFerr_Reco=0;
        float tmp_SFerr_Gen=0;

        eleSFHandler->getIdIsoSFAndError(tmp_SF_IdIso, tmp_SFerr_IdIso, electron, tree->isFastSim());
        eleSFHandler->getRecoSFAndError(tmp_SF_Reco, tmp_SFerr_Reco, electron);
        eleSFHandler->getGenSFAndError(tmp_SF_Gen, tmp_SFerr_Gen, electron, tmp_SF_IdIso, tmp_SFerr_IdIso);

        if (!(isfinite(tmp_SF_IdIso) && isfinite(tmp_SF_Reco) && isfinite(tmp_SF_Gen) && isfinite(tmp_SFerr_IdIso) && isfinite(tmp_SFerr_Reco) && isfinite(tmp_SFerr_Gen))){
          if (verbosity>=TVar::ERROR){
            MELAerr << "GenericEventAnalyzer::runEvent: Some electron scale factors are not finite!" << endl;
            MELAerr << "\t- ID+iso = " << tmp_SF_IdIso << " * (1 +- " << tmp_SFerr_IdIso << ")" << endl;
            MELAerr << "\t- ID+iso = " << tmp_SF_Reco << " * (1 +- " << tmp_SFerr_Reco << ")" << endl;
            MELAerr << "\t- ID+iso = " << tmp_SF_Gen << " * (1 +- " << tmp_SFerr_Gen << ")" << endl;

          }
          exit(1);

        }
        SF_IdIso *= std::min(1.f, std::max(0.f, tmp_SF_IdIso));
        SF_Reco *= std::min(1.f, std::max(0.f, tmp_SF_Reco));
        SF_Gen *= std::min(1.f, std::max(0.f, tmp_SF_Gen));
        SF_IdIso_Up *= std::min(1.f, std::max(0.f, tmp_SF_IdIso * (1.f + tmp_SFerr_IdIso)));
        SF_Reco_Up *= std::min(1.f, std::max(0.f, tmp_SF_Reco * (1.f + tmp_SFerr_Reco)));
        SF_Gen_Up *= std::min(1.f, std::max(0.f, tmp_SF_Gen * (1.f + tmp_SFerr_Gen)));
        SF_IdIso_Dn *= std::min(1.f, std::max(0.f, tmp_SF_IdIso * (1.f - tmp_SFerr_IdIso)));
        SF_Reco_Dn *= std::min(1.f, std::max(0.f, tmp_SF_Reco * (1.f - tmp_SFerr_Reco)));
        SF_Gen_Dn *= std::min(1.f, std::max(0.f, tmp_SF_Gen * (1.f - tmp_SFerr_Gen)));
      }
    }
    product.setNamedVal<std::vector<int>>("electrons_id", id);
    product.setNamedVal<std::vector<float>>("electrons_pt", pt);
    product.setNamedVal<std::vector<float>>("electrons_eta", eta);
    product.setNamedVal<std::vector<float>>("electrons_phi", phi);
    product.setNamedVal<std::vector<float>>("electrons_mass", mass);
    if (doWriteSelectionVariables){
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
    }
    product.setNamedVal<std::vector<long long>>("electrons_selectionBits", selectionBits);

    float electronSF = SF_IdIso*SF_Reco;
    float electronSFUp = SF_IdIso_Up*SF_Reco_Up;
    float electronSFDn = SF_IdIso_Dn*SF_Reco_Dn;
    product.setNamedVal<float>("weight_electrons", electronSF);
    product.setNamedVal<float>("weight_electrons_SFUp", electronSFUp);
    product.setNamedVal<float>("weight_electrons_SFDn", electronSFDn);
  }


  /*********************/
  /**  PHOTONS BLOCK  **/
  /*********************/
  std::vector<PhotonObject*> photons;
  validProducts &= (!doPhotons || photonHandler!=nullptr);
  if (!validProducts){
    MELAerr << "GenericEventAnalyzer::runEvent: Photon handle is invalid (Tree: " << tree->sampleIdentifier << ")." << endl;
    return validProducts;
  }
  if (photonHandler){
    validProducts &= photonHandler->constructPhotons();
    if (!validProducts){
      MELAerr << "GenericEventAnalyzer::runEvent: Photons could not be constructed (Tree: " << tree->sampleIdentifier << ")." << endl;
      tree->print();
      exit(1);
      return validProducts;
    }
    photons = photonHandler->getProducts();

    std::vector<float> pt;
    std::vector<float> eta;
    std::vector<float> phi;
    std::vector<float> mass;

    std::vector<float> etaSC;
    std::vector<float> recoChargedHadronIso;
    std::vector<float> recoNeutralHadronIso;
    std::vector<float> recoPhotonIso;
    std::vector<float> sigmaIEtaIEta_full5x5;
    std::vector<float> hOverE_full5x5;
    std::vector<long long> selectionBits;

    for (PhotonObject const* photon:photons){
      if (!photon) continue;
      if (!HelperFunctions::test_bit(photon->selectionBits, PhotonSelectionHelpers::kGenPtEta)) continue;

      PhotonVariables const& extras = photon->extras;

      pt.push_back(photon->pt());
      eta.push_back(photon->eta());
      phi.push_back(photon->phi());
      mass.push_back(photon->m());

      etaSC.push_back(extras.etaSC);
      recoChargedHadronIso.push_back(extras.recoChargedHadronIso);
      recoNeutralHadronIso.push_back(extras.recoNeutralHadronIso);
      recoPhotonIso.push_back(extras.recoPhotonIso);
      sigmaIEtaIEta_full5x5.push_back(extras.sigmaIEtaIEta_full5x5);
      hOverE_full5x5.push_back(extras.hOverE_full5x5);

      selectionBits.push_back(photon->selectionBits);
    }
    product.setNamedVal<std::vector<float>>("photons_pt", pt);
    product.setNamedVal<std::vector<float>>("photons_eta", eta);
    product.setNamedVal<std::vector<float>>("photons_phi", phi);
    product.setNamedVal<std::vector<float>>("photons_mass", mass);
    if (doWriteSelectionVariables){
      product.setNamedVal<std::vector<float>>("photons_etaSC", etaSC);
      product.setNamedVal<std::vector<float>>("photons_recoChargedHadronIso", recoChargedHadronIso);
      product.setNamedVal<std::vector<float>>("photons_recoNeutralHadronIso", recoNeutralHadronIso);
      product.setNamedVal<std::vector<float>>("photons_recoPhotonIso", recoPhotonIso);
      product.setNamedVal<std::vector<float>>("photons_sigmaIEtaIEta_full5x5", sigmaIEtaIEta_full5x5);
      product.setNamedVal<std::vector<float>>("photons_hOverE_full5x5", hOverE_full5x5);
    }
    product.setNamedVal<std::vector<long long>>("photons_selectionBits", selectionBits);
  }


  /********************/
  /**  JetMET BLOCK  **/
  /********************/
  validProducts &= (!doJetMET || jetHandler!=nullptr);
  if (!validProducts){
    MELAerr << "GenericEventAnalyzer::runEvent: JetMET handle is invalid (Tree: " << tree->sampleIdentifier << ")." << endl;
    return validProducts;
  }
  if (jetHandler){
    jetHandler->registerParticles(&muons, &electrons, &photons);
    validProducts &= jetHandler->constructJetMET();
    if (!validProducts){
      MELAerr << "GenericEventAnalyzer::runEvent: Jets or MET could not be constructed (Tree: " << tree->sampleIdentifier << ")." << endl;
      return validProducts;
    }

    std::vector<GenJetObject*> const& genjets = jetHandler->getGenJets();
    std::vector<AK4JetObject*> const& ak4jets = jetHandler->getAK4Jets();
    std::vector<AK8JetObject*> const& ak8jets = jetHandler->getAK8Jets();
    std::vector<TFTopObject*> const& tftops = jetHandler->getTFTops();
    METObject const* metobject = jetHandler->getMET();
    if (jetHandler->getMETFlag()){
      product.setNamedVal("pfmet_original", metobject->extras.met);
      product.setNamedVal("pfmetPhi_original", metobject->extras.phi);
      product.setNamedVal("pfmet", metobject->extras.met_JEC);
      product.setNamedVal("pfmetPhi", metobject->extras.phi_JEC);
      product.setNamedVal("pfmet_JECup", metobject->extras.met_JECup);
      product.setNamedVal("pfmetPhi_JECup", metobject->extras.phi_JECup);
      product.setNamedVal("pfmet_JECdn", metobject->extras.met_JECdn);
      product.setNamedVal("pfmetPhi_JECdn", metobject->extras.phi_JECdn);
    }

    // GenJets
    std::vector<float> genjets_pt;
    std::vector<float> genjets_eta;
    std::vector<float> genjets_phi;
    std::vector<float> genjets_mass;
    for (GenJetObject const* jet:genjets){
      if (!jet) continue;
      CMSLorentzVector const& momentum = jet->momentum;
      genjets_pt.push_back(momentum.Pt());
      genjets_eta.push_back(momentum.Eta());
      genjets_phi.push_back(momentum.Phi());
      genjets_mass.push_back(momentum.M());
    }
    if (jetHandler->getGenJetsFlag()){
      product.setNamedVal("genjets_pt", genjets_pt);
      product.setNamedVal("genjets_eta", genjets_eta);
      product.setNamedVal("genjets_phi", genjets_phi);
      product.setNamedVal("genjets_mass", genjets_mass);
    }

    // AK4Jets
    std::vector<float> ak4jets_pt;
    std::vector<float> ak4jets_eta;
    std::vector<float> ak4jets_phi;
    std::vector<float> ak4jets_mass;

    std::vector<int> ak4jets_npfcands;
    std::vector<int> ak4jets_parton_flavor;
    std::vector<int> ak4jets_hadron_flavor;
    std::vector<int> ak4jets_chargedHadronMultiplicity;
    std::vector<int> ak4jets_neutralHadronMultiplicity;
    std::vector<int> ak4jets_photonMultiplicity;
    std::vector<int> ak4jets_electronMultiplicity;
    std::vector<int> ak4jets_muonMultiplicity;
    std::vector<int> ak4jets_chargedMultiplicity;
    std::vector<int> ak4jets_neutralMultiplicity;
    std::vector<int> ak4jets_totalMultiplicity;

    std::vector<float> ak4jets_area;
    std::vector<float> ak4jets_undoJEC;
    std::vector<float> ak4jets_chargedHadronE;
    std::vector<float> ak4jets_chargedEmE;
    std::vector<float> ak4jets_neutralHadronE;
    std::vector<float> ak4jets_neutralEmE;
    std::vector<float> ak4jets_hfHadronE;
    std::vector<float> ak4jets_hfEmE;
    std::vector<float> ak4jets_photonE;
    std::vector<float> ak4jets_electronE;
    std::vector<float> ak4jets_muonE;

    std::vector<float> ak4jets_deepCSVb;
    std::vector<float> ak4jets_deepCSVc;
    std::vector<float> ak4jets_deepCSVl;
    std::vector<float> ak4jets_deepCSVbb;
    std::vector<float> ak4jets_deepCSVcc;
    std::vector<float> ak4jets_pfCombinedInclusiveSecondaryVertexV2BJetTag;
    std::vector<float> ak4jets_ptDistribution;
    std::vector<float> ak4jets_axis1;
    std::vector<float> ak4jets_axis2;

    std::vector<float> ak4jets_JEC;
    std::vector<float> ak4jets_JECup;
    std::vector<float> ak4jets_JECdn;

    std::vector<float> ak4jets_estimatedPtResolution;
    std::vector<float> ak4jets_JER;
    std::vector<float> ak4jets_JERup;
    std::vector<float> ak4jets_JERdn;

    std::vector<long long> ak4jets_selectionBits;

    std::vector<int> ak4jets_genjetIndex;
    std::vector<float> ak4jets_genjetDeltaR;

    for (AK4JetObject const* jet:ak4jets){
      if (!jet) continue;
      auto const& extras = jet->extras;

      ak4jets_selectionBits.push_back(jet->selectionBits);

      CMSLorentzVector finalMomentum = jet->getFinalMomentum();
      ak4jets_pt.push_back(finalMomentum.Pt());
      ak4jets_eta.push_back(finalMomentum.Eta());
      ak4jets_phi.push_back(finalMomentum.Phi());
      ak4jets_mass.push_back(finalMomentum.M());

      ak4jets_npfcands.push_back(extras.npfcands);
      ak4jets_parton_flavor.push_back(extras.parton_flavor);
      ak4jets_hadron_flavor.push_back(extras.hadron_flavor);
      ak4jets_chargedHadronMultiplicity.push_back(extras.chargedHadronMultiplicity);
      ak4jets_neutralHadronMultiplicity.push_back(extras.neutralHadronMultiplicity);
      ak4jets_photonMultiplicity.push_back(extras.photonMultiplicity);
      ak4jets_electronMultiplicity.push_back(extras.electronMultiplicity);
      ak4jets_muonMultiplicity.push_back(extras.muonMultiplicity);
      ak4jets_chargedMultiplicity.push_back(extras.chargedMultiplicity);
      ak4jets_neutralMultiplicity.push_back(extras.neutralMultiplicity);
      ak4jets_totalMultiplicity.push_back(extras.totalMultiplicity);

      ak4jets_area.push_back(extras.area);
      ak4jets_undoJEC.push_back(extras.undoJEC);
      ak4jets_chargedHadronE.push_back(extras.chargedHadronE);
      ak4jets_chargedEmE.push_back(extras.chargedEmE);
      ak4jets_neutralHadronE.push_back(extras.neutralHadronE);
      ak4jets_neutralEmE.push_back(extras.neutralEmE);
      ak4jets_hfHadronE.push_back(extras.hfHadronE);
      ak4jets_hfEmE.push_back(extras.hfEmE);
      ak4jets_photonE.push_back(extras.photonE);
      ak4jets_electronE.push_back(extras.electronE);
      ak4jets_muonE.push_back(extras.muonE);

      ak4jets_deepCSVb.push_back(extras.deepCSVb);
      ak4jets_deepCSVc.push_back(extras.deepCSVc);
      ak4jets_deepCSVl.push_back(extras.deepCSVl);
      ak4jets_deepCSVbb.push_back(extras.deepCSVbb);
      ak4jets_deepCSVcc.push_back(extras.deepCSVcc);
      ak4jets_pfCombinedInclusiveSecondaryVertexV2BJetTag.push_back(extras.pfCombinedInclusiveSecondaryVertexV2BJetTag);
      ak4jets_ptDistribution.push_back(extras.ptDistribution);
      ak4jets_axis1.push_back(extras.axis1);
      ak4jets_axis2.push_back(extras.axis2);

      ak4jets_estimatedPtResolution.push_back(extras.estimatedPtResolution);

      ak4jets_JEC.push_back((extras.JEC==extras.JECup && extras.JEC==extras.JECdn ? 1.f : extras.JEC/extras.undoJEC));
      ak4jets_JER.push_back(extras.JER);
      ak4jets_JECup.push_back(extras.JECup/extras.JEC);
      ak4jets_JECdn.push_back(extras.JECdn/extras.JEC);
      ak4jets_JERup.push_back(extras.JERup/extras.JER);
      ak4jets_JERdn.push_back(extras.JERdn/extras.JER);

      int genjetIndex=-1;
      float genjetDeltaR=-1;
      {
        unsigned int igenjet=0;
        for (GenJetObject const* tmp:genjets){
          if (tmp==jet->associatedGenJet){
            genjetIndex = igenjet;
            genjetDeltaR = reco::deltaR(finalMomentum, tmp->momentum);
            break;
          }
          igenjet++;
        }
      }
      ak4jets_genjetIndex.push_back(genjetIndex);
      ak4jets_genjetDeltaR.push_back(genjetDeltaR);
    }
    if (jetHandler->getAK4JetsFlag()){
      product.setNamedVal("ak4jets_pt", ak4jets_pt);
      product.setNamedVal("ak4jets_eta", ak4jets_eta);
      product.setNamedVal("ak4jets_phi", ak4jets_phi);
      product.setNamedVal("ak4jets_mass", ak4jets_mass);

      product.setNamedVal("ak4jets_npfcands", ak4jets_npfcands);
      product.setNamedVal("ak4jets_parton_flavor", ak4jets_parton_flavor);
      product.setNamedVal("ak4jets_hadron_flavor", ak4jets_hadron_flavor);
      product.setNamedVal("ak4jets_deepCSVb", ak4jets_deepCSVb);
      product.setNamedVal("ak4jets_deepCSVc", ak4jets_deepCSVc);
      product.setNamedVal("ak4jets_deepCSVl", ak4jets_deepCSVl);
      product.setNamedVal("ak4jets_deepCSVbb", ak4jets_deepCSVbb);
      product.setNamedVal("ak4jets_deepCSVcc", ak4jets_deepCSVcc);
      product.setNamedVal("ak4jets_pfCombinedInclusiveSecondaryVertexV2BJetTag", ak4jets_pfCombinedInclusiveSecondaryVertexV2BJetTag);
      product.setNamedVal("ak4jets_ptDistribution", ak4jets_ptDistribution);
      product.setNamedVal("ak4jets_axis1", ak4jets_axis1);
      product.setNamedVal("ak4jets_axis2", ak4jets_axis2);
      if (doWriteSelectionVariables){
        product.setNamedVal("ak4jets_chargedHadronMultiplicity", ak4jets_chargedHadronMultiplicity);
        product.setNamedVal("ak4jets_neutralHadronMultiplicity", ak4jets_neutralHadronMultiplicity);
        product.setNamedVal("ak4jets_photonMultiplicity", ak4jets_photonMultiplicity);
        product.setNamedVal("ak4jets_electronMultiplicity", ak4jets_electronMultiplicity);
        product.setNamedVal("ak4jets_muonMultiplicity", ak4jets_muonMultiplicity);
        product.setNamedVal("ak4jets_chargedMultiplicity", ak4jets_chargedMultiplicity);
        product.setNamedVal("ak4jets_neutralMultiplicity", ak4jets_neutralMultiplicity);
        product.setNamedVal("ak4jets_totalMultiplicity", ak4jets_totalMultiplicity);

        product.setNamedVal("ak4jets_area", ak4jets_area);
        product.setNamedVal("ak4jets_chargedHadronE", ak4jets_chargedHadronE);
        product.setNamedVal("ak4jets_chargedEmE", ak4jets_chargedEmE);
        product.setNamedVal("ak4jets_neutralHadronE", ak4jets_neutralHadronE);
        product.setNamedVal("ak4jets_neutralEmE", ak4jets_neutralEmE);
        product.setNamedVal("ak4jets_hfHadronE", ak4jets_hfHadronE);
        product.setNamedVal("ak4jets_hfEmE", ak4jets_hfEmE);
        product.setNamedVal("ak4jets_photonE", ak4jets_photonE);
        product.setNamedVal("ak4jets_electronE", ak4jets_electronE);
        product.setNamedVal("ak4jets_muonE", ak4jets_muonE);
      }
      product.setNamedVal("ak4jets_undoJEC", ak4jets_undoJEC);
      product.setNamedVal("ak4jets_JEC", ak4jets_JEC);
      product.setNamedVal("ak4jets_JECup", ak4jets_JECup);
      product.setNamedVal("ak4jets_JECdn", ak4jets_JECdn);

      product.setNamedVal("ak4jets_estimatedPtResolution", ak4jets_estimatedPtResolution);
      product.setNamedVal("ak4jets_JER", ak4jets_JER);
      product.setNamedVal("ak4jets_JERup", ak4jets_JERup);
      product.setNamedVal("ak4jets_JERdn", ak4jets_JERdn);

      product.setNamedVal("ak4jets_selectionBits", ak4jets_selectionBits);

      product.setNamedVal("ak4jets_genjetIndex", ak4jets_genjetIndex);
      product.setNamedVal("ak4jets_genjetDeltaR", ak4jets_genjetDeltaR);
    }

    // AK8Jets
    std::vector<float> ak8jets_pt;
    std::vector<float> ak8jets_eta;
    std::vector<float> ak8jets_phi;
    std::vector<float> ak8jets_mass;

    std::vector<int> ak8jets_parton_flavor;

    std::vector<float> ak8jets_area;
    std::vector<float> ak8jets_undoJEC;
    std::vector<float> ak8jets_tau1;
    std::vector<float> ak8jets_tau2;
    std::vector<float> ak8jets_tau3;
    std::vector<float> ak8jets_deepdisc_qcd;
    std::vector<float> ak8jets_deepdisc_top;
    std::vector<float> ak8jets_deepdisc_w;
    std::vector<float> ak8jets_deepdisc_z;
    std::vector<float> ak8jets_deepdisc_zbb;
    std::vector<float> ak8jets_deepdisc_hbb;
    std::vector<float> ak8jets_deepdisc_h4q;

    std::vector<float> ak8jets_JEC;
    std::vector<float> ak8jets_JECup;
    std::vector<float> ak8jets_JECdn;

    std::vector<float> ak8jets_estimatedPtResolution;
    std::vector<float> ak8jets_JER;
    std::vector<float> ak8jets_JERup;
    std::vector<float> ak8jets_JERdn;

    std::vector<long long> ak8jets_selectionBits;

    std::vector<int> ak8jets_genjetIndex;
    std::vector<float> ak8jets_genjetDeltaR;

    for (AK8JetObject const* jet:ak8jets){
      if (!jet) continue;
      auto const& extras = jet->extras;

      ak8jets_selectionBits.push_back(jet->selectionBits);

      CMSLorentzVector finalMomentum = jet->getFinalMomentum();
      ak8jets_pt.push_back(finalMomentum.Pt());
      ak8jets_eta.push_back(finalMomentum.Eta());
      ak8jets_phi.push_back(finalMomentum.Phi());
      ak8jets_mass.push_back(finalMomentum.M());

      ak8jets_parton_flavor.push_back(extras.parton_flavor);

      ak8jets_area.push_back(extras.area);
      ak8jets_undoJEC.push_back(extras.undoJEC);
      ak8jets_tau1.push_back(extras.tau1);
      ak8jets_tau2.push_back(extras.tau2);
      ak8jets_tau3.push_back(extras.tau3);
      ak8jets_deepdisc_qcd.push_back(extras.deepdisc_qcd);
      ak8jets_deepdisc_top.push_back(extras.deepdisc_top);
      ak8jets_deepdisc_w.push_back(extras.deepdisc_w);
      ak8jets_deepdisc_z.push_back(extras.deepdisc_z);
      ak8jets_deepdisc_zbb.push_back(extras.deepdisc_zbb);
      ak8jets_deepdisc_hbb.push_back(extras.deepdisc_hbb);
      ak8jets_deepdisc_h4q.push_back(extras.deepdisc_h4q);

      ak8jets_estimatedPtResolution.push_back(extras.estimatedPtResolution);

      ak8jets_JEC.push_back((extras.JEC==extras.JECup && extras.JEC==extras.JECdn ? 1.f : extras.JEC/extras.undoJEC));
      ak8jets_JER.push_back(extras.JER);
      ak8jets_JECup.push_back(extras.JECup/extras.JEC);
      ak8jets_JECdn.push_back(extras.JECdn/extras.JEC);
      ak8jets_JERup.push_back(extras.JERup/extras.JER);
      ak8jets_JERdn.push_back(extras.JERdn/extras.JER);

      int genjetIndex=-1;
      float genjetDeltaR=-1;
      {
        unsigned int igenjet=0;
        for (GenJetObject const* tmp:genjets){
          if (tmp==jet->associatedGenJet){
            genjetIndex = igenjet;
            genjetDeltaR = reco::deltaR(finalMomentum, tmp->momentum);
            break;
          }
          igenjet++;
        }
      }
      ak8jets_genjetIndex.push_back(genjetIndex);
      ak8jets_genjetDeltaR.push_back(genjetDeltaR);
    }
    if (jetHandler->getAK8JetsFlag()){
      product.setNamedVal("ak8jets_pt", ak8jets_pt);
      product.setNamedVal("ak8jets_eta", ak8jets_eta);
      product.setNamedVal("ak8jets_phi", ak8jets_phi);
      product.setNamedVal("ak8jets_mass", ak8jets_mass);

      product.setNamedVal("ak8jets_parton_flavor", ak8jets_parton_flavor);

      product.setNamedVal("ak8jets_tau1", ak8jets_tau1);
      product.setNamedVal("ak8jets_tau2", ak8jets_tau2);
      product.setNamedVal("ak8jets_tau3", ak8jets_tau3);
      product.setNamedVal("ak8jets_deepdisc_qcd", ak8jets_deepdisc_qcd);
      product.setNamedVal("ak8jets_deepdisc_top", ak8jets_deepdisc_top);
      product.setNamedVal("ak8jets_deepdisc_w", ak8jets_deepdisc_w);
      product.setNamedVal("ak8jets_deepdisc_z", ak8jets_deepdisc_z);
      product.setNamedVal("ak8jets_deepdisc_zbb", ak8jets_deepdisc_zbb);
      product.setNamedVal("ak8jets_deepdisc_hbb", ak8jets_deepdisc_hbb);
      product.setNamedVal("ak8jets_deepdisc_h4q", ak8jets_deepdisc_h4q);
      if (doWriteSelectionVariables){
        product.setNamedVal("ak8jets_area", ak8jets_area);
      }

      product.setNamedVal("ak8jets_undoJEC", ak8jets_undoJEC);
      product.setNamedVal("ak8jets_JEC", ak8jets_JEC);
      product.setNamedVal("ak8jets_JECup", ak8jets_JECup);
      product.setNamedVal("ak8jets_JECdn", ak8jets_JECdn);

      product.setNamedVal("ak8jets_estimatedPtResolution", ak8jets_estimatedPtResolution);
      product.setNamedVal("ak8jets_JER", ak8jets_JER);
      product.setNamedVal("ak8jets_JERup", ak8jets_JERup);
      product.setNamedVal("ak8jets_JERdn", ak8jets_JERdn);

      product.setNamedVal("ak8jets_selectionBits", ak8jets_selectionBits);

      product.setNamedVal("ak8jets_genjetIndex", ak8jets_genjetIndex);
      product.setNamedVal("ak8jets_genjetDeltaR", ak8jets_genjetDeltaR);
    }

    // TFTops
    std::vector<float> tftops_pt;
    std::vector<float> tftops_eta;
    std::vector<float> tftops_phi;
    std::vector<float> tftops_mass;

    std::vector<int> tftops_nSubjets;

    std::vector<float> tftops_disc;

    std::vector<long long> tftops_selectionBits;

    for (TFTopObject const* top:tftops){
      if (!top) continue;
      auto const& extras = top->extras;

      tftops_selectionBits.push_back(top->selectionBits);

      CMSLorentzVector finalMomentum = top->getFinalMomentum();
      tftops_pt.push_back(finalMomentum.Pt());
      tftops_eta.push_back(finalMomentum.Eta());
      tftops_phi.push_back(finalMomentum.Phi());
      tftops_mass.push_back(finalMomentum.M());

      tftops_nSubjets.push_back(extras.nSubjets);

      tftops_disc.push_back(extras.disc);
    }
    if (jetHandler->getTopsFlag()){
      product.setNamedVal("tftops_pt", tftops_pt);
      product.setNamedVal("tftops_eta", tftops_eta);
      product.setNamedVal("tftops_phi", tftops_phi);
      product.setNamedVal("tftops_mass", tftops_mass);

      product.setNamedVal("tftops_nSubjets", tftops_nSubjets);

      product.setNamedVal("tftops_disc", tftops_disc);

      product.setNamedVal("tftops_selectionBits", tftops_selectionBits);
    }
  }

  return validProducts;
}
