#include <cassert>
#include "ParticleObjectHelpers.h"
#include "PhotonHandler.h"
#include "PhotonSelectionHelpers.h"
#include "FrameworkVariables.hh"
#include "FrameworkTree.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


PhotonHandler::PhotonHandler() : IvyBase()
{
  this->addConsumed<float>(_photons_rho_);

#ifdef _photons_etaSC_
  this->addConsumed<std::vector<float>*>(_photons_etaSC_);
#endif
  this->addConsumed<std::vector<float>*>(_photons_recoChargedHadronIso_);
  this->addConsumed<std::vector<float>*>(_photons_recoNeutralHadronIso_);
  this->addConsumed<std::vector<float>*>(_photons_recoPhotonIso_);
  this->addConsumed<std::vector<float>*>(_photons_sigmaIEtaIEta_full5x5_);
  this->addConsumed<std::vector<float>*>(_photons_hOverE_full5x5_);

  this->addConsumed<std::vector<CMSLorentzVector>*>(_photons_momentum_);
}


bool PhotonHandler::constructPhotons(){
  clear();
  if (!currentTree) return false;

  float rho = 0;

#ifdef _photons_etaSC_
  std::vector<float>* etaSC = nullptr;
#endif
  std::vector<float>* recoChargedHadronIso = nullptr;
  std::vector<float>* recoNeutralHadronIso = nullptr;
  std::vector<float>* recoPhotonIso = nullptr;
  std::vector<float>* sigmaIEtaIEta_full5x5 = nullptr;
  std::vector<float>* hOverE_full5x5 = nullptr;

  std::vector<CMSLorentzVector>* momentum = nullptr;

  // Beyond this point starts checks and selection
  bool allVariablesPresent = (
    this->getConsumedValue(_photons_rho_, rho)
#ifdef _photons_etaSC_
    &&
    this->getConsumedValue(_photons_etaSC_, etaSC)
#endif
    &&
    this->getConsumedValue(_photons_recoChargedHadronIso_, recoChargedHadronIso)
    &&
    this->getConsumedValue(_photons_recoNeutralHadronIso_, recoNeutralHadronIso)
    &&
    this->getConsumedValue(_photons_recoPhotonIso_, recoPhotonIso)
    &&
    this->getConsumedValue(_photons_sigmaIEtaIEta_full5x5_, sigmaIEtaIEta_full5x5)
    &&
    this->getConsumedValue(_photons_hOverE_full5x5_, hOverE_full5x5)
    &&
    this->getConsumedValue(_photons_momentum_, momentum)
    );

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "PhotonHandler::constructPhotons: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "PhotonHandler::constructPhotons: All variables are set up!" << endl;
  if (!momentum){
    if (this->verbosity>=TVar::ERROR) MELAerr << "PhotonHandler::constructPhotons: Photons could not be linked! Pointer to " << _photons_momentum_ << " is null!" << endl;
    return false;
  }

  if (momentum->empty()) return true; // Construction is successful, it is just that no photons exist.

  unsigned int nProducts = momentum->size();
  productList.reserve(nProducts);
  for (unsigned int ip=0; ip<nProducts; ip++){
    if (this->verbosity>=TVar::DEBUG) MELAout << "PhotonHandler::constructPhotons: Attempting photon " << ip << "..." << endl;

    productList.push_back(new PhotonObject(22, momentum->at(ip)));
    PhotonObject*& obj = productList.back();

    obj->extras.rho = rho;

#ifdef _photons_etaSC_
    obj->extras.etaSC = etaSC->at(ip);
#else
    obj->extras.etaSC = obj->eta();
#endif
    obj->extras.recoChargedHadronIso = recoChargedHadronIso->at(ip);
    obj->extras.recoNeutralHadronIso = recoNeutralHadronIso->at(ip);
    obj->extras.recoPhotonIso = recoPhotonIso->at(ip);
    obj->extras.sigmaIEtaIEta_full5x5 = sigmaIEtaIEta_full5x5->at(ip);
    obj->extras.hOverE_full5x5 = hOverE_full5x5->at(ip);

    // Set the selection bits
    PhotonSelectionHelpers::setSelectionBits(*obj);

    if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  return true;
}

void PhotonHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(tree);
  if (!fwktree) return;
  SampleHelpers::setupUsingOptions(fwktree->getOptions());

  fwktree->bookEDMBranch<float>(_photons_rho_, 0);

#ifdef _photons_etaSC_
  fwktree->bookEDMBranch<std::vector<float>*>(_photons_etaSC_, nullptr);
#endif
  fwktree->bookEDMBranch<std::vector<float>*>(_photons_recoChargedHadronIso_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_photons_recoNeutralHadronIso_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_photons_recoPhotonIso_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_photons_sigmaIEtaIEta_full5x5_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_photons_hOverE_full5x5_, nullptr);

  fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_photons_momentum_, nullptr);
}

