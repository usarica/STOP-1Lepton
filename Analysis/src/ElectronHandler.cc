#include <cassert>
#include "ElectronHandler.h"
#include "ElectronSelectionHelpers.h"
#include "FrameworkVariables.hh"
#include "FrameworkTree.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


ElectronHandler::ElectronHandler() : IvyBase()
{
  this->addConsumed<float>(_electrons_rho_);

  this->addConsumed<std::vector<bool>*>(_electrons_conv_vtx_flag_);

  this->addConsumed<std::vector<int>*>(_electrons_charge_);
  this->addConsumed<std::vector<int>*>(_electrons_expectedMissingInnerHits_);

  this->addConsumed<std::vector<float>*>(_electrons_energySC_);
  this->addConsumed<std::vector<float>*>(_electrons_etaSC_);
  this->addConsumed<std::vector<float>*>(_electrons_etaSeedSC_);
  this->addConsumed<std::vector<float>*>(_electrons_sigmaIEtaIEta_full5x5_);
  this->addConsumed<std::vector<float>*>(_electrons_dEtaIn_);
  this->addConsumed<std::vector<float>*>(_electrons_dPhiIn_);
  this->addConsumed<std::vector<float>*>(_electrons_hOverE_);
  this->addConsumed<std::vector<float>*>(_electrons_ecalEnergy_);
  this->addConsumed<std::vector<float>*>(_electrons_eOverPIn_);
  this->addConsumed<std::vector<float>*>(_electrons_dxyPV_);
  this->addConsumed<std::vector<float>*>(_electrons_dzPV_);
  this->addConsumed<std::vector<float>*>(_electrons_miniIso_ch_);
  this->addConsumed<std::vector<float>*>(_electrons_miniIso_nh_);
  this->addConsumed<std::vector<float>*>(_electrons_miniIso_em_);

  this->addConsumed<std::vector<CMSLorentzVector>*>(_electrons_momentum_);
}


bool ElectronHandler::constructElectrons(){
  clear();
  if (!currentTree) return false;

  float rho = 0;

  vector<bool>* conv_vtx_flag = nullptr;

  vector<int>* charge = nullptr;
  vector<int>* expectedMissingInnerHits = nullptr;

  vector<float>* energySC = nullptr;
  vector<float>* etaSC = nullptr;
  vector<float>* etaSeedSC = nullptr;
  vector<float>* sigmaIEtaIEta_full5x5 = nullptr;
  vector<float>* dEtaIn = nullptr;
  vector<float>* dPhiIn = nullptr;
  vector<float>* hOverE = nullptr;
  vector<float>* ecalEnergy = nullptr;
  vector<float>* eOverPIn = nullptr;
  vector<float>* dxyPV = nullptr;
  vector<float>* dzPV = nullptr;
  vector<float>* miniIso_ch = nullptr;
  vector<float>* miniIso_nh = nullptr;
  vector<float>* miniIso_em = nullptr;

  vector<CMSLorentzVector>* momentum = nullptr;

  // Beyond this point starts checks and selection
  bool allVariablesPresent = (
    this->getConsumedValue(_electrons_rho_, rho)
    &&
    this->getConsumedValue(_electrons_conv_vtx_flag_, conv_vtx_flag)
    &&
    this->getConsumedValue(_electrons_charge_, charge)
    &&
    this->getConsumedValue(_electrons_expectedMissingInnerHits_, expectedMissingInnerHits)
    &&
    this->getConsumedValue(_electrons_energySC_, energySC)
    &&
    this->getConsumedValue(_electrons_etaSC_, etaSC)
    &&
    this->getConsumedValue(_electrons_etaSeedSC_, etaSeedSC)
    &&
    this->getConsumedValue(_electrons_sigmaIEtaIEta_full5x5_, sigmaIEtaIEta_full5x5)
    &&
    this->getConsumedValue(_electrons_dEtaIn_, dEtaIn)
    &&
    this->getConsumedValue(_electrons_dPhiIn_, dPhiIn)
    &&
    this->getConsumedValue(_electrons_hOverE_, hOverE)
    &&
    this->getConsumedValue(_electrons_ecalEnergy_, ecalEnergy)
    &&
    this->getConsumedValue(_electrons_eOverPIn_, eOverPIn)
    &&
    this->getConsumedValue(_electrons_dxyPV_, dxyPV)
    &&
    this->getConsumedValue(_electrons_dzPV_, dzPV)
    &&
    this->getConsumedValue(_electrons_miniIso_ch_, miniIso_ch)
    &&
    this->getConsumedValue(_electrons_miniIso_nh_, miniIso_nh)
    &&
    this->getConsumedValue(_electrons_miniIso_em_, miniIso_em)
    &&
    this->getConsumedValue(_electrons_momentum_, momentum)
    );

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "ElectronHandler::constructElectrons: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "ElectronHandler::constructElectrons: All variables are set up!" << endl;
  if (!charge){
    if (this->verbosity>=TVar::ERROR) MELAerr << "ElectronHandler::constructElectrons: Electrons could not be linked! Pointer to " << _electrons_charge_ << " is null!" << endl;
    return false;
  }

  if (charge->empty()) return true; // Construction is successful, it is just that no electrons exist.

  unsigned int nProducts = charge->size();
  productList.reserve(nProducts);
  for (unsigned int ip=0; ip<nProducts; ip++){
    if (this->verbosity>=TVar::DEBUG) MELAout << "ElectronHandler::constructElectrons: Attempting electron " << ip << "..." << endl;

    productList.push_back(new ElectronObject(-11*(charge->at(ip)>0 ? 1 : -1), momentum->at(ip)));
    ElectronObject*& obj = productList.back();

    obj->extras.rho = rho;

    obj->extras.conv_vtx_flag = conv_vtx_flag->at(ip);

    obj->extras.expectedMissingInnerHits = expectedMissingInnerHits->at(ip);

    obj->extras.energySC = energySC->at(ip);
    obj->extras.etaSC = etaSC->at(ip);
    obj->extras.etaSeedSC = etaSeedSC->at(ip);
    obj->extras.sigmaIEtaIEta_full5x5 = sigmaIEtaIEta_full5x5->at(ip);
    obj->extras.dEtaIn = dEtaIn->at(ip);
    obj->extras.dPhiIn = dPhiIn->at(ip);
    obj->extras.hOverE = hOverE->at(ip);
    obj->extras.ecalEnergy = ecalEnergy->at(ip);
    obj->extras.eOverPIn = eOverPIn->at(ip);
    obj->extras.dxyPV = dxyPV->at(ip);
    obj->extras.dzPV = dzPV->at(ip);
    obj->extras.miniIso_ch = miniIso_ch->at(ip);
    obj->extras.miniIso_nh = miniIso_nh->at(ip);
    obj->extras.miniIso_em = miniIso_em->at(ip);

    // Set the selection bits
    ElectronSelectionHelpers::setSelectionBits(*obj);

    if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;
  }

  return true;
}

void ElectronHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(tree);
  if (!fwktree) return;

  fwktree->bookEDMBranch<float>(_electrons_rho_, 0);

  fwktree->bookEDMBranch<std::vector<bool>*>(_electrons_conv_vtx_flag_, nullptr);

  fwktree->bookEDMBranch<std::vector<int>*>(_electrons_charge_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_electrons_expectedMissingInnerHits_, nullptr);

  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_energySC_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_etaSC_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_etaSeedSC_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_sigmaIEtaIEta_full5x5_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_dEtaIn_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_dPhiIn_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_hOverE_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_ecalEnergy_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_eOverPIn_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_dxyPV_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_dzPV_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_miniIso_ch_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_miniIso_nh_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_electrons_miniIso_em_, nullptr);

  fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_electrons_momentum_, nullptr);
}

