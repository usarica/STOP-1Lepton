#include "ElectronHandler.h"
#include "FrameworkVariables.hh"
#include "FrameworkTree.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


ElectronHandler::ElectronHandler() : IvyBase()
{
  this->addConsumed<float>(_electrons_rho_);

  this->addConsumed<std::vector<bool>* const>(_electrons_conv_vtx_flag_);

  this->addConsumed<std::vector<int>* const>(_electrons_charge_);
  this->addConsumed<std::vector<int>* const>(_electrons_expectedMissingInnerHits_);

  this->addConsumed<std::vector<float>* const>(_electrons_energySC_);
  this->addConsumed<std::vector<float>* const>(_electrons_etaSC_);
  this->addConsumed<std::vector<float>* const>(_electrons_etaSeedSC_);
  this->addConsumed<std::vector<float>* const>(_electrons_sigmaIEtaIEta_full5x5_);
  this->addConsumed<std::vector<float>* const>(_electrons_dEtaIn_);
  this->addConsumed<std::vector<float>* const>(_electrons_dPhiIn_);
  this->addConsumed<std::vector<float>* const>(_electrons_hOverE_);
  this->addConsumed<std::vector<float>* const>(_electrons_ecalEnergy_);
  this->addConsumed<std::vector<float>* const>(_electrons_eOverPIn_);
  this->addConsumed<std::vector<float>* const>(_electrons_dxyPV_);
  this->addConsumed<std::vector<float>* const>(_electrons_dzPV_);
  this->addConsumed<std::vector<float>* const>(_electrons_miniIso_ch_);
  this->addConsumed<std::vector<float>* const>(_electrons_miniIso_nh_);
  this->addConsumed<std::vector<float>* const>(_electrons_miniIso_em_);

  this->addConsumed<std::vector<CMSLorentzVector>* const>(_electrons_momentum_);
}


bool ElectronHandler::constructElectrons(){
  clear();
  if (!currentTree) return false;

  vector<int> const* charge = nullptr;
  vector<int>* const* chargePtr = valVints[_electrons_charge_];
  vector<CMSLorentzVector> const* momentum = nullptr;
  vector<CMSLorentzVector>* const* momentumPtr = valVCMSLorentzVectors[_electrons_momentum_];

  float const* rho = valfloats[_electrons_rho_];

  vector<bool> const* conv_vtx_flag = nullptr;
  vector<bool>* const* conv_vtx_flagPtr = valVbools[_electrons_conv_vtx_flag_];

  vector<int> const* expectedMissingInnerHits = nullptr;
  vector<int>* const* expectedMissingInnerHitsPtr = valVints[_electrons_expectedMissingInnerHits_];

  vector<float> const* energySC = nullptr;
  vector<float>* const* energySCPtr = valVfloats[_electrons_energySC_];
  vector<float> const* etaSC = nullptr;
  vector<float>* const* etaSCPtr = valVfloats[_electrons_etaSC_];
  vector<float> const* etaSeedSC = nullptr;
  vector<float>* const* etaSeedSCPtr = valVfloats[_electrons_etaSeedSC_];
  vector<float> const* sigmaIEtaIEta_full5x5 = nullptr;
  vector<float>* const* sigmaIEtaIEta_full5x5Ptr = valVfloats[_electrons_sigmaIEtaIEta_full5x5_];
  vector<float> const* dEtaIn = nullptr;
  vector<float>* const* dEtaInPtr = valVfloats[_electrons_dEtaIn_];
  vector<float> const* dPhiIn = nullptr;
  vector<float>* const* dPhiInPtr = valVfloats[_electrons_dPhiIn_];
  vector<float> const* hOverE = nullptr;
  vector<float>* const* hOverEPtr = valVfloats[_electrons_hOverE_];
  vector<float> const* ecalEnergy = nullptr;
  vector<float>* const* ecalEnergyPtr = valVfloats[_electrons_ecalEnergy_];
  vector<float> const* eOverPIn = nullptr;
  vector<float>* const* eOverPInPtr = valVfloats[_electrons_eOverPIn_];
  vector<float> const* dxyPV = nullptr;
  vector<float>* const* dxyPVPtr = valVfloats[_electrons_dxyPV_];
  vector<float> const* dzPV = nullptr;
  vector<float>* const* dzPVPtr = valVfloats[_electrons_dzPV_];
  vector<float> const* miniIso_ch = nullptr;
  vector<float>* const* miniIso_chPtr = valVfloats[_electrons_miniIso_ch_];
  vector<float> const* miniIso_nh = nullptr;
  vector<float>* const* miniIso_nhPtr = valVfloats[_electrons_miniIso_nh_];
  vector<float> const* miniIso_em = nullptr;
  vector<float>* const* miniIso_emPtr = valVfloats[_electrons_miniIso_em_];

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

    obj->extras.rho = *rho;

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

