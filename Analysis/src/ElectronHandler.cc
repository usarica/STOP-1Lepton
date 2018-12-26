#include "ElectronHandler.h"
#include "FrameworkVariables.hh"


using namespace std;


ElectronHandler::ElectronHandler() : IvyBase()
{
  this->addConsumed<float>(_electrons_rho_);

  this->addConsumed<std::vector<bool>>(_electrons_conv_vtx_flag_);

  this->addConsumed<std::vector<int>>(_electrons_charge_);
  this->addConsumed<std::vector<int>>(_electrons_expectedMissingInnerHits_);

  this->addConsumed<std::vector<float>>(_electrons_energySC_);
  this->addConsumed<std::vector<float>>(_electrons_etaSC_);
  this->addConsumed<std::vector<float>>(_electrons_etaSeedSC_);
  this->addConsumed<std::vector<float>>(_electrons_sigmaIEtaIEta_full5x5_);
  this->addConsumed<std::vector<float>>(_electrons_dEtaIn_);
  this->addConsumed<std::vector<float>>(_electrons_dPhiIn_);
  this->addConsumed<std::vector<float>>(_electrons_hOverE_);
  this->addConsumed<std::vector<float>>(_electrons_ecalEnergy_);
  this->addConsumed<std::vector<float>>(_electrons_eOverPIn_);
  this->addConsumed<std::vector<float>>(_electrons_dxyPV_);
  this->addConsumed<std::vector<float>>(_electrons_dzPV_);
  this->addConsumed<std::vector<float>>(_electrons_miniIso_ch_);
  this->addConsumed<std::vector<float>>(_electrons_miniIso_nh_);
  this->addConsumed<std::vector<float>>(_electrons_miniIso_em_);

  this->addConsumed<std::vector<CMSLorentzVector>>(_electrons_momentum_);
}


bool ElectronHandler::constructElectrons(){
  clear();
  if (!currentTree) return false;

  vector<int> const* charge = valVints[_electrons_charge_];
  vector<CMSLorentzVector> const* momentum = valVCMSLorentzVectors[_electrons_momentum_];

  float const* rho = valfloats[_electrons_rho_];

  vector<bool> const* conv_vtx_flag = valVbools[_electrons_conv_vtx_flag_];

  vector<int> const* expectedMissingInnerHits = valVints[_electrons_expectedMissingInnerHits_];

  vector<float> const* energySC = valVfloats[_electrons_energySC_];
  vector<float> const* etaSC = valVfloats[_electrons_etaSC_];
  vector<float> const* etaSeedSC = valVfloats[_electrons_etaSeedSC_];
  vector<float> const* sigmaIEtaIEta_full5x5 = valVfloats[_electrons_sigmaIEtaIEta_full5x5_];
  vector<float> const* dEtaIn = valVfloats[_electrons_dEtaIn_];
  vector<float> const* dPhiIn = valVfloats[_electrons_dPhiIn_];
  vector<float> const* hOverE = valVfloats[_electrons_hOverE_];
  vector<float> const* ecalEnergy = valVfloats[_electrons_ecalEnergy_];
  vector<float> const* eOverPIn = valVfloats[_electrons_eOverPIn_];
  vector<float> const* dxyPV = valVfloats[_electrons_dxyPV_];
  vector<float> const* dzPV = valVfloats[_electrons_dzPV_];
  vector<float> const* miniIso_ch = valVfloats[_electrons_miniIso_ch_];
  vector<float> const* miniIso_nh = valVfloats[_electrons_miniIso_nh_];
  vector<float> const* miniIso_em = valVfloats[_electrons_miniIso_em_];

  if (charge->empty()) return false;

  unsigned int nProducts = charge->size();
  productList.reserve(nProducts);
  for (unsigned int ip=0; ip<nProducts; ip++){
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
  }

  return true;
}

void ElectronHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  tree->bookBranch<std::vector<bool>*>(_electrons_conv_vtx_flag_, nullptr);

  tree->bookBranch<std::vector<int>*>(_electrons_charge_, nullptr);
  tree->bookBranch<std::vector<int>*>(_electrons_expectedMissingInnerHits_, nullptr);

  tree->bookBranch<std::vector<float>*>(_electrons_energySC_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_etaSC_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_etaSeedSC_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_rho_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_sigmaIEtaIEta_full5x5_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_dEtaIn_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_dPhiIn_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_hOverE_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_ecalEnergy_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_eOverPIn_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_dxyPV_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_dzPV_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_miniIso_ch_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_miniIso_nh_, nullptr);
  tree->bookBranch<std::vector<float>*>(_electrons_miniIso_em_, nullptr);

  tree->bookBranch<std::vector<CMSLorentzVector>*>(_electrons_momentum_, nullptr);
}

