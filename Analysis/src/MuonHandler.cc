#include <cassert>
#include "ParticleObjectHelpers.h"
#include "MuonHandler.h"
#include "MuonSelectionHelpers.h"
#include "FrameworkVariables.hh"
#include "FrameworkTree.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


MuonHandler::MuonHandler() : IvyBase()
{
  this->addConsumed<float>(_muons_rho_);

  this->addConsumed<std::vector<unsigned int>*>(_muons_POGSelectorBit_);

  this->addConsumed<std::vector<int>*>(_muons_charge_);
  this->addConsumed<std::vector<int>*>(_muons_isPFMuon_);
  this->addConsumed<std::vector<int>*>(_muons_type_);
  this->addConsumed<std::vector<int>*>(_muons_validHits_);
  this->addConsumed<std::vector<int>*>(_muons_lostHits_);
  this->addConsumed<std::vector<int>*>(_muons_expectedMissingInnerHits_);
  this->addConsumed<std::vector<int>*>(_muons_expectedMissingOuterHits_);
  this->addConsumed<std::vector<int>*>(_muons_GlobalFit_Ndof_);

  this->addConsumed<std::vector<float>*>(_muons_GlobalFit_Chisq_);
  this->addConsumed<std::vector<float>*>(_muons_LocalPos_Chisq_);
  this->addConsumed<std::vector<float>*>(_muons_TrkKink_);
  this->addConsumed<std::vector<float>*>(_muons_SegComp_);
  this->addConsumed<std::vector<float>*>(_muons_dxyPV_);
  this->addConsumed<std::vector<float>*>(_muons_dzPV_);
  this->addConsumed<std::vector<float>*>(_muons_miniIso_ch_);
  this->addConsumed<std::vector<float>*>(_muons_miniIso_nh_);
  this->addConsumed<std::vector<float>*>(_muons_miniIso_em_);

  this->addConsumed<std::vector<CMSLorentzVector>*>(_muons_momentum_);
}


bool MuonHandler::constructMuons(){
  clear();
  if (!currentTree) return false;

  float rho = 0;

  std::vector<unsigned int>* POGSelectorBit = nullptr;

  std::vector<int>* charge = nullptr;
  std::vector<int>* isPFMuon = nullptr;
  std::vector<int>* type = nullptr;
  std::vector<int>* validHits = nullptr;
  std::vector<int>* lostHits = nullptr;
  std::vector<int>* expectedMissingInnerHits = nullptr;
  std::vector<int>* expectedMissingOuterHits = nullptr;
  std::vector<int>* GlobalFit_Ndof = nullptr;

  std::vector<float>* GlobalFit_Chisq = nullptr;
  std::vector<float>* LocalPos_Chisq = nullptr;
  std::vector<float>* TrkKink = nullptr;
  std::vector<float>* SegComp = nullptr;
  std::vector<float>* dxyPV = nullptr;
  std::vector<float>* dzPV = nullptr;
  std::vector<float>* miniIso_ch = nullptr;
  std::vector<float>* miniIso_nh = nullptr;
  std::vector<float>* miniIso_em = nullptr;

  std::vector<CMSLorentzVector>* momentum = nullptr;

  // Beyond this point starts checks and selection
  bool allVariablesPresent = (
    this->getConsumedValue(_muons_rho_, rho)
    &&
    this->getConsumedValue(_muons_POGSelectorBit_, POGSelectorBit)
    &&
    this->getConsumedValue(_muons_charge_, charge)
    &&
    this->getConsumedValue(_muons_isPFMuon_, isPFMuon)
    &&
    this->getConsumedValue(_muons_type_, type)
    &&
    this->getConsumedValue(_muons_validHits_, validHits)
    &&
    this->getConsumedValue(_muons_lostHits_, lostHits)
    &&
    this->getConsumedValue(_muons_expectedMissingInnerHits_, expectedMissingInnerHits)
    &&
    this->getConsumedValue(_muons_expectedMissingOuterHits_, expectedMissingOuterHits)
    &&
    this->getConsumedValue(_muons_GlobalFit_Ndof_, GlobalFit_Ndof)
    &&
    this->getConsumedValue(_muons_GlobalFit_Chisq_, GlobalFit_Chisq)
    &&
    this->getConsumedValue(_muons_LocalPos_Chisq_, LocalPos_Chisq)
    &&
    this->getConsumedValue(_muons_TrkKink_, TrkKink)
    &&
    this->getConsumedValue(_muons_SegComp_, SegComp)
    &&
    this->getConsumedValue(_muons_dxyPV_, dxyPV)
    &&
    this->getConsumedValue(_muons_dzPV_, dzPV)
    &&
    this->getConsumedValue(_muons_miniIso_ch_, miniIso_ch)
    &&
    this->getConsumedValue(_muons_miniIso_nh_, miniIso_nh)
    &&
    this->getConsumedValue(_muons_miniIso_em_, miniIso_em)
    &&
    this->getConsumedValue(_muons_momentum_, momentum)
    );

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "MuonHandler::constructMuons: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "MuonHandler::constructMuons: All variables are set up!" << endl;
  if (!charge){
    if (this->verbosity>=TVar::ERROR) MELAerr << "MuonHandler::constructMuons: Muons could not be linked! Pointer to " << _muons_charge_ << " is null!" << endl;
    return false;
  }

  if (charge->empty()) return true; // Construction is successful, it is just that no muons exist.

  unsigned int nProducts = charge->size();
  productList.reserve(nProducts);
  for (unsigned int ip=0; ip<nProducts; ip++){
    if (this->verbosity>=TVar::DEBUG) MELAout << "MuonHandler::constructMuons: Attempting muon " << ip << "..." << endl;

    productList.push_back(new MuonObject(-13*(charge->at(ip)>0 ? 1 : -1), momentum->at(ip)));
    MuonObject*& obj = productList.back();

    obj->extras.rho = rho;

    obj->extras.POGSelectorBit = (long long) POGSelectorBit->at(ip);

    obj->extras.isPFMuon = (bool) isPFMuon->at(ip);
    obj->extras.type = type->at(ip);
    obj->extras.validHits = validHits->at(ip);
    obj->extras.lostHits = lostHits->at(ip);
    obj->extras.expectedMissingInnerHits = expectedMissingInnerHits->at(ip);
    obj->extras.expectedMissingOuterHits = expectedMissingOuterHits->at(ip);
    obj->extras.GlobalFit_Ndof = GlobalFit_Ndof->at(ip);

    obj->extras.GlobalFit_Chisq = GlobalFit_Chisq->at(ip);
    obj->extras.LocalPos_Chisq = LocalPos_Chisq->at(ip);
    obj->extras.TrkKink = TrkKink->at(ip);
    obj->extras.SegComp = SegComp->at(ip);
    obj->extras.dxyPV = dxyPV->at(ip);
    obj->extras.dzPV = dzPV->at(ip);
    obj->extras.miniIso_ch = miniIso_ch->at(ip);
    obj->extras.miniIso_nh = miniIso_nh->at(ip);
    obj->extras.miniIso_em = miniIso_em->at(ip);

    // Set the selection bits
    MuonSelectionHelpers::setSelectionBits(*obj);

    if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  return true;
}

void MuonHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(tree);
  if (!fwktree) return;
  SampleHelpers::setupUsingOptions(fwktree->getOptions());

  fwktree->bookEDMBranch<float>(_muons_rho_, 0);

  fwktree->bookEDMBranch<std::vector<unsigned int>*>(_muons_POGSelectorBit_, nullptr);

  fwktree->bookEDMBranch<std::vector<int>*>(_muons_charge_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_muons_isPFMuon_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_muons_type_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_muons_validHits_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_muons_lostHits_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_muons_expectedMissingInnerHits_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_muons_expectedMissingOuterHits_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_muons_GlobalFit_Ndof_, nullptr);

  fwktree->bookEDMBranch<std::vector<float>*>(_muons_GlobalFit_Chisq_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_muons_LocalPos_Chisq_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_muons_TrkKink_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_muons_SegComp_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_muons_dxyPV_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_muons_dzPV_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_muons_miniIso_ch_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_muons_miniIso_nh_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_muons_miniIso_em_, nullptr);

  fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_muons_momentum_, nullptr);
}

