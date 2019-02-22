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

#define VECTOR_ITERATOR_HANDLER_DIRECTIVES \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<unsigned int>, POGSelectorBit) \
\
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, charge) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, isPFMuon) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, type) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, validHits) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, lostHits) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, expectedMissingInnerHits) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, expectedMissingOuterHits) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<int>, GlobalFit_Ndof) \
\
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, GlobalFit_Chisq) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, LocalPos_Chisq) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, TrkKink) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, SegComp) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, dxyPV) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, dzPV) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, miniIso_ch) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, miniIso_nh) \
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<float>, miniIso_em) \
\
VECTOR_ITERATOR_HANDLER_DIRECTIVE(std::vector<CMSLorentzVector>, momentum)

#define VECTOR_ITERATOR_HANDLER_DIRECTIVE(TYPE, NAME) TYPE::const_iterator itBegin_##NAME, itEnd_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef VECTOR_ITERATOR_HANDLER_DIRECTIVE

#define VECTOR_ITERATOR_HANDLER_DIRECTIVE(TYPE, NAME) && this->getConsumedCIterators<TYPE>(_muons_##NAME##_, &itBegin_##NAME, &itEnd_##NAME)
  // Beyond this point starts checks and selection
  bool allVariablesPresent = (
    this->getConsumedValue(_muons_rho_, rho)
    VECTOR_ITERATOR_HANDLER_DIRECTIVES
    );
#undef VECTOR_ITERATOR_HANDLER_DIRECTIVE

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "MuonHandler::constructMuons: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "MuonHandler::constructMuons: All variables are set up!" << endl;

  if (itBegin_charge == itEnd_charge) return true; // Construction is successful, it is just that no muons exist.

  size_t nProducts = (itEnd_charge - itBegin_charge);
  productList.reserve(nProducts);
#define VECTOR_ITERATOR_HANDLER_DIRECTIVE(TYPE, NAME) auto it_##NAME = itBegin_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef VECTOR_ITERATOR_HANDLER_DIRECTIVE
#define VECTOR_ITERATOR_HANDLER_DIRECTIVE(TYPE, NAME) it_##NAME++;
  {
    size_t ip=0;
    while (it_charge != itEnd_charge){
      if (this->verbosity>=TVar::DEBUG) MELAout << "MuonHandler::constructMuons: Attempting muon " << ip << "..." << endl;

      productList.push_back(new MuonObject(-13*(*it_charge>0 ? 1 : -1), *it_momentum));
      MuonObject*& obj = productList.back();

      obj->extras.rho = rho;

      obj->extras.POGSelectorBit = (long long) *it_POGSelectorBit;

      obj->extras.isPFMuon = (bool) *it_isPFMuon;
      obj->extras.type = *it_type;
      obj->extras.validHits = *it_validHits;
      obj->extras.lostHits = *it_lostHits;
      obj->extras.expectedMissingInnerHits = *it_expectedMissingInnerHits;
      obj->extras.expectedMissingOuterHits = *it_expectedMissingOuterHits;
      obj->extras.GlobalFit_Ndof = *it_GlobalFit_Ndof;

      obj->extras.GlobalFit_Chisq = *it_GlobalFit_Chisq;
      obj->extras.LocalPos_Chisq = *it_LocalPos_Chisq;
      obj->extras.TrkKink = *it_TrkKink;
      obj->extras.SegComp = *it_SegComp;
      obj->extras.dxyPV = *it_dxyPV;
      obj->extras.dzPV = *it_dzPV;
      obj->extras.miniIso_ch = *it_miniIso_ch;
      obj->extras.miniIso_nh = *it_miniIso_nh;
      obj->extras.miniIso_em = *it_miniIso_em;

      // Set the selection bits
      MuonSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES
    }
  }
#undef VECTOR_ITERATOR_HANDLER_DIRECTIVE
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
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

