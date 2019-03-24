#include <cassert>
#include "ParticleObjectHelpers.h"
#include "PFCandHandler.h"
#include "FrameworkVariables.hh"
#include "FrameworkTree.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


PFCandHandler::PFCandHandler() : IvyBase()
{
  this->addConsumed<std::vector<bool>*>(_pfcands_trackHighPurity_);

  this->addConsumed<std::vector<int>*>(_pfcands_charge_);
  this->addConsumed<std::vector<int>*>(_pfcands_id_);

  this->addConsumed<std::vector<float>*>(_pfcands_dxy_);
  this->addConsumed<std::vector<float>*>(_pfcands_dz_);
  this->addConsumed<std::vector<float>*>(_pfcands_dxyError_);
  this->addConsumed<std::vector<float>*>(_pfcands_dzError_);

  this->addConsumed<std::vector<CMSLorentzVector>*>(_pfcands_momentum_);
}


bool PFCandHandler::constructPFCands(){
  clear();
  if (!currentTree) return false;

  std::vector<bool>* trackHighPurity = nullptr;

  std::vector<int>* charge = nullptr;
  std::vector<int>* id = nullptr;

  std::vector<float>* dxy = nullptr;
  std::vector<float>* dz = nullptr;
  std::vector<float>* dxyError = nullptr;
  std::vector<float>* dzError = nullptr;

  std::vector<CMSLorentzVector>* momentum = nullptr;

  // Beyond this point starts checks and selection
  bool allVariablesPresent = (
    this->getConsumedValue(_pfcands_trackHighPurity_, trackHighPurity)
    &&
    this->getConsumedValue(_pfcands_charge_, charge)
    &&
    this->getConsumedValue(_pfcands_id_, id)
    &&
    this->getConsumedValue(_pfcands_dxy_, dxy)
    &&
    this->getConsumedValue(_pfcands_dz_, dz)
    &&
    this->getConsumedValue(_pfcands_dxyError_, dxyError)
    &&
    this->getConsumedValue(_pfcands_dzError_, dzError)
    &&
    this->getConsumedValue(_pfcands_momentum_, momentum)
    );

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "PFCandHandler::constructPFCands: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "PFCandHandler::constructPFCands: All variables are set up!" << endl;
  if (!charge){
    if (this->verbosity>=TVar::ERROR) MELAerr << "PFCandHandler::constructPFCands: PFCands could not be linked! Pointer to " << _pfcands_charge_ << " is null!" << endl;
    return false;
  }

  if (id->empty()) return true; // Construction is successful, it is just that no pfcands exist.

  unsigned int nProducts = id->size();
  productList.reserve(nProducts);
  for (unsigned int ip=0; ip<nProducts; ip++){
    if (this->verbosity>=TVar::DEBUG) MELAout << "PFCandHandler::constructPFCands: Attempting PFCand " << ip << "..." << endl;

    productList.push_back(new PFCandObject(id->at(ip), momentum->at(ip)));
    PFCandObject*& obj = productList.back();

    obj->extras.trackHighPurity = trackHighPurity->at(ip);

    obj->extras.charge = charge->at(ip);

    obj->extras.dxy = dxy->at(ip);
    obj->extras.dz = dz->at(ip);
    obj->extras.dxyError = dxyError->at(ip);
    obj->extras.dzError = dzError->at(ip);

    if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  return true;
}

void PFCandHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(tree);
  if (!fwktree) return;
  SampleHelpers::setupUsingOptions(fwktree->getOptions());


  fwktree->bookEDMBranch<std::vector<bool>*>(_pfcands_trackHighPurity_, nullptr);

  fwktree->bookEDMBranch<std::vector<int>*>(_pfcands_charge_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_pfcands_id_, nullptr);

  fwktree->bookEDMBranch<std::vector<float>*>(_pfcands_dxy_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_pfcands_dz_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_pfcands_dxyError_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_pfcands_dzError_, nullptr);

  fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_pfcands_momentum_, nullptr);
}

