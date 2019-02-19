#include <cassert>
#include "ParticleObjectHelpers.h"
#include "TauHandler.h"
#include "MuonSelectionHelpers.h"
#include "ElectronSelectionHelpers.h"
#include "TauSelectionHelpers.h"
#include "FrameworkVariables.hh"
#include "FrameworkTree.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


TauHandler::TauHandler() :
  IvyBase(),

  registeredMuons(nullptr),
  registeredElectrons(nullptr)
{
  this->addConsumed<std::vector<int>*>(_taus_charge_);

  this->addConsumed<std::vector<TString>*>(_taus_pfidnames_);
  this->addConsumed<std::vector<std::vector<float>>*>(_taus_pfids_);

  this->addConsumed<std::vector<CMSLorentzVector>*>(_taus_momentum_);
}


bool TauHandler::constructTaus(){
  clear();
  if (!currentTree) return false;

  std::vector<int>* charge = nullptr;

  std::vector<TString>* pfidnames = nullptr;
  std::vector<std::vector<float>>* pfids = nullptr;

  std::vector<CMSLorentzVector>* momentum = nullptr;

  // Beyond this point starts checks and selection
  bool allVariablesPresent = (
    this->getConsumedValue(_taus_charge_, charge)
    &&
    this->getConsumedValue(_taus_pfidnames_, pfidnames)
    &&
    this->getConsumedValue(_taus_pfids_, pfids)
    &&
    this->getConsumedValue(_taus_momentum_, momentum)
    );

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "TauHandler::constructTaus: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "TauHandler::constructTaus: All variables are set up!" << endl;
  if (!charge){
    if (this->verbosity>=TVar::ERROR) MELAerr << "TauHandler::constructTaus: Taus could not be linked! Pointer to " << _taus_charge_ << " is null!" << endl;
    return false;
  }

  if (charge->empty() || pfidnames->empty()) return true; // Construction is successful, it is just that no taus exist.
  auto it_begin = pfidnames->cbegin();
  auto it_end = pfidnames->cend();
  auto it_decaymode = std::find(it_begin, it_end, "decayModeFinding");
  auto it_pfiso = std::find(it_begin, it_end, TauSelectionHelpers::strPFIsoVar);
  if (this->verbosity>=TVar::ERROR){
    if (it_decaymode==pfidnames->end()){
      MELAerr << "TauHandler::constructTaus: Decay mode variable could not be found!" << endl;
      assert(0);
    }
    if (it_pfiso==pfidnames->end()){
      MELAerr << "TauHandler::constructTaus: Isolation variable could not be found!" << endl;
      assert(0);
    }
  }
  const size_t iDecayMode = (it_decaymode-it_begin);
  const size_t iPFIso = (it_pfiso-it_begin);

  unsigned int nProducts = charge->size();
  productList.reserve(nProducts);
  for (unsigned int ip=0; ip<nProducts; ip++){
    if (this->verbosity>=TVar::DEBUG) MELAout << "TauHandler::constructTaus: Attempting tau " << ip << "..." << endl;

    CMSLorentzVector const& mom = momentum->at(ip);

    // Check whether there is a deltaR matching and clean the iso. tracks
    // Note that the actual STOP-1L event veto also requires OS requirement with the leading-pT lepton
    bool doSkip = false;
    if (registeredMuons && !doSkip){
      for (auto const* part:(*registeredMuons)){
        if (!(part->testSelection(MuonSelectionHelpers::kVetoIDReco) && part->testSelection(MuonSelectionHelpers::kSkimPtEta))) continue;
        if (part->deltaR(mom)<TauSelectionHelpers::deltaR_veto_comparison){ doSkip = true; break; }
      }
    }
    if (registeredElectrons && !doSkip){
      for (auto const* part:(*registeredElectrons)){
        if (!(part->testSelection(ElectronSelectionHelpers::kVetoIDReco) && part->testSelection(ElectronSelectionHelpers::kSkimPtEta))) continue;
        if (part->deltaR(mom)<TauSelectionHelpers::deltaR_veto_comparison){ doSkip = true; break; }
      }
    }
    if (doSkip) continue;

    productList.push_back(new TauObject(-15*(charge->at(ip)<0 ? 1 : -1), mom));
    TauObject*& obj = productList.back();

    obj->extras.charge = charge->at(ip);

    obj->extras.pfDecayModeFinding = (pfids->at(ip).at(iDecayMode)>=1.f);
    obj->extras.pfIso = (pfids->at(ip).at(iPFIso)>=1.f);

    // Set the selection bits
    TauSelectionHelpers::setSelectionBits(*obj);

    if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  // Reset registered particles
  registeredMuons = nullptr;
  registeredElectrons = nullptr;

  return true;
}

void TauHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(tree);
  if (!fwktree) return;
  SampleHelpers::setupUsingOptions(fwktree->getOptions());

  fwktree->bookEDMBranch<std::vector<int>*>(_taus_charge_, nullptr);

  fwktree->bookEDMBranch<std::vector<TString>*>(_taus_pfidnames_, nullptr);
  fwktree->bookEDMBranch<std::vector<std::vector<float>>*>(_taus_pfids_, nullptr);

  fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_taus_momentum_, nullptr);
}
