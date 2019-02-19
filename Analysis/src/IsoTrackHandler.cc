#include <cassert>
#include "ParticleObjectHelpers.h"
#include "IsoTrackHandler.h"
#include "MuonSelectionHelpers.h"
#include "ElectronSelectionHelpers.h"
#include "IsoTrackSelectionHelpers.h"
#include "FrameworkVariables.hh"
#include "FrameworkTree.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


IsoTrackHandler::IsoTrackHandler() :
  IvyBase(),

  registeredMuons(nullptr),
  registeredElectrons(nullptr)
{
  this->addConsumed<std::vector<bool>*>(_isotracks_isPFCand_);
  this->addConsumed<std::vector<bool>*>(_isotracks_hasLepOverlap_);

  this->addConsumed<std::vector<int>*>(_isotracks_charge_);
  this->addConsumed<std::vector<int>*>(_isotracks_id_);

  this->addConsumed<std::vector<float>*>(_isotracks_pfIso_ch_);
  this->addConsumed<std::vector<float>*>(_isotracks_dz_);

  this->addConsumed<std::vector<CMSLorentzVector>*>(_isotracks_momentum_);
}


bool IsoTrackHandler::constructIsoTracks(){
  clear();
  if (!currentTree) return false;

  std::vector<bool>* isPFCand = nullptr;
  std::vector<bool>* hasLepOverlap = nullptr;

  std::vector<int>* charge = nullptr;
  std::vector<int>* id = nullptr;

  std::vector<float>* pfIso_ch = nullptr;
  std::vector<float>* dz = nullptr;

  std::vector<CMSLorentzVector>* momentum = nullptr;

  // Beyond this point starts checks and selection
  bool allVariablesPresent = (
    this->getConsumedValue(_isotracks_isPFCand_, isPFCand)
    &&
    this->getConsumedValue(_isotracks_hasLepOverlap_, hasLepOverlap)
    &&
    this->getConsumedValue(_isotracks_charge_, charge)
    &&
    this->getConsumedValue(_isotracks_id_, id)
    &&
    this->getConsumedValue(_isotracks_pfIso_ch_, pfIso_ch)
    &&
    this->getConsumedValue(_isotracks_dz_, dz)
    &&
    this->getConsumedValue(_isotracks_momentum_, momentum)
    );

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "IsoTrackHandler::constructIsoTracks: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "IsoTrackHandler::constructIsoTracks: All variables are set up!" << endl;
  if (!charge){
    if (this->verbosity>=TVar::ERROR) MELAerr << "IsoTrackHandler::constructIsoTracks: IsoTracks could not be linked! Pointer to " << _isotracks_charge_ << " is null!" << endl;
    return false;
  }

  if (charge->empty()) return true; // Construction is successful, it is just that no isotracks exist.

  unsigned int nProducts = charge->size();
  productList.reserve(nProducts);
  for (unsigned int ip=0; ip<nProducts; ip++){
    if (this->verbosity>=TVar::DEBUG) MELAout << "IsoTrackHandler::constructIsoTracks: Attempting isotrack " << ip << "..." << endl;

    CMSLorentzVector const& mom = momentum->at(ip);

    // Check whether there is a deltaR matching and clean the iso. tracks
    // Note that the actual STOP-1L event veto also requires OS requirement with the leading-pT lepton
    bool doSkip = false;
    if (registeredMuons && !doSkip){
      for (auto const* part:(*registeredMuons)){
        if (!(part->testSelection(MuonSelectionHelpers::kVetoIDReco) && part->testSelection(MuonSelectionHelpers::kSkimPtEta))) continue;
        if (part->deltaR(mom)<IsoTrackSelectionHelpers::deltaR_veto_comparison){ doSkip = true; break; }
      }
    }
    if (registeredElectrons && !doSkip){
      for (auto const* part:(*registeredElectrons)){
        if (!(part->testSelection(ElectronSelectionHelpers::kVetoIDReco) && part->testSelection(ElectronSelectionHelpers::kSkimPtEta))) continue;
        if (part->deltaR(mom)<IsoTrackSelectionHelpers::deltaR_veto_comparison){ doSkip = true; break; }
      }
    }
    if (doSkip) continue;

    productList.push_back(new IsoTrackObject(id->at(ip), mom));
    IsoTrackObject*& obj = productList.back();

    obj->extras.isPFCand = isPFCand->at(ip);
    obj->extras.hasLepOverlap = hasLepOverlap->at(ip);

    obj->extras.charge = charge->at(ip);

    obj->extras.pfIso_ch = pfIso_ch->at(ip);
    obj->extras.dz = dz->at(ip);

    // Set the selection bits
    IsoTrackSelectionHelpers::setSelectionBits(*obj);

    if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  // Reset registered particles
  registeredMuons = nullptr;
  registeredElectrons = nullptr;

  return true;
}

void IsoTrackHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(tree);
  if (!fwktree) return;
  SampleHelpers::setupUsingOptions(fwktree->getOptions());

  fwktree->bookEDMBranch<std::vector<bool>*>(_isotracks_isPFCand_, nullptr);
  fwktree->bookEDMBranch<std::vector<bool>*>(_isotracks_hasLepOverlap_, nullptr);

  fwktree->bookEDMBranch<std::vector<int>*>(_isotracks_charge_, nullptr);
  fwktree->bookEDMBranch<std::vector<int>*>(_isotracks_id_, nullptr);

  fwktree->bookEDMBranch<std::vector<float>*>(_isotracks_pfIso_ch_, nullptr);
  fwktree->bookEDMBranch<std::vector<float>*>(_isotracks_dz_, nullptr);

  fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_isotracks_momentum_, nullptr);
}
