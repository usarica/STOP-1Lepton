#include <cassert>
#include "Samples.h"
#include "FrameworkVariables.hh"
#include "FrameworkTree.h"
#include "VertexPUHandler.h"
#include "VertexSelectionHelpers.h"
#include "SecondaryVertexSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


VertexPUHandler::VertexPUHandler() :
  IvyBase(),

  doVertices(true),
  doSecondaryVertices(true),
  doPUInfos(true)
{}


bool VertexPUHandler::constructVertices(){
  if (!doVertices) return true;

  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree) return false;

  std::vector<int>* isValid = nullptr;
  std::vector<int>* isFake = nullptr;
  std::vector<float>* ndof = nullptr;
  std::vector<CMSLorentzVector>* position = nullptr;

  bool allVariablesPresent = (
    this->getConsumedValue(_vtxs_isValid_, isValid)
    && this->getConsumedValue(_vtxs_isFake_, isFake)
    && this->getConsumedValue(_vtxs_ndof_, ndof)
    && this->getConsumedValue(_vtxs_position_, position)
    );
  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "VertexPUHandler::constructVertices: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  size_t nProducts = ndof->size();
  vertices.reserve(nProducts);
  for (size_t ivtx=0; ivtx<nProducts; ivtx++){
    if (this->verbosity>=TVar::DEBUG) MELAout << "VertexPUHandler::constructVertices: Attempting vertex " << ivtx << "..." << endl;

    vertices.push_back(new VertexObject(position->at(ivtx)));
    VertexObject*& obj = vertices.back();

    obj->isValid = isValid->at(ivtx);
    obj->isFake = isFake->at(ivtx);
    obj->ndof = ndof->at(ivtx);

    VertexSelectionHelpers::setSelectionBits(*obj); // Set the selection bits here directly

    if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;
  }

  return true;
}
bool VertexPUHandler::constructSecondaryVertices(){
  if (!doSecondaryVertices) return true;

  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree) return false;

  std::vector<int>* nTracks = nullptr;
  std::vector<float>* IP2D = nullptr;
  std::vector<float>* SIP2D = nullptr;
  std::vector<float>* IP3D = nullptr;
  std::vector<float>* SIP3D = nullptr;
  std::vector<float>* anglePV = nullptr;
  std::vector<CMSLorentzVector>* position = nullptr;
  std::vector<CMSLorentzVector>* momentum = nullptr;

  bool allVariablesPresent = (
    this->getConsumedValue(_svs_nTracks_, nTracks)
    && this->getConsumedValue(_svs_IP2D_, IP2D)
    && this->getConsumedValue(_svs_SIP2D_, SIP2D)
    && this->getConsumedValue(_svs_IP3D_, IP3D)
    && this->getConsumedValue(_svs_SIP3D_, SIP3D)
    && this->getConsumedValue(_svs_anglePV_, anglePV)
    && this->getConsumedValue(_svs_position_, position)
    && this->getConsumedValue(_svs_momentum_, momentum)
    );
  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "VertexPUHandler::constructSecondaryVertices: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  size_t nProducts = momentum->size();
  vertices.reserve(nProducts);
  for (size_t ivtx=0; ivtx<nProducts; ivtx++){
    if (this->verbosity>=TVar::DEBUG) MELAout << "VertexPUHandler::constructSecondaryVertices: Attempting vertex " << ivtx << "..." << endl;

    secondary_vertices.push_back(new SecondaryVertexObject());
    SecondaryVertexObject*& obj = secondary_vertices.back();

    CMSLorentzVector const& vec_pos = position->at(ivtx);

    obj->nTracks = nTracks->at(ivtx);
    obj->IP2D = IP2D->at(ivtx);
    obj->SIP2D = SIP2D->at(ivtx);
    obj->IP3D = IP3D->at(ivtx);
    obj->SIP3D = SIP3D->at(ivtx);
    obj->anglePV = anglePV->at(ivtx);
    obj->position = TVector3(vec_pos.x(), vec_pos.y(), vec_pos.z());
    obj->momentum = momentum->at(ivtx);

    SecondaryVertexSelectionHelpers::setSelectionBits(*obj); // Set the selection bits here directly

    if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;
  }

  return true;
}
bool VertexPUHandler::constructPUInfos(){
  if (!doPUInfos) return true;

  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree) return false;
  if (fwktree->isData()) return true; // Data trees do not contain PU info branches

  std::vector<int>* bunchCrossing = nullptr;
  std::vector<int>* nPUVertices = nullptr;
  std::vector<float>* nTrueVertices = nullptr;

  bool allVariablesPresent = (
    this->getConsumedValue(_puinfos_bunchCrossing_, bunchCrossing)
    && this->getConsumedValue(_puinfos_nPUVtxs_, nPUVertices)
    && this->getConsumedValue(_puinfos_nTrueVtxs_, nTrueVertices)
    );
  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "VertexPUHandler::constructPUInfos: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  size_t nProducts = nTrueVertices->size();
  puinfos.reserve(nProducts);
  for (size_t ipuinfo=0; ipuinfo<nProducts; ipuinfo++){
    if (this->verbosity>=TVar::DEBUG) MELAout << "VertexPUHandler::constructPUInfos: Attempting PU info " << ipuinfo << "..." << endl;

    puinfos.push_back(new PUInfoObject());
    PUInfoObject*& obj = puinfos.back();

    obj->bunchCrossing = bunchCrossing->at(ipuinfo);
    obj->nPUVertices = nPUVertices->at(ipuinfo);
    obj->nTrueVertices = nTrueVertices->at(ipuinfo);

    if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;
  }

  return true;
}
bool VertexPUHandler::constructVertexPUInfos(){
  clear();
  if (!currentTree) return false;

  return (constructVertices() && constructSecondaryVertices() && constructPUInfos());
}

void VertexPUHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(tree);
  if (!fwktree) return;
  SampleHelpers::setupUsingOptions(fwktree->getOptions());

  if (doVertices){
    this->addConsumed<std::vector<int>*>(_vtxs_isValid_);
    this->addConsumed<std::vector<int>*>(_vtxs_isFake_);
    this->addConsumed<std::vector<float>*>(_vtxs_ndof_);
    this->addConsumed<std::vector<CMSLorentzVector>*>(_vtxs_position_);

    fwktree->bookEDMBranch<std::vector<int>*>(_vtxs_isValid_, nullptr);
    fwktree->bookEDMBranch<std::vector<int>*>(_vtxs_isFake_, nullptr);
    fwktree->bookEDMBranch<std::vector<float>*>(_vtxs_ndof_, nullptr);
    fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_vtxs_position_, nullptr);
  }
  if (doSecondaryVertices){
    this->addConsumed<std::vector<int>*>(_svs_nTracks_);
    this->addConsumed<std::vector<float>*>(_svs_IP2D_);
    this->addConsumed<std::vector<float>*>(_svs_SIP2D_);
    this->addConsumed<std::vector<float>*>(_svs_IP3D_);
    this->addConsumed<std::vector<float>*>(_svs_SIP3D_);
    this->addConsumed<std::vector<float>*>(_svs_anglePV_);
    this->addConsumed<std::vector<CMSLorentzVector>*>(_svs_position_);
    this->addConsumed<std::vector<CMSLorentzVector>*>(_svs_momentum_);

    fwktree->bookEDMBranch<std::vector<int>*>(_svs_nTracks_, nullptr);
    fwktree->bookEDMBranch<std::vector<float>*>(_svs_IP2D_, nullptr);
    fwktree->bookEDMBranch<std::vector<float>*>(_svs_SIP2D_, nullptr);
    fwktree->bookEDMBranch<std::vector<float>*>(_svs_IP3D_, nullptr);
    fwktree->bookEDMBranch<std::vector<float>*>(_svs_SIP3D_, nullptr);
    fwktree->bookEDMBranch<std::vector<float>*>(_svs_anglePV_, nullptr);
    fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_svs_position_, nullptr);
    fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_svs_momentum_, nullptr);
  }
  if (doPUInfos && fwktree->isMC()){
    this->addConsumed<std::vector<int>*>(_puinfos_bunchCrossing_);
    this->addConsumed<std::vector<int>*>(_puinfos_nPUVtxs_);
    this->addConsumed<std::vector<float>*>(_puinfos_nTrueVtxs_);

    this->defineConsumedSloppy(_puinfos_bunchCrossing_); // Define this sloppy for data will not have it
    this->defineConsumedSloppy(_puinfos_nPUVtxs_); // Define this sloppy for data will not have it
    this->defineConsumedSloppy(_puinfos_nTrueVtxs_); // Define this sloppy for data will not have it

    fwktree->bookEDMBranch<std::vector<int>*>(_puinfos_bunchCrossing_, nullptr);
    fwktree->bookEDMBranch<std::vector<int>*>(_puinfos_nPUVtxs_, nullptr);
    fwktree->bookEDMBranch<std::vector<float>*>(_puinfos_nTrueVtxs_, nullptr);
  }
}
