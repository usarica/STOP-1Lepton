#include <cassert>
#include "Samples.h"
#include "SampleHelpers.h"
#include "ParticleObjectHelpers.h"
#include "FrameworkVariables.hh"
#include "GenInfoHandler.h"
#include "FrameworkTree.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


GenInfoHandler::GenInfoHandler() :
  IvyBase(),
  doEventInfo(true),
  doParticleInfo(true),
  geninfo(nullptr)
{
  // Add consumed for gen. info
  if (doEventInfo){
    this->addConsumed<unsigned int>(_geninfo_processID_);
    this->addConsumed<float>(_geninfo_qscale_);
    this->addConsumed<float>(_geninfo_alphaS_);
    this->addConsumed<float>(_gen_met_);
    this->addConsumed<float>(_gen_metPhi_);
    // Define as sloppy
    this->defineConsumedSloppy(_geninfo_processID_);
    this->defineConsumedSloppy(_geninfo_qscale_);
    this->defineConsumedSloppy(_geninfo_alphaS_);
    this->defineConsumedSloppy(_gen_met_);
    this->defineConsumedSloppy(_gen_metPhi_);
  }

  // Add consumed for gen. particles
  if (doParticleInfo){
    this->addConsumed<std::vector<bool>*>(_genparticles_isPromptFinalState_);
    this->addConsumed<std::vector<bool>*>(_genparticles_isPromptDecayed_);
    this->addConsumed<std::vector<bool>*>(_genparticles_isDirectPromptTauDecayProductFinalState_);
    this->addConsumed<std::vector<bool>*>(_genparticles_isHardProcess_);
    this->addConsumed<std::vector<bool>*>(_genparticles_fromHardProcessFinalState_);
    this->addConsumed<std::vector<bool>*>(_genparticles_fromHardProcessDecayed_);
    this->addConsumed<std::vector<bool>*>(_genparticles_isDirectHardProcessTauDecayProductFinalState_);
    this->addConsumed<std::vector<bool>*>(_genparticles_fromHardProcessBeforeFSR_);
    this->addConsumed<std::vector<bool>*>(_genparticles_isLastCopy_);
    this->addConsumed<std::vector<bool>*>(_genparticles_isLastCopyBeforeFSR_);

    this->addConsumed<std::vector<int>*>(_genparticles_id_);
    this->addConsumed<std::vector<int>*>(_genparticles_status_);

    this->addConsumed<std::vector<CMSLorentzVector>*>(_genparticles_p4_);
    // Define as sloppy
    this->defineConsumedSloppy(_genparticles_isPromptFinalState_);
    this->defineConsumedSloppy(_genparticles_isPromptDecayed_);
    this->defineConsumedSloppy(_genparticles_isDirectPromptTauDecayProductFinalState_);
    this->defineConsumedSloppy(_genparticles_isHardProcess_);
    this->defineConsumedSloppy(_genparticles_fromHardProcessFinalState_);
    this->defineConsumedSloppy(_genparticles_fromHardProcessDecayed_);
    this->defineConsumedSloppy(_genparticles_isDirectHardProcessTauDecayProductFinalState_);
    this->defineConsumedSloppy(_genparticles_fromHardProcessBeforeFSR_);
    this->defineConsumedSloppy(_genparticles_isLastCopy_);
    this->defineConsumedSloppy(_genparticles_isLastCopyBeforeFSR_);

    this->defineConsumedSloppy(_genparticles_id_);
    this->defineConsumedSloppy(_genparticles_status_);

    this->defineConsumedSloppy(_genparticles_p4_);
  }
}

void GenInfoHandler::clear(){
  // Clear everything in reverse order of creation
  for (auto& obj:genparticles) delete obj;
  genparticles.clear();

  delete geninfo; geninfo=nullptr;
}


bool GenInfoHandler::constructGenInfo(){
  clear();
  if (!currentTree) return false;

  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree) return false;
  if (!fwktree->isMC()) return true;

  unsigned int geninfo_processID = 0;
  float geninfo_qscale = 0;
  float geninfo_alphaS = 0;
  float gen_met = 0;
  float gen_metPhi = 0;

  std::vector<bool>* genparticles_isPromptFinalState = nullptr;
  std::vector<bool>* genparticles_isPromptDecayed = nullptr;
  std::vector<bool>* genparticles_isDirectPromptTauDecayProductFinalState = nullptr;
  std::vector<bool>* genparticles_isHardProcess = nullptr;
  std::vector<bool>* genparticles_fromHardProcessFinalState = nullptr;
  std::vector<bool>* genparticles_fromHardProcessDecayed = nullptr;
  std::vector<bool>* genparticles_isDirectHardProcessTauDecayProductFinalState = nullptr;
  std::vector<bool>* genparticles_fromHardProcessBeforeFSR = nullptr;
  std::vector<bool>* genparticles_isLastCopy = nullptr;
  std::vector<bool>* genparticles_isLastCopyBeforeFSR = nullptr;

  std::vector<int>* genparticles_id = nullptr;
  std::vector<int>* genparticles_status = nullptr;

  std::vector<CMSLorentzVector>* genparticles_p4 = nullptr;

  bool allVariablesPresent = (
    (
      !doEventInfo || (
        this->getConsumedValue(_geninfo_processID_, geninfo_processID)
        && this->getConsumedValue(_geninfo_qscale_, geninfo_qscale)
        && this->getConsumedValue(_geninfo_alphaS_, geninfo_alphaS)
        && this->getConsumedValue(_gen_met_, gen_met)
        && this->getConsumedValue(_gen_metPhi_, gen_metPhi)
        )
      )
    &&
    (
      !doParticleInfo || (
        this->getConsumedValue(_genparticles_isPromptFinalState_, genparticles_isPromptFinalState)
        && this->getConsumedValue(_genparticles_isPromptDecayed_, genparticles_isPromptDecayed)
        && this->getConsumedValue(_genparticles_isDirectPromptTauDecayProductFinalState_, genparticles_isDirectPromptTauDecayProductFinalState)
        && this->getConsumedValue(_genparticles_isHardProcess_, genparticles_isHardProcess)
        && this->getConsumedValue(_genparticles_fromHardProcessFinalState_, genparticles_fromHardProcessFinalState)
        && this->getConsumedValue(_genparticles_fromHardProcessDecayed_, genparticles_fromHardProcessDecayed)
        && this->getConsumedValue(_genparticles_isDirectHardProcessTauDecayProductFinalState_, genparticles_isDirectHardProcessTauDecayProductFinalState)
        && this->getConsumedValue(_genparticles_fromHardProcessBeforeFSR_, genparticles_fromHardProcessBeforeFSR)
        && this->getConsumedValue(_genparticles_isLastCopy_, genparticles_isLastCopy)
        && this->getConsumedValue(_genparticles_isLastCopyBeforeFSR_, genparticles_isLastCopyBeforeFSR)

        && this->getConsumedValue(_genparticles_id_, genparticles_id)
        && this->getConsumedValue(_genparticles_status_, genparticles_status)

        && this->getConsumedValue(_genparticles_p4_, genparticles_p4)
        )
      )
    );
  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "JetMETHandler::constructGenInfo: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  // Construct gen. info
  if (doEventInfo){
    geninfo = new GenEventInfo;
    geninfo->processID = geninfo_processID;
    geninfo->qscale = geninfo_qscale;
    geninfo->alphaS = geninfo_alphaS;
    geninfo->genMET = gen_met;
    geninfo->genMETPhi = gen_metPhi;
    geninfo->xsec = SampleHelpers::datasetInfoExtractor.getXsecFromFile(fwktree->getOptions().sampleName(), fwktree->getTag().getRawTag());
  }

  // Construct gen. particles
  if (doParticleInfo){
    const size_t nGenParticles = genparticles_p4->size();
    for (size_t ip=0; ip<nGenParticles; ip++){
      genparticles.push_back(new GenParticleObject(genparticles_id->at(ip), genparticles_p4->at(ip)));
      auto& extras = genparticles.back()->extras;

      extras.isPromptFinalState = genparticles_isPromptFinalState->at(ip);
      extras.isPromptDecayed = genparticles_isPromptDecayed->at(ip);
      extras.isDirectPromptTauDecayProductFinalState = genparticles_isDirectPromptTauDecayProductFinalState->at(ip);
      extras.isHardProcess = genparticles_isHardProcess->at(ip);
      extras.fromHardProcessFinalState = genparticles_fromHardProcessFinalState->at(ip);
      extras.fromHardProcessDecayed = genparticles_fromHardProcessDecayed->at(ip);
      extras.isDirectHardProcessTauDecayProductFinalState = genparticles_isDirectHardProcessTauDecayProductFinalState->at(ip);
      extras.fromHardProcessBeforeFSR = genparticles_fromHardProcessBeforeFSR->at(ip);
      extras.isLastCopy = genparticles_isLastCopy->at(ip);
      extras.isLastCopyBeforeFSR = genparticles_isLastCopyBeforeFSR->at(ip);

      extras.status = genparticles_status->at(ip);
    }
  }

  return true;
}

void GenInfoHandler::bookBranches(BaseTree* tree){
  if (!tree) return;
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(tree);
  if (!fwktree) return;
  if (!fwktree->isMC()) return;

  if (doParticleInfo || doEventInfo) SampleHelpers::setupUsingOptions(fwktree->getOptions());

  if (doEventInfo){
    fwktree->bookEDMBranch<unsigned int>(_geninfo_processID_, 0);
    fwktree->bookEDMBranch<float>(_geninfo_qscale_, 0);
    fwktree->bookEDMBranch<float>(_geninfo_alphaS_, 0);
    fwktree->bookEDMBranch<float>(_gen_met_, 0);
    fwktree->bookEDMBranch<float>(_gen_metPhi_, 0);
  }
  //
  if (doParticleInfo){
    fwktree->bookEDMBranch<std::vector<bool>*>(_genparticles_isPromptFinalState_, nullptr);
    fwktree->bookEDMBranch<std::vector<bool>*>(_genparticles_isPromptDecayed_, nullptr);
    fwktree->bookEDMBranch<std::vector<bool>*>(_genparticles_isDirectPromptTauDecayProductFinalState_, nullptr);
    fwktree->bookEDMBranch<std::vector<bool>*>(_genparticles_isHardProcess_, nullptr);
    fwktree->bookEDMBranch<std::vector<bool>*>(_genparticles_fromHardProcessFinalState_, nullptr);
    fwktree->bookEDMBranch<std::vector<bool>*>(_genparticles_fromHardProcessDecayed_, nullptr);
    fwktree->bookEDMBranch<std::vector<bool>*>(_genparticles_isDirectHardProcessTauDecayProductFinalState_, nullptr);
    fwktree->bookEDMBranch<std::vector<bool>*>(_genparticles_fromHardProcessBeforeFSR_, nullptr);
    fwktree->bookEDMBranch<std::vector<bool>*>(_genparticles_isLastCopy_, nullptr);
    fwktree->bookEDMBranch<std::vector<bool>*>(_genparticles_isLastCopyBeforeFSR_, nullptr);

    fwktree->bookEDMBranch<std::vector<int>*>(_genparticles_id_, nullptr);
    fwktree->bookEDMBranch<std::vector<int>*>(_genparticles_status_, nullptr);

    fwktree->bookEDMBranch<std::vector<CMSLorentzVector>*>(_genparticles_p4_, nullptr);
  }
}
