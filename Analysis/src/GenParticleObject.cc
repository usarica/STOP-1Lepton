#include <algorithm>
#include <utility>
#include "Samples.h"
#include "GenParticleObject.h"


GenParticleVariables::GenParticleVariables() :
  isPromptFinalState(nullptr),
  isPromptDecayed(nullptr),
  isDirectPromptTauDecayProductFinalState(nullptr),
  isHardProcess(nullptr),
  fromHardProcessFinalState(nullptr),
  fromHardProcessDecayed(nullptr),
  isDirectHardProcessTauDecayProductFinalState(nullptr),
  fromHardProcessBeforeFSR(nullptr),
  isLastCopy(nullptr),
  isLastCopyBeforeFSR(nullptr),
  status(-9000)
{}
GenParticleVariables::GenParticleVariables(GenParticleVariables const& other) :
  isPromptFinalState(other.isPromptFinalState),
  isPromptDecayed(other.isPromptDecayed),
  isDirectPromptTauDecayProductFinalState(other.isDirectPromptTauDecayProductFinalState),
  isHardProcess(other.isHardProcess),
  fromHardProcessFinalState(other.fromHardProcessFinalState),
  fromHardProcessDecayed(other.fromHardProcessDecayed),
  isDirectHardProcessTauDecayProductFinalState(other.isDirectHardProcessTauDecayProductFinalState),
  fromHardProcessBeforeFSR(other.fromHardProcessBeforeFSR),
  isLastCopy(other.isLastCopy),
  isLastCopyBeforeFSR(other.isLastCopyBeforeFSR),
  status(other.status)
{}
void GenParticleVariables::swap(GenParticleVariables& other){
  std::swap(isPromptFinalState, other.isPromptFinalState);
  std::swap(isPromptDecayed, other.isPromptDecayed);
  std::swap(isDirectPromptTauDecayProductFinalState, other.isDirectPromptTauDecayProductFinalState);
  std::swap(isHardProcess, other.isHardProcess);
  std::swap(fromHardProcessFinalState, other.fromHardProcessFinalState);
  std::swap(fromHardProcessDecayed, other.fromHardProcessDecayed);
  std::swap(isDirectHardProcessTauDecayProductFinalState, other.isDirectHardProcessTauDecayProductFinalState);
  std::swap(fromHardProcessBeforeFSR, other.fromHardProcessBeforeFSR);
  std::swap(isLastCopy, other.isLastCopy);
  std::swap(isLastCopyBeforeFSR, other.isLastCopyBeforeFSR);
  std::swap(status, other.status);
}
GenParticleVariables& GenParticleVariables::operator=(const GenParticleVariables& other){
  GenParticleVariables tmp(other);
  swap(tmp);
  return *this;
}


GenParticleObject::GenParticleObject() :
  ParticleObject(),
  extras()
{}
GenParticleObject::GenParticleObject(int id_) :
  ParticleObject(id_),
  extras()
{}
GenParticleObject::GenParticleObject(int id_, CMSLorentzVector const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
GenParticleObject::GenParticleObject(const GenParticleObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void GenParticleObject::swap(GenParticleObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
GenParticleObject& GenParticleObject::operator=(const GenParticleObject& other){
  GenParticleObject tmp(other);
  swap(tmp);
  return *this;
}
GenParticleObject::~GenParticleObject(){}
