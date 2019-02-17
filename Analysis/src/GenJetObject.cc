#include <algorithm>
#include <utility>
#include "GenJetObject.h"


GenJetObject::GenJetObject() : ParticleObject()
{}
GenJetObject::GenJetObject(CMSLorentzVector const& momentum_) : ParticleObject(0, momentum_)
{}
GenJetObject::GenJetObject(const GenJetObject& other) :
  ParticleObject(other)
{}
void GenJetObject::swap(GenJetObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
}
GenJetObject& GenJetObject::operator=(const GenJetObject& other){
  GenJetObject tmp(other);
  swap(tmp);
  return *this;
}
GenJetObject::~GenJetObject(){}
