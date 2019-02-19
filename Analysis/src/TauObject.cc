#include <algorithm>
#include <utility>
#include "TauObject.h"


TauVariables::TauVariables() :
  charge(0),

  pfDecayModeFinding(false),
  pfIso(false)
{}
TauVariables::TauVariables(TauVariables const& other) :
  charge(other.charge),
  
  pfDecayModeFinding(other.pfDecayModeFinding),
  pfIso(other.pfIso)
{}
void TauVariables::swap(TauVariables& other){
  std::swap(charge, other.charge);

  std::swap(pfDecayModeFinding, other.pfDecayModeFinding);
  std::swap(pfIso, other.pfIso);
}
TauVariables& TauVariables::operator=(const TauVariables& other){
  TauVariables tmp(other);
  swap(tmp);
  return *this;
}


TauObject::TauObject() :
  ParticleObject(),
  extras()
{}
TauObject::TauObject(int id_, CMSLorentzVector const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
TauObject::TauObject(const TauObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void TauObject::swap(TauObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
TauObject& TauObject::operator=(const TauObject& other){
  TauObject tmp(other);
  swap(tmp);
  return *this;
}
TauObject::~TauObject(){}
