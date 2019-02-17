#include <algorithm>
#include <utility>
#include "PFCandObject.h"


PFCandVariables::PFCandVariables() :
  trackHighPurity(0),

  charge(0),

  dxy(0),
  dz(0),
  dxyError(0),
  dzError(0)
{}
PFCandVariables::PFCandVariables(PFCandVariables const& other) :
  trackHighPurity(other.trackHighPurity),

  charge(other.charge),
  
  dxy(other.dxy),
  dz(other.dz),
  dxyError(other.dxyError),
  dzError(other.dzError)
{}
void PFCandVariables::swap(PFCandVariables& other){
  std::swap(trackHighPurity, other.trackHighPurity);

  std::swap(charge, other.charge);

  std::swap(dxy, other.dxy);
  std::swap(dz, other.dz);
  std::swap(dxyError, other.dxyError);
  std::swap(dzError, other.dzError);
}
PFCandVariables& PFCandVariables::operator=(const PFCandVariables& other){
  PFCandVariables tmp(other);
  swap(tmp);
  return *this;
}


PFCandObject::PFCandObject() :
  ParticleObject(),
  extras()
{}
PFCandObject::PFCandObject(int id_) :
  ParticleObject(id_),
  extras()
{}
PFCandObject::PFCandObject(int id_, CMSLorentzVector const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
PFCandObject::PFCandObject(const PFCandObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void PFCandObject::swap(PFCandObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
PFCandObject& PFCandObject::operator=(const PFCandObject& other){
  PFCandObject tmp(other);
  swap(tmp);
  return *this;
}
PFCandObject::~PFCandObject(){}
