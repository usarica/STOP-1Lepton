#include <algorithm>
#include <utility>
#include "IsoTrackObject.h"


IsoTrackVariables::IsoTrackVariables() :
  isPFCand(false),
  hasLepOverlap(false),

  charge(0),

  pfIso_ch(0),
  dz(0)
{}
IsoTrackVariables::IsoTrackVariables(IsoTrackVariables const& other) :
  isPFCand(other.isPFCand),
  hasLepOverlap(other.hasLepOverlap),

  charge(other.charge),
  
  pfIso_ch(other.pfIso_ch),
  dz(other.dz)
{}
void IsoTrackVariables::swap(IsoTrackVariables& other){
  std::swap(isPFCand, other.isPFCand);
  std::swap(hasLepOverlap, other.hasLepOverlap);

  std::swap(charge, other.charge);

  std::swap(pfIso_ch, other.pfIso_ch);
  std::swap(dz, other.dz);
}
IsoTrackVariables& IsoTrackVariables::operator=(const IsoTrackVariables& other){
  IsoTrackVariables tmp(other);
  swap(tmp);
  return *this;
}


IsoTrackObject::IsoTrackObject() :
  ParticleObject(),
  extras()
{}
IsoTrackObject::IsoTrackObject(int id_, CMSLorentzVector const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
IsoTrackObject::IsoTrackObject(const IsoTrackObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void IsoTrackObject::swap(IsoTrackObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
IsoTrackObject& IsoTrackObject::operator=(const IsoTrackObject& other){
  IsoTrackObject tmp(other);
  swap(tmp);
  return *this;
}
IsoTrackObject::~IsoTrackObject(){}
