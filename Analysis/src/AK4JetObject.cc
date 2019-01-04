#include <algorithm>
#include <utility>
#include "Samples.h"
#include "AK4JetObject.h"


AK4JetVariables::AK4JetVariables() :
  rho(0),

  npfcands(0),
  chargedMultiplicity(0),
  neutralMultiplicity(0),

  undoJEC(0),
  chargedHadronE(0),
  chargedEmE(0),
  neutralHadronE(0),
  neutralEmE(0),
  muonE(0)
{}
AK4JetVariables::AK4JetVariables(AK4JetVariables const& other) :
  rho(other.rho),

  npfcands(other.npfcands),
  chargedMultiplicity(other.chargedMultiplicity),
  neutralMultiplicity(other.neutralMultiplicity),

  undoJEC(other.undoJEC),
  chargedHadronE(other.chargedHadronE),
  chargedEmE(other.chargedEmE),
  neutralHadronE(other.neutralHadronE),
  neutralEmE(other.neutralEmE),
  muonE(other.muonE)
{}
void AK4JetVariables::swap(AK4JetVariables& other){
  std::swap(rho, other.rho);

  std::swap(npfcands, other.npfcands);
  std::swap(chargedMultiplicity, other.chargedMultiplicity);
  std::swap(neutralMultiplicity, other.neutralMultiplicity);

  std::swap(undoJEC, other.undoJEC);
  std::swap(chargedHadronE, other.chargedHadronE);
  std::swap(chargedEmE, other.chargedEmE);
  std::swap(neutralHadronE, other.neutralHadronE);
  std::swap(neutralEmE, other.neutralEmE);
  std::swap(muonE, other.muonE);
}
AK4JetVariables& AK4JetVariables::operator=(const AK4JetVariables& other){
  AK4JetVariables tmp(other);
  swap(tmp);
  return *this;
}


AK4JetObject::AK4JetObject() :
  ParticleObject(),
  extras()
{}
AK4JetObject::AK4JetObject(int id_) :
  ParticleObject(id_),
  extras()
{}
AK4JetObject::AK4JetObject(int id_, CMSLorentzVector momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
AK4JetObject::AK4JetObject(const AK4JetObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void AK4JetObject::swap(AK4JetObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
AK4JetObject& AK4JetObject::operator=(const AK4JetObject& other){
  AK4JetObject tmp(other);
  swap(tmp);
  return *this;
}
AK4JetObject::~AK4JetObject(){}
