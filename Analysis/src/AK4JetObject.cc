#include <algorithm>
#include <utility>
#include "Samples.h"
#include "AK4JetObject.h"


AK4JetVariables::AK4JetVariables() :
  rho(0),

  npfcands(0),
  chargedHadronMultiplicity(0),
  neutralHadronMultiplicity(0),
  photonMultiplicity(0),
  electronMultiplicity(0),
  muonMultiplicity(0),
  chargedMultiplicity(0),
  neutralMultiplicity(0),
  totalMultiplicity(0),

  undoJEC(0),
  chargedHadronE(0),
  chargedEmE(0),
  neutralHadronE(0),
  neutralEmE(0),
  hfHadronE(0),
  hfEmE(0),
  photonE(0),
  electronE(0),
  muonE(0),

  pfCombinedInclusiveSecondaryVertexV2BJetTag(0),
  ptDistribution(0),
  axis1(0),
  axis2(0)
{}
AK4JetVariables::AK4JetVariables(AK4JetVariables const& other) :
  rho(other.rho),

  npfcands(other.npfcands),
  chargedHadronMultiplicity(other.chargedHadronMultiplicity),
  neutralHadronMultiplicity(other.neutralHadronMultiplicity),
  photonMultiplicity(other.photonMultiplicity),
  electronMultiplicity(other.electronMultiplicity),
  muonMultiplicity(other.muonMultiplicity),
  chargedMultiplicity(other.chargedMultiplicity),
  neutralMultiplicity(other.neutralMultiplicity),
  totalMultiplicity(other.totalMultiplicity),

  undoJEC(other.undoJEC),
  chargedHadronE(other.chargedHadronE),
  chargedEmE(other.chargedEmE),
  neutralHadronE(other.neutralHadronE),
  neutralEmE(other.neutralEmE),
  hfHadronE(other.hfHadronE),
  hfEmE(other.hfEmE),
  photonE(other.photonE),
  electronE(other.electronE),
  muonE(other.muonE),

  pfCombinedInclusiveSecondaryVertexV2BJetTag(other.pfCombinedInclusiveSecondaryVertexV2BJetTag),
  ptDistribution(other.ptDistribution),
  axis1(other.axis1),
  axis2(other.axis2)
{}
void AK4JetVariables::swap(AK4JetVariables& other){
  std::swap(rho, other.rho);

  std::swap(npfcands, other.npfcands);
  std::swap(chargedHadronMultiplicity, other.chargedHadronMultiplicity);
  std::swap(neutralHadronMultiplicity, other.neutralHadronMultiplicity);
  std::swap(photonMultiplicity, other.photonMultiplicity);
  std::swap(electronMultiplicity, other.electronMultiplicity);
  std::swap(muonMultiplicity, other.muonMultiplicity);
  std::swap(chargedMultiplicity, other.chargedMultiplicity);
  std::swap(neutralMultiplicity, other.neutralMultiplicity);
  std::swap(totalMultiplicity, other.totalMultiplicity);

  std::swap(undoJEC, other.undoJEC);
  std::swap(chargedHadronE, other.chargedHadronE);
  std::swap(chargedEmE, other.chargedEmE);
  std::swap(neutralHadronE, other.neutralHadronE);
  std::swap(neutralEmE, other.neutralEmE);
  std::swap(hfHadronE, other.hfHadronE);
  std::swap(hfEmE, other.hfEmE);
  std::swap(photonE, other.photonE);
  std::swap(electronE, other.electronE);
  std::swap(muonE, other.muonE);

  std::swap(pfCombinedInclusiveSecondaryVertexV2BJetTag, other.pfCombinedInclusiveSecondaryVertexV2BJetTag);
  std::swap(ptDistribution, other.ptDistribution);
  std::swap(axis1, other.axis1);
  std::swap(axis2, other.axis2);
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
