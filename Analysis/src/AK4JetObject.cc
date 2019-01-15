#include <algorithm>
#include <utility>
#include "Samples.h"
#include "AK4JetObject.h"


AK4JetVariables::AK4JetVariables() :
  rho(0),

  npfcands(0),
  parton_flavor(0),
  hadron_flavor(0),
  chargedHadronMultiplicity(0),
  neutralHadronMultiplicity(0),
  photonMultiplicity(0),
  electronMultiplicity(0),
  muonMultiplicity(0),
  chargedMultiplicity(0),
  neutralMultiplicity(0),
  totalMultiplicity(0),

  area(0),
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

  deepCSVb(-1),
  deepCSVc(-1),
  deepCSVl(-1),
  deepCSVbb(-1),
  deepCSVcc(-1),
  pfCombinedInclusiveSecondaryVertexV2BJetTag(-1),
  ptDistribution(0),
  axis1(0),
  axis2(0),

  JEC(1),
  JECup(1),
  JECdn(1),
  JER(1),
  JERup(1),
  JERdn(1)
{}
AK4JetVariables::AK4JetVariables(AK4JetVariables const& other) :
  rho(other.rho),

  npfcands(other.npfcands),
  parton_flavor(other.parton_flavor),
  hadron_flavor(other.hadron_flavor),
  chargedHadronMultiplicity(other.chargedHadronMultiplicity),
  neutralHadronMultiplicity(other.neutralHadronMultiplicity),
  photonMultiplicity(other.photonMultiplicity),
  electronMultiplicity(other.electronMultiplicity),
  muonMultiplicity(other.muonMultiplicity),
  chargedMultiplicity(other.chargedMultiplicity),
  neutralMultiplicity(other.neutralMultiplicity),
  totalMultiplicity(other.totalMultiplicity),

  area(other.area),
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

  deepCSVb(other.deepCSVb),
  deepCSVc(other.deepCSVc),
  deepCSVl(other.deepCSVl),
  deepCSVbb(other.deepCSVbb),
  deepCSVcc(other.deepCSVcc),
  pfCombinedInclusiveSecondaryVertexV2BJetTag(other.pfCombinedInclusiveSecondaryVertexV2BJetTag),
  ptDistribution(other.ptDistribution),
  axis1(other.axis1),
  axis2(other.axis2),

  JEC(other.JEC),
  JECup(other.JECup),
  JECdn(other.JECdn),
  JER(other.JER),
  JERup(other.JERup),
  JERdn(other.JERdn)
{}
void AK4JetVariables::swap(AK4JetVariables& other){
  std::swap(rho, other.rho);

  std::swap(npfcands, other.npfcands);
  std::swap(parton_flavor, other.parton_flavor);
  std::swap(hadron_flavor, other.hadron_flavor);
  std::swap(chargedHadronMultiplicity, other.chargedHadronMultiplicity);
  std::swap(neutralHadronMultiplicity, other.neutralHadronMultiplicity);
  std::swap(photonMultiplicity, other.photonMultiplicity);
  std::swap(electronMultiplicity, other.electronMultiplicity);
  std::swap(muonMultiplicity, other.muonMultiplicity);
  std::swap(chargedMultiplicity, other.chargedMultiplicity);
  std::swap(neutralMultiplicity, other.neutralMultiplicity);
  std::swap(totalMultiplicity, other.totalMultiplicity);

  std::swap(area, other.area);
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

  std::swap(deepCSVb, other.deepCSVb);
  std::swap(deepCSVc, other.deepCSVc);
  std::swap(deepCSVl, other.deepCSVl);
  std::swap(deepCSVbb, other.deepCSVbb);
  std::swap(deepCSVcc, other.deepCSVcc);
  std::swap(pfCombinedInclusiveSecondaryVertexV2BJetTag, other.pfCombinedInclusiveSecondaryVertexV2BJetTag);
  std::swap(ptDistribution, other.ptDistribution);
  std::swap(axis1, other.axis1);
  std::swap(axis2, other.axis2);

  std::swap(JEC, other.JEC);
  std::swap(JECup, other.JECup);
  std::swap(JECdn, other.JECdn);
  std::swap(JER, other.JER);
  std::swap(JERup, other.JERup);
  std::swap(JERdn, other.JERdn);
}
AK4JetVariables& AK4JetVariables::operator=(const AK4JetVariables& other){
  AK4JetVariables tmp(other);
  swap(tmp);
  return *this;
}


AK4JetObject::AK4JetObject() :
  ParticleObject(),
  extras(),
  associatedGenJet(nullptr)
{}
AK4JetObject::AK4JetObject(int id_) :
  ParticleObject(id_),
  extras(),
  associatedGenJet(nullptr)
{}
AK4JetObject::AK4JetObject(int id_, CMSLorentzVector momentum_) :
  ParticleObject(id_, momentum_),
  extras(),
  associatedGenJet(nullptr)
{}
AK4JetObject::AK4JetObject(const AK4JetObject& other) :
  ParticleObject(other),
  extras(other.extras),
  associatedGenJet(other.associatedGenJet)
{}
void AK4JetObject::swap(AK4JetObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  std::swap(associatedGenJet, other.associatedGenJet);
  extras.swap(other.extras);
}
AK4JetObject& AK4JetObject::operator=(const AK4JetObject& other){
  AK4JetObject tmp(other);
  swap(tmp);
  return *this;
}
AK4JetObject::~AK4JetObject(){}

CMSLorentzVector AK4JetObject::getCorrectedMomentum(int icorr) const{
  const float& JEC = extras.JEC;
  const float& JER = extras.JER;
  const float& JECup = extras.JECup;
  const float& JERup = extras.JERup;
  const float& JECdn = extras.JECdn;
  const float& JERdn = extras.JERdn;
  switch (icorr){
  case 1:
    return momentum*JECup*JER;
  case -1:
    return momentum*JECdn*JER;
  case 2:
    return momentum*JEC*JERup;
  case -2:
    return momentum*JEC*JERdn;
  default:
    return momentum*JEC*JER;
  }
}
