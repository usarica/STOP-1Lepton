#include <algorithm>
#include <utility>
#include "Samples.h"
#include "AK8JetObject.h"


AK8JetVariables::AK8JetVariables() :
  rho(0),

  parton_flavor(0),

  area(0),
  undoJEC(0),
  tau1(0),
  tau2(0),
  tau3(0),
  deepdisc_qcd(0),
  deepdisc_top(0),
  deepdisc_w(0),
  deepdisc_z(0),
  deepdisc_zbb(0),
  deepdisc_hbb(0),
  deepdisc_h4q(0),

  JEC(1),
  JECup(1),
  JECdn(1),

  estimatedPtResolution(0),
  JER(1),
  JERup(1),
  JERdn(1)
{}
AK8JetVariables::AK8JetVariables(AK8JetVariables const& other) :
  rho(other.rho),

  parton_flavor(other.parton_flavor),

  area(other.area),
  undoJEC(other.undoJEC),
  tau1(other.tau1),
  tau2(other.tau2),
  tau3(other.tau3),
  deepdisc_qcd(other.deepdisc_qcd),
  deepdisc_top(other.deepdisc_top),
  deepdisc_w(other.deepdisc_w),
  deepdisc_z(other.deepdisc_z),
  deepdisc_zbb(other.deepdisc_zbb),
  deepdisc_hbb(other.deepdisc_hbb),
  deepdisc_h4q(other.deepdisc_h4q),

  JEC(other.JEC),
  JECup(other.JECup),
  JECdn(other.JECdn),

  estimatedPtResolution(other.estimatedPtResolution),
  JER(other.JER),
  JERup(other.JERup),
  JERdn(other.JERdn)
{}
void AK8JetVariables::swap(AK8JetVariables& other){
  std::swap(rho, other.rho);

  std::swap(parton_flavor, other.parton_flavor);

  std::swap(area, other.area);
  std::swap(undoJEC, other.undoJEC);
  std::swap(tau1, other.tau1);
  std::swap(tau2, other.tau2);
  std::swap(tau3, other.tau3);
  std::swap(deepdisc_qcd, other.deepdisc_qcd);
  std::swap(deepdisc_top, other.deepdisc_top);
  std::swap(deepdisc_w, other.deepdisc_w);
  std::swap(deepdisc_z, other.deepdisc_z);
  std::swap(deepdisc_zbb, other.deepdisc_zbb);
  std::swap(deepdisc_hbb, other.deepdisc_hbb);
  std::swap(deepdisc_h4q, other.deepdisc_h4q);

  std::swap(JEC, other.JEC);
  std::swap(JECup, other.JECup);
  std::swap(JECdn, other.JECdn);

  std::swap(estimatedPtResolution, other.estimatedPtResolution);
  std::swap(JER, other.JER);
  std::swap(JERup, other.JERup);
  std::swap(JERdn, other.JERdn);
}
AK8JetVariables& AK8JetVariables::operator=(const AK8JetVariables& other){
  AK8JetVariables tmp(other);
  swap(tmp);
  return *this;
}


AK8JetObject::AK8JetObject() :
  ParticleObject(),
  extras(),
  associatedGenJet(nullptr)
{}
AK8JetObject::AK8JetObject(int id_) :
  ParticleObject(id_),
  extras(),
  associatedGenJet(nullptr)
{}
AK8JetObject::AK8JetObject(int id_, CMSLorentzVector momentum_) :
  ParticleObject(id_, momentum_),
  extras(),
  associatedGenJet(nullptr)
{}
AK8JetObject::AK8JetObject(const AK8JetObject& other) :
  ParticleObject(other),
  extras(other.extras),
  associatedGenJet(other.associatedGenJet)
{}
void AK8JetObject::swap(AK8JetObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  std::swap(associatedGenJet, other.associatedGenJet);
  extras.swap(other.extras);
}
AK8JetObject& AK8JetObject::operator=(const AK8JetObject& other){
  AK8JetObject tmp(other);
  swap(tmp);
  return *this;
}
AK8JetObject::~AK8JetObject(){}

CMSLorentzVector AK8JetObject::getCorrectedMomentum(int icorr) const{
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
