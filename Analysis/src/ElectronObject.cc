#include <algorithm>
#include <utility>
#include "Samples.h"
#include "ElectronObject.h"


ElectronVariables::ElectronVariables() :
  conv_vtx_flag(false),
  expectedMissingInnerHits(0),
  energySC(0),
  etaSC(0),
  etaSeedSC(0),
  rho(0),
  sigmaIEtaIEta_full5x5(0),
  dEtaIn(0),
  dPhiIn(0),
  hOverE(0),
  ecalEnergy(0),
  eOverPIn(0),
  dxyPV(0),
  dzPV(0),
  miniIso_ch(0),
  miniIso_nh(0),
  miniIso_em(0)
{}
ElectronVariables::ElectronVariables(ElectronVariables const& other) :
  conv_vtx_flag(other.conv_vtx_flag),
  expectedMissingInnerHits(other.expectedMissingInnerHits),
  energySC(other.energySC),
  etaSC(other.etaSC),
  etaSeedSC(other.etaSeedSC),
  rho(other.rho),
  sigmaIEtaIEta_full5x5(other.sigmaIEtaIEta_full5x5),
  dEtaIn(other.dEtaIn),
  dPhiIn(other.dPhiIn),
  hOverE(other.hOverE),
  ecalEnergy(other.ecalEnergy),
  eOverPIn(other.eOverPIn),
  dxyPV(other.dxyPV),
  dzPV(other.dzPV),
  miniIso_ch(other.miniIso_ch),
  miniIso_nh(other.miniIso_nh),
  miniIso_em(other.miniIso_em)
{}
void ElectronVariables::swap(ElectronVariables& other){
  std::swap(conv_vtx_flag, other.conv_vtx_flag);
  std::swap(expectedMissingInnerHits, other.expectedMissingInnerHits);
  std::swap(energySC, other.energySC);
  std::swap(etaSC, other.etaSC);
  std::swap(etaSeedSC, other.etaSeedSC);
  std::swap(rho, other.rho);
  std::swap(sigmaIEtaIEta_full5x5, other.sigmaIEtaIEta_full5x5);
  std::swap(dEtaIn, other.dEtaIn);
  std::swap(dPhiIn, other.dPhiIn);
  std::swap(hOverE, other.hOverE);
  std::swap(ecalEnergy, other.ecalEnergy);
  std::swap(eOverPIn, other.eOverPIn);
  std::swap(dxyPV, other.dxyPV);
  std::swap(dzPV, other.dzPV);
  std::swap(miniIso_ch, other.miniIso_ch);
  std::swap(miniIso_nh, other.miniIso_nh);
  std::swap(miniIso_em, other.miniIso_em);
}
ElectronVariables& ElectronVariables::operator=(const ElectronVariables& other){
  ElectronVariables tmp(other);
  swap(tmp);
  return *this;
}


ElectronObject::ElectronObject() :
  ParticleObject(),
  extras()
{}
ElectronObject::ElectronObject(int id_) :
  ParticleObject(id_),
  extras()
{}
ElectronObject::ElectronObject(int id_, CMSLorentzVector momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
ElectronObject::ElectronObject(const ElectronObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void ElectronObject::swap(ElectronObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
ElectronObject& ElectronObject::operator=(const ElectronObject& other){
  ElectronObject tmp(other);
  swap(tmp);
  return *this;
}
ElectronObject::~ElectronObject(){}

float ElectronObject::EinvMinusPinv() const{
  float const& ecalEnergy = extras.ecalEnergy;
  float const& EoverP = extras.eOverPIn;
  return (1.f-EoverP)/ecalEnergy;
}
