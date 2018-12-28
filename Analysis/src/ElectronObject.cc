#include <algorithm>
#include <utility>
#include "Samples.h"
#include "ElectronObject.h"
#include "PDGHelpers.h"


using namespace PDGHelpers;


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
  id(-9000),
  momentum(0, 0, 0, 0),
  extras()
{}
ElectronObject::ElectronObject(int id_) :
  id(id_),
  momentum(0, 0, 0, 0),
  extras()
{}
ElectronObject::ElectronObject(int id_, CMSLorentzVector momentum_) :
  id(id_),
  momentum(momentum_),
  extras()
{}
ElectronObject::ElectronObject(const ElectronObject& other) :
  id(other.id),
  momentum(other.momentum),
  extras(other.extras)
{}
void ElectronObject::swap(ElectronObject& other){
  std::swap(id, other.id);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
ElectronObject& ElectronObject::operator=(const ElectronObject& other){
  ElectronObject tmp(other);
  swap(tmp);
  return *this;
}
ElectronObject::~ElectronObject(){}

float ElectronObject::EinvMinusPinv()const{
  float const& ecalEnergy = extras.ecalEnergy;
  float const& EoverP = extras.eOverPIn;
  return (1.f-EoverP)/ecalEnergy;
}
float ElectronObject::charge()const{
  float cpos=0;
  if (isAWBoson(id) || abs(id)==37 || abs(id)==2212 || abs(id)==211 || abs(id)==321 || abs(id)==411 || abs(id)==521) cpos = 1.;
  else if (isALepton(id)) cpos = -1.;
  else if (isUpTypeQuark(id)) cpos = 2./3.;
  else if (isDownTypeQuark(id)) cpos = -1./3.;
  if (id<0) cpos *= -1.;
  return cpos;
}
