#include <algorithm>
#include <utility>
#include "Samples.h"
#include "PhotonObject.h"


PhotonVariables::PhotonVariables() :
  rho(0),
  etaSC(0),
  recoChargedHadronIso(0),
  recoNeutralHadronIso(0),
  recoPhotonIso(0),
  sigmaIEtaIEta_full5x5(0),
  hOverE_full5x5(0)
{}
PhotonVariables::PhotonVariables(PhotonVariables const& other) :
  rho(other.rho),
  etaSC(other.etaSC),
  recoChargedHadronIso(other.recoChargedHadronIso),
  recoNeutralHadronIso(other.recoNeutralHadronIso),
  recoPhotonIso(other.recoPhotonIso),
  sigmaIEtaIEta_full5x5(other.sigmaIEtaIEta_full5x5),
  hOverE_full5x5(other.hOverE_full5x5)
{}
void PhotonVariables::swap(PhotonVariables& other){
  std::swap(rho, other.rho);
  std::swap(etaSC, other.etaSC);
  std::swap(recoChargedHadronIso, other.recoChargedHadronIso);
  std::swap(recoNeutralHadronIso, other.recoNeutralHadronIso);
  std::swap(recoPhotonIso, other.recoPhotonIso);
  std::swap(sigmaIEtaIEta_full5x5, other.sigmaIEtaIEta_full5x5);
  std::swap(hOverE_full5x5, other.hOverE_full5x5);
}
PhotonVariables& PhotonVariables::operator=(const PhotonVariables& other){
  PhotonVariables tmp(other);
  swap(tmp);
  return *this;
}


PhotonObject::PhotonObject() :
  ParticleObject(),
  extras()
{}
PhotonObject::PhotonObject(int id_) :
  ParticleObject(id_),
  extras()
{}
PhotonObject::PhotonObject(int id_, CMSLorentzVector const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
PhotonObject::PhotonObject(const PhotonObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void PhotonObject::swap(PhotonObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
PhotonObject& PhotonObject::operator=(const PhotonObject& other){
  PhotonObject tmp(other);
  swap(tmp);
  return *this;
}
