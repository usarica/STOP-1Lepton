#ifndef PHOTONOBJECT_H
#define PHOTONOBJECT_H

#include "ParticleObject.h"


class PhotonVariables{
public:
  float rho;
  float etaSC; // Supercluster eta
  float recoChargedHadronIso;
  float recoNeutralHadronIso;
  float recoPhotonIso;
  float sigmaIEtaIEta_full5x5;
  float hOverE_full5x5;

  PhotonVariables();
  PhotonVariables(PhotonVariables const& other);
  PhotonVariables& operator=(const PhotonVariables& other);

  void swap(PhotonVariables& other);

};

class PhotonObject : public ParticleObject{
public:
  PhotonVariables extras;

  PhotonObject();
  PhotonObject(int id_);
  PhotonObject(int id_, CMSLorentzVector mom_);
  PhotonObject(const PhotonObject& other);
  PhotonObject& operator=(const PhotonObject& other);
  ~PhotonObject(){}

  void swap(PhotonObject& other);

};

#endif
