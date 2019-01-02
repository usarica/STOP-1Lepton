#ifndef ELECTRONOBJECT_H
#define ELECTRONOBJECT_H

#include "ParticleObject.h"


class ElectronVariables{
public:
  bool conv_vtx_flag; // Conversion veto

  int expectedMissingInnerHits; // expectedMissingInnerHits

  float energySC; // Supercluster energy, eSC
  float etaSC; // Supercluster eta
  float etaSeedSC; // Supercluster seed eta
  float rho; // evt_fixgridfastjet_all_rho, actually a constant over the event
  float sigmaIEtaIEta_full5x5; // full5x5_sigmaIetaIeta
  float dEtaIn; // abs(dEtaIn)
  float dPhiIn; // abs(dPhiIn)
  float hOverE; // hOverE
  float ecalEnergy;
  float eOverPIn;
  float dxyPV; // abs(d0)
  float dzPV; // abs(dz)
  float miniIso_ch; // Charged min. iso.
  float miniIso_nh; // Neutral min. iso.
  float miniIso_em; // EM min. iso.

  ElectronVariables();
  ElectronVariables(ElectronVariables const& other);
  ElectronVariables& operator=(const ElectronVariables& other);

  void swap(ElectronVariables& other);

};

class ElectronObject : public ParticleObject{
public:
  ElectronVariables extras;

  ElectronObject();
  ElectronObject(int id_);
  ElectronObject(int id_, CMSLorentzVector mom_);
  ElectronObject(const ElectronObject& other);
  ElectronObject& operator=(const ElectronObject& other);
  ~ElectronObject();

  void swap(ElectronObject& other);

  float EinvMinusPinv() const;

};

#endif
