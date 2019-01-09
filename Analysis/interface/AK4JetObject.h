#ifndef AK4JETOBJECT_H
#define AK4JETOBJECT_H

#include "ParticleObject.h"


class AK4JetVariables{
public:
  float rho;

  int npfcands;
  int parton_flavor;
  int hadron_flavor;
  int chargedHadronMultiplicity;
  int neutralHadronMultiplicity;
  int photonMultiplicity;
  int electronMultiplicity;
  int muonMultiplicity;
  int chargedMultiplicity;
  int neutralMultiplicity;
  int totalMultiplicity;

  float area;
  float undoJEC;
  float chargedHadronE;
  float chargedEmE;
  float neutralHadronE;
  float neutralEmE;
  float hfHadronE;
  float hfEmE;
  float photonE;
  float electronE;
  float muonE;

  float deepCSVb;
  float deepCSVc;
  float deepCSVl;
  float deepCSVbb;
  float deepCSVcc;
  float pfCombinedInclusiveSecondaryVertexV2BJetTag; // CSV b-tag
  float ptDistribution;
  float axis1;
  float axis2;

  float JEC;
  float JECup;
  float JECdn;
  float JER;
  float JERup;
  float JERdn;

  AK4JetVariables();
  AK4JetVariables(AK4JetVariables const& other);
  AK4JetVariables& operator=(const AK4JetVariables& other);

  void swap(AK4JetVariables& other);

};

class AK4JetObject : public ParticleObject{
public:
  constexpr static float ConeRadiusConstant = 0.4;
  AK4JetVariables extras;

  AK4JetObject();
  AK4JetObject(int id_);
  AK4JetObject(int id_, CMSLorentzVector mom_);
  AK4JetObject(const AK4JetObject& other);
  AK4JetObject& operator=(const AK4JetObject& other);
  ~AK4JetObject();

  void swap(AK4JetObject& other);

  CMSLorentzVector getCorrectedMomentum(int icorr) const; // icorr = 0 for nominal, +-1 for JEC up/dn, +-1 for JER up/dn
  CMSLorentzVector getFinalMomentum() const{ return getCorrectedMomentum(0); }

};

#endif
