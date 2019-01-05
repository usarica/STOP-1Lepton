#ifndef AK4JETOBJECT_H
#define AK4JETOBJECT_H

#include "ParticleObject.h"


class AK4JetVariables{
public:
  float rho;

  int npfcands;
  int chargedHadronMultiplicity;
  int neutralHadronMultiplicity;
  int photonMultiplicity;
  int electronMultiplicity;
  int muonMultiplicity;
  int chargedMultiplicity;
  int neutralMultiplicity;
  int totalMultiplicity;

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

  float pfCombinedInclusiveSecondaryVertexV2BJetTag;
  float ptDistribution;
  float axis1;
  float axis2;

  AK4JetVariables();
  AK4JetVariables(AK4JetVariables const& other);
  AK4JetVariables& operator=(const AK4JetVariables& other);

  void swap(AK4JetVariables& other);

};

class AK4JetObject : public ParticleObject{
public:
  AK4JetVariables extras;

  AK4JetObject();
  AK4JetObject(int id_);
  AK4JetObject(int id_, CMSLorentzVector mom_);
  AK4JetObject(const AK4JetObject& other);
  AK4JetObject& operator=(const AK4JetObject& other);
  ~AK4JetObject();

  void swap(AK4JetObject& other);

};

#endif
