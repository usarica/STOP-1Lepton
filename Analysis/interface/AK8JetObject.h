#ifndef AK8JETOBJECT_H
#define AK8JETOBJECT_H

#include "ParticleObject.h"
#include "GenJetObject.h"


class AK8JetVariables{
public:
  float rho;

  int parton_flavor;

  float area;
  float undoJEC;
  float tau1;
  float tau2;
  float tau3;
  float deepdisc_qcd;
  float deepdisc_top;
  float deepdisc_w;
  float deepdisc_z;
  float deepdisc_zbb;
  float deepdisc_hbb;
  float deepdisc_h4q;

  float JEC;
  float JECup;
  float JECdn;

  float estimatedPtResolution;
  float JER;
  float JERup;
  float JERdn;

  AK8JetVariables();
  AK8JetVariables(AK8JetVariables const& other);
  AK8JetVariables& operator=(const AK8JetVariables& other);

  void swap(AK8JetVariables& other);

};

class AK8JetObject : public ParticleObject{
public:
  constexpr static float ConeRadiusConstant = 0.8;
  AK8JetVariables extras;
  GenJetObject* associatedGenJet;

  AK8JetObject();
  AK8JetObject(int id_);
  AK8JetObject(int id_, CMSLorentzVector mom_);
  AK8JetObject(const AK8JetObject& other);
  AK8JetObject& operator=(const AK8JetObject& other);
  ~AK8JetObject();

  void swap(AK8JetObject& other);

  CMSLorentzVector getCorrectedMomentum(int icorr) const; // icorr = 0 for nominal, +-1 for JEC up/dn, +-1 for JER up/dn
  CMSLorentzVector getFinalMomentum() const{ return getCorrectedMomentum(0); }

};

#endif
