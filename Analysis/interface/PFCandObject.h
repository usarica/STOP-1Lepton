#ifndef PFCANDOBJECT_H
#define PFCANDOBJECT_H

#include "ParticleObject.h"


class PFCandVariables{
public:
  bool trackHighPurity;

  int charge;

  float dxy;
  float dz;
  float dxyError;
  float dzError;

  PFCandVariables();
  PFCandVariables(PFCandVariables const& other);
  PFCandVariables& operator=(const PFCandVariables& other);

  void swap(PFCandVariables& other);

};

class PFCandObject : public ParticleObject{
public:
  PFCandVariables extras;

  PFCandObject();
  PFCandObject(int id_);
  PFCandObject(int id_, CMSLorentzVector const& mom_);
  PFCandObject(const PFCandObject& other);
  PFCandObject& operator=(const PFCandObject& other);
  ~PFCandObject();

  void swap(PFCandObject& other);

};

#endif
