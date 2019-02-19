#ifndef TAUOBJECT_H
#define TAUOBJECT_H

#include "ParticleObject.h"


class TauVariables{
public:
  int charge;

  bool pfDecayModeFinding;
  bool pfIso;
  
  TauVariables();
  TauVariables(TauVariables const& other);
  TauVariables& operator=(const TauVariables& other);

  void swap(TauVariables& other);

};

class TauObject : public ParticleObject{
public:
  TauVariables extras;

  TauObject();
  TauObject(int id_, CMSLorentzVector const& mom_);
  TauObject(const TauObject& other);
  TauObject& operator=(const TauObject& other);
  ~TauObject();

  void swap(TauObject& other);

};

#endif
