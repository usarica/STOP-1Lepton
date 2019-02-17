#ifndef GENJETOBJECT_H
#define GENJETOBJECT_H

#include "ParticleObject.h"


class GenJetObject : public ParticleObject{
public:
  GenJetObject();
  GenJetObject(CMSLorentzVector const& mom_);
  GenJetObject(const GenJetObject& other);
  GenJetObject& operator=(const GenJetObject& other);
  ~GenJetObject();

  void swap(GenJetObject& other);

};

#endif
