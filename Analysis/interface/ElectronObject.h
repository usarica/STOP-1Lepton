#ifndef ELECTRONOBJECT_H
#define ELECTRONOBJECT_H

#include "MELAParticle.h"
#include "SimpleEntry.h"


class ElectronObject : public MELAParticle{
protected:
  void setup_extras();

public:
  SimpleEntry extras;

  ElectronObject();
  ElectronObject(int id_);
  ElectronObject(int id_, TLorentzVector p4_);
  ElectronObject(const ElectronObject& other);
  ElectronObject& operator=(const ElectronObject& other);
  ~ElectronObject();

  float EinvOverPinv();

};

#endif
