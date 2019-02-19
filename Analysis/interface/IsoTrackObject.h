#ifndef ISOTRACKOBJECT_H
#define ISOTRACKOBJECT_H

#include "ParticleObject.h"


class IsoTrackVariables{
public:
  bool isPFCand;
  bool hasLepOverlap;

  int charge;

  float pfIso_ch;
  float dz;

  IsoTrackVariables();
  IsoTrackVariables(IsoTrackVariables const& other);
  IsoTrackVariables& operator=(const IsoTrackVariables& other);

  void swap(IsoTrackVariables& other);

};

class IsoTrackObject : public ParticleObject{
public:
  IsoTrackVariables extras;

  IsoTrackObject();
  IsoTrackObject(int id_, CMSLorentzVector const& mom_);
  IsoTrackObject(const IsoTrackObject& other);
  IsoTrackObject& operator=(const IsoTrackObject& other);
  ~IsoTrackObject();

  void swap(IsoTrackObject& other);

};

#endif
