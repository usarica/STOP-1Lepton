#ifndef TFTOPOBJECT_H
#define TFTOPOBJECT_H

#include <vector>
#include "ParticleObject.h"


class TFTopVariables{
public:
  int nSubjets;
  float disc;
  std::vector<TLorentzVector> subjets_momentum;

  TFTopVariables();
  TFTopVariables(TFTopVariables const& other);
  TFTopVariables& operator=(const TFTopVariables& other);

  void swap(TFTopVariables& other);

};

class TFTopObject : public ParticleObject{
public:
  TFTopVariables extras;

  TFTopObject();
  TFTopObject(int id_);
  TFTopObject(int id_, CMSLorentzVector mom_);
  TFTopObject(const TFTopObject& other);
  TFTopObject& operator=(const TFTopObject& other);
  ~TFTopObject();

  void swap(TFTopObject& other);

};

#endif
