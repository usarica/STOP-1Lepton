#ifndef MUONOBJECT_H
#define MUONOBJECT_H

#include "ParticleObject.h"


class MuonVariables{
public:
  bool isPFMuon; // This is an int in the framework, just convert it to a bool (extra bits are not needed at all).

  int type;
  int validHits;
  int lostHits;
  int expectedMissingInnerHits;
  int expectedMissingOuterHits;
  int GlobalFit_Ndof;

  long long POGSelectorBit; // This is an unsigned int in the framework. (Somebody should tell signed vs unsigned does not matter in bit tests...)

  float rho;
  float GlobalFit_Chisq;
  float LocalPos_Chisq;
  float TrkKink;
  float SegComp;
  float dxyPV;
  float dzPV;
  float miniIso_ch;
  float miniIso_nh;
  float miniIso_em;

  MuonVariables();
  MuonVariables(MuonVariables const& other);
  MuonVariables& operator=(const MuonVariables& other);

  void swap(MuonVariables& other);

};

class MuonObject : public ParticleObject{
public:
  MuonVariables extras;

  MuonObject();
  MuonObject(int id_);
  MuonObject(int id_, CMSLorentzVector mom_);
  MuonObject(const MuonObject& other);
  MuonObject& operator=(const MuonObject& other);
  ~MuonObject();

  void swap(MuonObject& other);

};

#endif
