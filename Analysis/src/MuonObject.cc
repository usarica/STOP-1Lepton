#include <algorithm>
#include <utility>
#include "Samples.h"
#include "MuonObject.h"


MuonVariables::MuonVariables() :
  isPFMuon(false),

  type(0),
  validHits(0),
  lostHits(0),
  expectedMissingInnerHits(0),
  expectedMissingOuterHits(0),
  GlobalFit_Ndof(0),

  POGSelectorBit(0),

  rho(0),
  GlobalFit_Chisq(0),
  LocalPos_Chisq(0),
  TrkKink(0),
  SegComp(0),
  dxyPV(0),
  dzPV(0),
  miniIso_ch(0),
  miniIso_nh(0),
  miniIso_em(0)
{}
MuonVariables::MuonVariables(MuonVariables const& other) :
  isPFMuon(other.isPFMuon),

  type(other.type),
  validHits(other.validHits),
  lostHits(other.lostHits),
  expectedMissingInnerHits(other.expectedMissingInnerHits),
  expectedMissingOuterHits(other.expectedMissingOuterHits),
  GlobalFit_Ndof(other.GlobalFit_Ndof),

  POGSelectorBit(other.POGSelectorBit),

  rho(other.rho),
  GlobalFit_Chisq(other.GlobalFit_Chisq),
  LocalPos_Chisq(other.LocalPos_Chisq),
  TrkKink(other.TrkKink),
  SegComp(other.SegComp),
  dxyPV(other.dxyPV),
  dzPV(other.dzPV),
  miniIso_ch(other.miniIso_ch),
  miniIso_nh(other.miniIso_nh),
  miniIso_em(other.miniIso_em)
{}
void MuonVariables::swap(MuonVariables& other){
  std::swap(isPFMuon, other.isPFMuon);

  std::swap(type, other.type);
  std::swap(validHits, other.validHits);
  std::swap(lostHits, other.lostHits);
  std::swap(expectedMissingInnerHits, other.expectedMissingInnerHits);
  std::swap(expectedMissingOuterHits, other.expectedMissingOuterHits);
  std::swap(GlobalFit_Ndof, other.GlobalFit_Ndof);

  std::swap(POGSelectorBit, other.POGSelectorBit);

  std::swap(rho, other.rho);
  std::swap(GlobalFit_Chisq, other.GlobalFit_Chisq);
  std::swap(LocalPos_Chisq, other.LocalPos_Chisq);
  std::swap(TrkKink, other.TrkKink);
  std::swap(SegComp, other.SegComp);
  std::swap(dxyPV, other.dxyPV);
  std::swap(dzPV, other.dzPV);
  std::swap(miniIso_ch, other.miniIso_ch);
  std::swap(miniIso_nh, other.miniIso_nh);
  std::swap(miniIso_em, other.miniIso_em);
}
MuonVariables& MuonVariables::operator=(const MuonVariables& other){
  MuonVariables tmp(other);
  swap(tmp);
  return *this;
}


MuonObject::MuonObject() :
  ParticleObject(),
  extras()
{}
MuonObject::MuonObject(int id_) :
  ParticleObject(id_),
  extras()
{}
MuonObject::MuonObject(int id_, CMSLorentzVector momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
MuonObject::MuonObject(const MuonObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void MuonObject::swap(MuonObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
MuonObject& MuonObject::operator=(const MuonObject& other){
  MuonObject tmp(other);
  swap(tmp);
  return *this;
}
MuonObject::~MuonObject(){}
