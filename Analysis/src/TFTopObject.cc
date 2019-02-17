#include <algorithm>
#include <utility>
#include "Samples.h"
#include "TFTopObject.h"


TFTopVariables::TFTopVariables() :
  nSubjets(0),
  disc(-1)
{}
TFTopVariables::TFTopVariables(TFTopVariables const& other) :
  nSubjets(other.nSubjets),
  disc(other.disc),
  subjets_momentum(other.subjets_momentum)
{}
void TFTopVariables::swap(TFTopVariables& other){
  std::swap(nSubjets, other.nSubjets);
  std::swap(disc, other.disc);
  std::swap(subjets_momentum, other.subjets_momentum);
}
TFTopVariables& TFTopVariables::operator=(const TFTopVariables& other){
  TFTopVariables tmp(other);
  swap(tmp);
  return *this;
}


TFTopObject::TFTopObject() :
  ParticleObject(),
  extras()
{}
TFTopObject::TFTopObject(int id_) :
  ParticleObject(id_),
  extras()
{}
TFTopObject::TFTopObject(int id_, CMSLorentzVector const& momentum_) :
  ParticleObject(id_, momentum_),
  extras()
{}
TFTopObject::TFTopObject(const TFTopObject& other) :
  ParticleObject(other),
  extras(other.extras)
{}
void TFTopObject::swap(TFTopObject& other){
  std::swap(id, other.id);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  extras.swap(other.extras);
}
TFTopObject& TFTopObject::operator=(const TFTopObject& other){
  TFTopObject tmp(other);
  swap(tmp);
  return *this;
}
TFTopObject::~TFTopObject(){}
