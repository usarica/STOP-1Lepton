#ifndef PARTICLEOBJECT_H
#define PARTICLEOBJECT_H

#include "TLorentzVector.h"
#include "CMSLorentzVector.h"
#include <DataFormats/Math/interface/deltaR.h>


class ParticleObject{
public:
  int id;
  long long selectionBits;
  CMSLorentzVector momentum;

  ParticleObject();
  ParticleObject(int id_);
  ParticleObject(int id_, CMSLorentzVector const& mom_);
  ParticleObject(const ParticleObject& other);
  virtual ~ParticleObject(){}

  virtual CMSLorentzVector getFinalMomentum() const{ return momentum; } // To be overloaded in daughter classes

  // Swap and assignment operators are not virtual; they bring more complication than necessary, so they are implemented in the derived classes.

  void resetSelectionBits(){ selectionBits=0; }
  void setSelectionBit(unsigned int ibit);
  bool testSelection(unsigned int ibit) const;

  float charge() const;
  float m() const{ return momentum.M(); }
  float x() const{ return momentum.X(); }
  float y() const{ return momentum.Y(); }
  float z() const{ return momentum.Z(); }
  float t() const{ return momentum.T(); }
  float energy() const{ return this->t(); }
  float p() const{ return momentum.P(); }
  float pt() const{ return momentum.Pt(); }
  float eta() const{ return momentum.Eta(); }
  float phi() const{ return momentum.Phi(); }
  float rapidity() const{ return momentum.Rapidity(); }
  float dot(const TLorentzVector& v) const{ return (momentum.T()*v.T()-(momentum.X()*v.X()+momentum.Y()*v.Y()+momentum.Z()*v.Z())); }
  float dot(const CMSLorentzVector& v) const{ return (momentum.T()*v.T()-(momentum.X()*v.X()+momentum.Y()*v.Y()+momentum.Z()*v.Z())); }
  float dot(const ParticleObject& part) const{ return dot(part.momentum); }
  float dot(const ParticleObject* part) const{ if (part!=0) return dot(*part); else return 0; }
  float deltaR(const TLorentzVector& v) const{ TLorentzVector tmp(momentum.X(), momentum.Y(), momentum.Z(), momentum.T()); return tmp.DeltaR(v); }
  float deltaR(const CMSLorentzVector& v) const{ return reco::deltaR(momentum, v); }
  float deltaR(const ParticleObject& part) const{ return deltaR(part.momentum); }
  float deltaR(const ParticleObject* part) const{ if (part!=0) return deltaR(*part); else return -1; }
  float deltaPhi(float phi_) const;
  TVector3 vect() const{ TLorentzVector tmp(momentum.X(), momentum.Y(), momentum.Z(), momentum.T()); return tmp.Vect(); }

};

#endif
