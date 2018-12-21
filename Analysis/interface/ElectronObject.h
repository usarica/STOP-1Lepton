#ifndef ELECTRONOBJECT_H
#define ELECTRONOBJECT_H

#include "TLorentzVector.h"
#include "CMSLorentzVector.h"


class ElectronVariables{
public:
  bool conv_vtx_flag; // Conversion veto

  int expectedMissingInnerHits; // expectedMissingInnerHits

  float energySC; // Supercluster energy, eSC
  float etaSC; // Supercluster eta
  float rho; // evt_fixgridfastjet_all_rho, actually a constant over the event
  float sigmaIEtaIEta_full5x5; // full5x5_sigmaIetaIeta
  float dEtaIn; // abs(dEtaIn)
  float dPhiIn; // abs(dPhiIn)
  float hOverE; // hOverE
  float ecalEnergy;
  float eOverPIn;
  float dxyPV; // abs(d0)
  float dzPV; // abs(dz)

  ElectronVariables();
  ElectronVariables(ElectronVariables const& other);
  ElectronVariables& operator=(const ElectronVariables& other);

  void swap(ElectronVariables& other);

};

class ElectronObject{
public:
  int id;
  CMSLorentzVector momentum;
  ElectronVariables extras;

  ElectronObject();
  ElectronObject(int id_);
  ElectronObject(int id_, CMSLorentzVector mom_);
  ElectronObject(const ElectronObject& other);
  ElectronObject& operator=(const ElectronObject& other);
  ~ElectronObject();

  void swap(ElectronObject& other);

  float EinvMinusPinv();

  float charge()const;
  float m()const{ return momentum.M(); }
  float x()const{ return momentum.X(); }
  float y()const{ return momentum.Y(); }
  float z()const{ return momentum.Z(); }
  float t()const{ return momentum.T(); }
  float p()const{ return momentum.P(); }
  float pt()const{ return momentum.Pt(); }
  float eta()const{ return momentum.Eta(); }
  float phi()const{ return momentum.Phi(); }
  float rapidity()const{ return momentum.Rapidity(); }
  float dot(const TLorentzVector& v)const{ return (momentum.T()*v.T()-(momentum.X()*v.X()+momentum.Y()*v.Y()+momentum.Z()*v.Z())); }
  float dot(const CMSLorentzVector& v)const{ return (momentum.T()*v.T()-(momentum.X()*v.X()+momentum.Y()*v.Y()+momentum.Z()*v.Z())); }
  float dot(const ElectronObject& part)const{ return dot(part.momentum); }
  float dot(const ElectronObject* part)const{ if (part!=0) return dot(*part); else return 0; }
  float deltaR(const TLorentzVector& v)const{ TLorentzVector tmp(momentum.X(), momentum.Y(), momentum.Z(), momentum.T()); return tmp.DeltaR(v); }
  float deltaR(const CMSLorentzVector& v)const{ TLorentzVector tmp(v.X(), v.Y(), v.Z(), v.T()); return deltaR(tmp); }
  float deltaR(const ElectronObject& part)const{ return deltaR(part.momentum); }
  float deltaR(const ElectronObject* part)const{ if (part!=0) return deltaR(*part); else return -1; }
  TVector3 vect()const{ TLorentzVector tmp(momentum.X(), momentum.Y(), momentum.Z(), momentum.T()); return tmp.Vect(); }

};

#endif
