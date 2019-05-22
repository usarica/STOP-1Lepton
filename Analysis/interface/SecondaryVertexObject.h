#ifndef SECONDARYVERTEXOBJECT_H
#define SECONDARYVERTEXOBJECT_H

#include "CMSLorentzVector.h"
#include "TLorentzVector.h"


class SecondaryVertexObject{
public:
  unsigned int nTracks;
  float IP2D;
  float SIP2D;
  float IP3D;
  float SIP3D;
  float anglePV;
  long long selectionBits;
  TVector3 position;
  CMSLorentzVector momentum;

  SecondaryVertexObject();
  SecondaryVertexObject(SecondaryVertexObject const& other);
  ~SecondaryVertexObject(){}

  SecondaryVertexObject& operator=(const SecondaryVertexObject& other);
  void swap(SecondaryVertexObject& other);

  void setSelectionBit(unsigned int ibit);
  bool testSelection(unsigned int ibit) const;

  float x() const{ return position.X(); }
  float y() const{ return position.Y(); }
  float z() const{ return position.Z(); }
  float mag() const{ return position.Mag(); }
  float rho() const{ return position.Perp(); }
  float posEta() const{ return position.Eta(); }
  float posPhi() const{ return position.Phi(); }
  float px() const{ return momentum.X(); }
  float py() const{ return momentum.Y(); }
  float pz() const{ return momentum.Z(); }
  float pt() const{ return momentum.Pt(); }
  float momEta() const{ return momentum.Eta(); }
  float momPhi() const{ return momentum.Phi(); }
  float energy() const{ return momentum.T(); }
  float m() const{ return momentum.M(); }
};

#endif
