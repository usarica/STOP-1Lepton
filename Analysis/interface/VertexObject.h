#ifndef VERTEXOBJECT_H
#define VERTEXOBJECT_H

#include "CMSLorentzVector.h"
#include "TVector3.h"


class VertexObject{
public:
  bool isValid;
  bool isFake;
  float ndof;
  long long selectionBits;
  TVector3 position;


  VertexObject();
  VertexObject(CMSLorentzVector const& position_);
  VertexObject(TVector3 const& position_);
  VertexObject(VertexObject const& other);
  ~VertexObject(){}

  VertexObject& operator=(const VertexObject& other);
  void swap(VertexObject& other);

  void setSelectionBit(unsigned int ibit);
  bool testSelection(unsigned int ibit) const;

  float x() const{ return position.X(); }
  float y() const{ return position.Y(); }
  float z() const{ return position.Z(); }
  float mag() const{ return position.Mag(); }
  float rho() const{ return position.Perp(); }
  float eta() const{ return position.Eta(); }
  float phi() const{ return position.Phi(); }
  float deltaPhi(float phi_) const;
};

#endif
