#include <algorithm>
#include <utility>
#include "HelperFunctions.h"
#include "VertexObject.h"


VertexObject::VertexObject() :
  isValid(false),
  isFake(false),
  ndof(-1),
  selectionBits(0),
  position(0, 0, 0)
{}
VertexObject::VertexObject(CMSLorentzVector const& position_) :
  isValid(false),
  isFake(false),
  ndof(-1),
  selectionBits(0),
  position(position_.x(), position_.y(), position_.z())
{}
VertexObject::VertexObject(TVector3 const& position_) :
  isValid(false),
  isFake(false),
  ndof(-1),
  selectionBits(0),
  position(position_)
{}
VertexObject::VertexObject(const VertexObject& other) :
  isValid(other.isValid),
  isFake(other.isFake),
  ndof(other.ndof),
  selectionBits(other.selectionBits),
  position(other.position)
{}

void VertexObject::setSelectionBit(unsigned int ibit){ HelperFunctions::set_bit(this->selectionBits, ibit); }
bool VertexObject::testSelection(unsigned int ibit) const{ return HelperFunctions::test_bit(this->selectionBits, ibit); }

VertexObject& VertexObject::operator=(const VertexObject& other){
  VertexObject tmp(other);
  swap(tmp);
  return *this;
}
void VertexObject::swap(VertexObject& other){
  std::swap(isValid, other.isValid);
  std::swap(isFake, other.isFake);
  std::swap(ndof, other.ndof);
  std::swap(selectionBits, other.selectionBits);
  std::swap(position, other.position);
}


float VertexObject::deltaPhi(float phi_) const{
  float dPhi = phi_-phi();
  if (dPhi>TMath::Pi()) dPhi = 2.*TMath::Pi() - dPhi;
  else if (dPhi<-TMath::Pi()) dPhi = 2.*TMath::Pi() + dPhi;
  return dPhi;
}
