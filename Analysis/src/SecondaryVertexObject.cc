#include <algorithm>
#include <utility>
#include "HelperFunctions.h"
#include "SecondaryVertexObject.h"


SecondaryVertexObject::SecondaryVertexObject() :
  nTracks(0),
  IP2D(-1),
  SIP2D(-1),
  IP3D(-1),
  SIP3D(-1),
  anglePV(0),
  selectionBits(0),
  position(0, 0, 0),
  momentum(0, 0, 0, 0)
{}
SecondaryVertexObject::SecondaryVertexObject(const SecondaryVertexObject& other) :
  nTracks(other.nTracks),
  IP2D(other.IP2D),
  SIP2D(other.SIP2D),
  IP3D(other.IP3D),
  SIP3D(other.SIP3D),
  anglePV(other.anglePV),
  selectionBits(other.selectionBits),
  position(other.position),
  momentum(other.momentum)
{}

void SecondaryVertexObject::setSelectionBit(unsigned int ibit){ HelperFunctions::set_bit(this->selectionBits, ibit); }
bool SecondaryVertexObject::testSelection(unsigned int ibit) const{ return HelperFunctions::test_bit(this->selectionBits, ibit); }

SecondaryVertexObject& SecondaryVertexObject::operator=(const SecondaryVertexObject& other){
  SecondaryVertexObject tmp(other);
  swap(tmp);
  return *this;
}
void SecondaryVertexObject::swap(SecondaryVertexObject& other){
  std::swap(nTracks, other.nTracks);
  std::swap(IP2D, other.IP2D);
  std::swap(SIP2D, other.SIP2D);
  std::swap(IP3D, other.IP3D);
  std::swap(SIP3D, other.SIP3D);
  std::swap(anglePV, other.anglePV);
  std::swap(selectionBits, other.selectionBits);
  std::swap(position, other.position);
  std::swap(momentum, other.momentum);
}
