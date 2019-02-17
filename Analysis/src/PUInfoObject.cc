#include <algorithm>
#include <utility>
#include "HelperFunctions.h"
#include "PUInfoObject.h"


PUInfoObject::PUInfoObject() :
  bunchCrossing(0),
  nPUVertices(0),
  nTrueVertices(0)
{}
PUInfoObject::PUInfoObject(const PUInfoObject& other) :
  bunchCrossing(other.bunchCrossing),
  nPUVertices(other.nPUVertices),
  nTrueVertices(other.nTrueVertices)
{}

PUInfoObject& PUInfoObject::operator=(const PUInfoObject& other){
  PUInfoObject tmp(other);
  swap(tmp);
  return *this;
}
void PUInfoObject::swap(PUInfoObject& other){
  std::swap(bunchCrossing, other.bunchCrossing);
  std::swap(nPUVertices, other.nPUVertices);
  std::swap(nTrueVertices, other.nTrueVertices);
}
