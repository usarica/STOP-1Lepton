#include <algorithm>
#include <utility>
#include "Samples.h"
#include "METObject.h"


METVariables::METVariables() :
  met(0),
  phi(0)
{}
METVariables::METVariables(METVariables const& other) :
  met(other.met),
  phi(other.phi)
{}
void METVariables::swap(METVariables& other){
  std::swap(met, other.met);
  std::swap(phi, other.phi);
}
METVariables& METVariables::operator=(const METVariables& other){
  METVariables tmp(other);
  swap(tmp);
  return *this;
}


METObject::METObject() :
  extras()
{}
METObject::METObject(const METObject& other) :
  extras(other.extras)
{}
void METObject::swap(METObject& other){
  extras.swap(other.extras);
}
METObject& METObject::operator=(const METObject& other){
  METObject tmp(other);
  swap(tmp);
  return *this;
}
METObject::~METObject(){}
