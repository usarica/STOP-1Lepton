#include <algorithm>
#include <utility>
#include "Samples.h"
#include "METObject.h"


METVariables::METVariables() :
  met_raw(0),
  phi_raw(0),
  met(0),
  phi(0),
  met_JEC(0),
  phi_JEC(0),
  met_JECup(0),
  phi_JECup(0),
  met_JECdn(0),
  phi_JECdn(0)
{}
METVariables::METVariables(METVariables const& other) :
  met_raw(other.met_raw),
  phi_raw(other.phi_raw),
  met(other.met),
  phi(other.phi),
  met_JEC(other.met_JEC),
  phi_JEC(other.phi_JEC),
  met_JECup(other.met_JECup),
  phi_JECup(other.phi_JECup),
  met_JECdn(other.met_JECdn),
  phi_JECdn(other.phi_JECdn)
{}
void METVariables::swap(METVariables& other){
  std::swap(met_raw, other.met_raw);
  std::swap(phi_raw, other.phi_raw);
  std::swap(met, other.met);
  std::swap(phi, other.phi);
  std::swap(met_JEC, other.met_JEC);
  std::swap(phi_JEC, other.phi_JEC);
  std::swap(met_JECup, other.met_JECup);
  std::swap(phi_JECup, other.phi_JECup);
  std::swap(met_JECdn, other.met_JECdn);
  std::swap(phi_JECdn, other.phi_JECdn);
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
