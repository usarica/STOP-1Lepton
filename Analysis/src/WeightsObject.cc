#include <cassert>
#include <algorithm>
#include <utility>
#include "WeightsObject.h"
#include "MELAStreamHelpers.hh"


WeightVariables::WeightVariables() :
  wgt_central(0),
  wgt_muF2(0),
  wgt_muF0p5(0),
  wgt_muR2(0),
  wgt_muR0p5(0),
  wgt_PDFVariationUp(0),
  wgt_PDFVariationDn(0),
  wgt_AsMZUp(0),
  wgt_AsMZDn(0),
  wgt_PSUp(0),
  wgt_PSDn(0)
{}
WeightVariables::WeightVariables(WeightVariables const& other) :
  wgt_central(other.wgt_central),
  wgt_muF2(other.wgt_muF2),
  wgt_muF0p5(other.wgt_muF0p5),
  wgt_muR2(other.wgt_muR2),
  wgt_muR0p5(other.wgt_muR0p5),
  wgt_PDFVariationUp(other.wgt_PDFVariationUp),
  wgt_PDFVariationDn(other.wgt_PDFVariationDn),
  wgt_AsMZUp(other.wgt_AsMZUp),
  wgt_AsMZDn(other.wgt_AsMZDn),
  wgt_PSUp(other.wgt_PSUp),
  wgt_PSDn(other.wgt_PSDn)

{}
void WeightVariables::swap(WeightVariables& other){
  std::swap(wgt_central, other.wgt_central);
  std::swap(wgt_muF2, other.wgt_muF2);
  std::swap(wgt_muF0p5, other.wgt_muF0p5);
  std::swap(wgt_muR2, other.wgt_muR2);
  std::swap(wgt_muR0p5, other.wgt_muR0p5);
  std::swap(wgt_PDFVariationUp, other.wgt_PDFVariationUp);
  std::swap(wgt_PDFVariationDn, other.wgt_PDFVariationDn);
  std::swap(wgt_AsMZUp, other.wgt_AsMZUp);
  std::swap(wgt_AsMZDn, other.wgt_AsMZDn);
  std::swap(wgt_PSUp, other.wgt_PSUp);
  std::swap(wgt_PSDn, other.wgt_PSDn);
}
WeightVariables& WeightVariables::operator=(const WeightVariables& other){
  WeightVariables tmp(other);
  swap(tmp);
  return *this;
}
TString WeightVariables::getWeightName(WeightVariables::WeightType type){
  using MELAStreamHelpers::MELAerr;
  using std::endl;

  TString prefix = "weight_";
  switch (type){
  case wCentral:
    return prefix + "Central";
  case wFacScaleUp:
    return prefix + "FacScaleUp";
  case wFacScaleDn:
    return prefix + "FacScaleDn";
  case wRenScaleUp:
    return prefix + "RenScaleUp";
  case wRenScaleDn:
    return prefix + "RenScaleDn";
  case wPDFUp:
    return prefix + "PDFUp";
  case wPDFDn:
    return prefix + "PDFDn";
  case wAsMZUp:
    return prefix + "AsMZUp";
  case wAsMZDn:
    return prefix + "AsMZDn";
  case wPSUp:
    return prefix + "PSUp";
  case wPSDn:
    return prefix + "PSDn";
  default:
    MELAerr << "WeightVariables::getWeightName: Weight type " << (int) type << " is unknown. Please modify this function." << endl;
    assert(0);
    return "";
  }
}
float WeightVariables::getWeight(WeightVariables::WeightType type) const{
  using MELAStreamHelpers::MELAerr;
  using std::endl;

  switch (type){
  case wCentral:
    return wgt_central;
  case wFacScaleUp:
    return wgt_muF2;
  case wFacScaleDn:
    return wgt_muF0p5;
  case wRenScaleUp:
    return wgt_muR2;
  case wRenScaleDn:
    return wgt_muR0p5;
  case wPDFUp:
    return wgt_PDFVariationUp;
  case wPDFDn:
    return wgt_PDFVariationDn;
  case wAsMZUp:
    return wgt_AsMZUp;
  case wAsMZDn:
    return wgt_AsMZDn;
  case wPSUp:
    return wgt_PSUp;
  case wPSDn:
    return wgt_PSDn;
  default:
    MELAerr << "WeightVariables::getWeight: Weight type " << (int) type << " is unknown. Please modify this function." << endl;
    assert(0);
    return 0;
  }
}



WeightsObject::WeightsObject() :
  extras()
{}
WeightsObject::WeightsObject(const WeightsObject& other) :
  extras(other.extras)
{}
void WeightsObject::swap(WeightsObject& other){
  extras.swap(other.extras);
}
WeightsObject& WeightsObject::operator=(const WeightsObject& other){
  WeightsObject tmp(other);
  swap(tmp);
  return *this;
}
WeightsObject::~WeightsObject(){}
