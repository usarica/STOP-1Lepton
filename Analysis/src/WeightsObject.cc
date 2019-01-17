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
  wgt_PSDn(0),

  wgt_central_default(0),
  wgt_PDFVariationUp_default(0),
  wgt_PDFVariationDn_default(0),
  wgt_ISRUp(0),
  wgt_ISRDn(0),
  wgt_FSRUp(0),
  wgt_FSRDn(0)
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
  wgt_PSDn(other.wgt_PSDn),

  wgt_central_default(other.wgt_central_default),
  wgt_PDFVariationUp_default(other.wgt_PDFVariationUp_default),
  wgt_PDFVariationDn_default(other.wgt_PDFVariationDn_default),
  wgt_ISRUp(other.wgt_ISRUp),
  wgt_ISRDn(other.wgt_ISRDn),
  wgt_FSRUp(other.wgt_FSRUp),
  wgt_FSRDn(other.wgt_FSRDn)
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


  std::swap(wgt_central_default, other.wgt_central_default);
  std::swap(wgt_PDFVariationUp_default, other.wgt_PDFVariationUp_default);
  std::swap(wgt_PDFVariationDn_default, other.wgt_PDFVariationDn_default);
  std::swap(wgt_ISRUp, other.wgt_ISRUp);
  std::swap(wgt_ISRDn, other.wgt_ISRDn);
  std::swap(wgt_FSRUp, other.wgt_FSRUp);
  std::swap(wgt_FSRDn, other.wgt_FSRDn);
}
WeightVariables& WeightVariables::operator=(const WeightVariables& other){
  WeightVariables tmp(other);
  swap(tmp);
  return *this;
}

TString WeightVariables::getWeightName(WeightVariables::WeightType type){
  using MELAStreamHelpers::MELAerr;
  using std::endl;

  TString prefix = "genweight";
  switch (type){
  case wCentral:
    return prefix;
  case wFacScaleUp:
    return prefix + "_FacScaleUp";
  case wFacScaleDn:
    return prefix + "_FacScaleDn";
  case wRenScaleUp:
    return prefix + "_RenScaleUp";
  case wRenScaleDn:
    return prefix + "_RenScaleDn";
  case wPDFUp:
    return prefix + "_PDFUp";
  case wPDFDn:
    return prefix + "_PDFDn";
  case wAsMZUp:
    return prefix + "_AsMZUp";
  case wAsMZDn:
    return prefix + "_AsMZDn";
  case wPSUp:
    return prefix + "_PSUp";
  case wPSDn:
    return prefix + "_PSDn";

  case wCentral_Default:
    return prefix + "_Default";
  case wPDFUp_Default:
    return prefix + "_Default_PDFUp";
  case wPDFDn_Default:
    return prefix + "_Default_PDFDn";
  case wISRUp:
    return prefix + "_ISRUp";
  case wISRDn:
    return prefix + "_ISRDn";
  case wFSRUp:
    return prefix + "_FSRUp";
  case wFSRDn:
    return prefix + "_FSRDn";

  default:
    MELAerr << "WeightVariables::getWeightName: Weight type " << (int) type << " is unknown. Please modify this function." << endl;
    assert(0);
    return "";
  }
}
TString WeightVariables::getWeightLabel(WeightVariables::WeightType type, bool use2016){
  using MELAStreamHelpers::MELAerr;
  using std::endl;

  switch (type){
  case wCentral:
    return (use2016 ? "Nominal (2016)" : "Nominal (2017+)");
  case wFacScaleUp:
    return "#mu_{F} up";
  case wFacScaleDn:
    return "#mu_{F} dn";
  case wRenScaleUp:
    return "#mu_{R} up";
  case wRenScaleDn:
    return "#mu_{R} dn";
  case wPDFUp:
    return (use2016 ? "PDF (2016) up" : "PDF up");
  case wPDFDn:
    return (use2016 ? "PDF (2016) dn" : "PDF dn");
  case wAsMZUp:
    return "#alpha_{s}(m_{Z}) up";
  case wAsMZDn:
    return "#alpha_{s}(m_{Z}) dn";
  case wPSUp:
    return "ISR*FSR up";
  case wPSDn:
    return "ISR*FSR dn";

  case wCentral_Default:
    return "Nominal (def. PDF)";
  case wPDFUp_Default:
    return "Def. PDF up";
  case wPDFDn_Default:
    return "Def. PDF dn";
  case wISRUp:
    return "ISR up";
  case wISRDn:
    return "ISR dn";
  case wFSRUp:
    return "FSR up";
  case wFSRDn:
    return "FSR dn";

  default:
    MELAerr << "WeightVariables::getWeightLabel: Weight type " << (int) type << " is unknown. Please modify this function." << endl;
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

  case wCentral_Default:
    return wgt_central_default;
  case wPDFUp_Default:
    return wgt_PDFVariationUp_default;
  case wPDFDn_Default:
    return wgt_PDFVariationDn_default;
  case wISRUp:
    return wgt_ISRUp;
  case wISRDn:
    return wgt_ISRDn;
  case wFSRUp:
    return wgt_FSRUp;
  case wFSRDn:
    return wgt_FSRDn;

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
