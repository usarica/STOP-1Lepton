#include <FWCore/MessageLogger/interface/MessageLogger.h>
#include <FWCore/Utilities/interface/EDMException.h>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include "TLorentzVector.h"
#include "TString.h"
#include "LHEWeightHandler.h"


using namespace std;


LHEWeightHandler::LHEWeightHandler(int year_, LHEWeightHandler::PDFChoice pdfChoice_, LHEWeightHandler::QCDOrderChoice orderChoice_) :
  exceptionalCases(),
  year(year_),
  pdfChoice(pdfChoice_),
  orderChoice(orderChoice_)
{
  if (year==2016 && pdfChoice==tryNNPDF31){
    edm::LogWarning warning("InvalidPDFChoice");
    warning << "MC in 2016 only includes NNPDF 3.0. Defaulting the pdfChoice flag to keepDefaultPDF";
    pdfChoice=keepDefaultPDF;
  }
  if (pdfChoice==tryNNPDF30 && orderChoice==tryNNLO && year==2016){
    edm::LogWarning warning("InvalidPDFChoice");
    warning << "No NNPDF 3.0 NNLO weights are present in " << year << " MC. Will default orderChoice to tryNLO.";
    orderChoice=tryNLO;
  }
  clear();
}
LHEWeightHandler::~LHEWeightHandler(){ clear(); }

void LHEWeightHandler::clear(){
  LHEWeight.clear();
  LHEWeight_PDFVariationUpDn.clear();
  LHEWeight_AsMZUpDn.clear();
  LHEOriginalWeight=1;
  defaultMemberZeroWeight=1;
  defaultWeightScale=1;
}

float LHEWeightHandler::safeDivide(float numerator, float denominator){ return (!(std::isfinite(numerator) && std::isfinite(denominator)) || denominator==0.f ? 0.f : numerator/denominator); }

float const& LHEWeightHandler::getLHEOriginalWeight() const{ return this->LHEOriginalWeight; }
float const& LHEWeightHandler::getMemberZeroWeight() const{ return this->defaultMemberZeroWeight; }
float const& LHEWeightHandler::getWeightRescale() const{
  return defaultWeightScale; // Note this is already divided by defaultMemberZeroWeight
}
float LHEWeightHandler::getLHEWeight(unsigned int whichWeight, float defaultValue) const{
  if (whichWeight<LHEWeight.size()) return LHEWeight.at(whichWeight);
  else return defaultValue;
}
float LHEWeightHandler::getLHEWeight_PDFVariationUpDn(int whichUpDn, float defaultValue) const{
  if (whichUpDn<0 && LHEWeight_PDFVariationUpDn.size()>1) return LHEWeight_PDFVariationUpDn.at(1);
  else if (whichUpDn>0 && LHEWeight_PDFVariationUpDn.size()>0) return LHEWeight_PDFVariationUpDn.at(0);
  else return defaultValue;
}
float LHEWeightHandler::getLHEWeigh_AsMZUpDn(int whichUpDn, float defaultValue) const{
  if (whichUpDn<0 && LHEWeight_AsMZUpDn.size()>1) return LHEWeight_AsMZUpDn.at(1);
  else if (whichUpDn>0 && LHEWeight_AsMZUpDn.size()>0) return LHEWeight_AsMZUpDn.at(0);
  else return defaultValue;
}


void LHEWeightHandler::extract(float wgt_central, std::vector<float> const& weights, std::vector<std::string> const& weightIds){
  // PDF information
  bool const& specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1 = exceptionalCases.specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1;
  bool const& specialPDF_NNPDF31_NNLO_as_0118_nf_4 = exceptionalCases.specialPDF_NNPDF31_NNLO_as_0118_nf_4;
  bool const& specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 = exceptionalCases.specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1;

  // LHE weights (Maximum 9 or size of weights array for QCD variations, 2 for PDF variations and 2 for alphas(mZ) variations)
  LHEOriginalWeight = wgt_central;
  vector<float> LHEPDFAlphaSMZWgt;
  vector<float> LHEPDFVariationWgt;
  AlternateWeightsType weightstype = unknown;
  bool found_defaultMemberZeroWeight=false;
  bool found_defaultWeightScale=false;

  PDFVariationMode pdfVarMode = useNone;
  std::vector<std::string>::const_iterator it_weightid = weightIds.cbegin();
  for (const float& weight : weights){
    std::string const& weight_id = *it_weightid;
    if (weight_id=="") continue;
    int wgtid;
    // We don't use non-numerical indices, but they exist in some 2016 MC samples
    try{ wgtid = stoi(weight_id.c_str()); }
    catch (std::invalid_argument& e){ continue; }

    float wgtval = weight;
    //cout << "PDF id = " << PDFid.at(0) << " " << wgtid << " -> " << wgtval << endl;
    if (year == 2016){
      // Do not assign defaultMemberZeroWeight or defaultWeightScale. This will be handled correctly later on via LHEWeight[0].
      if (weightstype == unknown && 1001 <= wgtid && wgtid <= 1009){ // This is the NNPDF 3.0 NLO central value and muR/muF variations
        LHEWeight.push_back(wgtval);
      }
      else if (weightstype == unknown && (wgtid == 2000 || wgtid == 2001)){ // This is the NNPDF 3.0 NLO variations, but do nothing if wgtid == 2000 since it is the central value
        weightstype = powheg;
        if (wgtid == 2001){
          LHEPDFVariationWgt.push_back(wgtval);
          if (pdfVarMode == useNone) pdfVarMode = useMC;
        }
      }
      else if (weightstype == powheg && 2001 <= wgtid && wgtid <= 2102){ // This is the NNPDF 3.0 NLO variations
        LHEPDFVariationWgt.push_back(wgtval);
        if (pdfVarMode == useNone) pdfVarMode = useMC;
      }

      else if (weightstype == unknown && wgtid == 1){ // This is the NNPDF 3.0 LO a_s=0.130 central value
        weightstype = madgraph_0offset;
        LHEWeight.push_back(wgtval);
      }
      else if (weightstype == madgraph_0offset && 2 <= wgtid && wgtid <= 9){ // This is the NNPDF 3.0 LO a_s=0.130 muR/muF variations
        LHEWeight.push_back(wgtval);
      }
      else if (weightstype == madgraph_0offset && 10 <= wgtid && wgtid <= 110){ // This is the NNPDF 3.0 LO a_s=0.130 PDF variations
        // Notice that there could be 101 variations. This is handled later.
        LHEPDFVariationWgt.push_back(wgtval);
        if (pdfVarMode == useNone) pdfVarMode = useMC;
      }

      else if (weightstype == unknown && wgtid == 1010){ // Begins NNPDF 3.0 NLO a_s=0.118. This is member 0, so do nothing.
        weightstype = madgraph_1000offset;
      }
      else if (weightstype == madgraph_1000offset && 1011 <= wgtid && wgtid <= 1112){ // These are the NNPDF 3.0 NLO PDF and a_s variations
        LHEPDFVariationWgt.push_back(wgtval);
        if (pdfVarMode == useNone) pdfVarMode = useMC;
      }
    }
    else if (year == 2017 || year == 2018){
      //Madgraph 0 offset
      if (weightstype == unknown && wgtid == 1){
        weightstype = madgraph_0offset;
        LHEWeight.push_back(wgtval);
      }
      else if (weightstype == madgraph_0offset && 2 <= wgtid && wgtid <= 9) LHEWeight.push_back(wgtval);
      else if (weightstype == madgraph_0offset && wgtid == 10){ // Begins NNPDF 3.1 NNLO a_s=0.118
        defaultMemberZeroWeight = wgtval; found_defaultMemberZeroWeight=true;
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (weightstype == madgraph_0offset && 11 <= wgtid && wgtid <= 112){ // These are the NNPDF 3.1 NNLO PDF (<=110) and a_s (111, 112) variations (Hessian-type)
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){
          LHEPDFVariationWgt.push_back(wgtval);
          if (pdfVarMode == useNone) pdfVarMode = useHessian;
        }
      }
      else if (weightstype == madgraph_0offset && (wgtid == 117 || wgtid == 118)){/*do nothing, these are the NNLO a_s=0.117/0.119 variations from a slightly different pdf set*/ }
      else if (weightstype == madgraph_0offset && 113 <= wgtid && wgtid <= 120){/*do nothing, these are other NNLO a_s value variations*/ }
      else if (weightstype == madgraph_0offset && 121 <= wgtid && wgtid <= 223){
        /*These are the NNPDF 3.1 NLO PDF and a_s variations (Hessian-type)*/
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && orderChoice==tryNLO){
          if (wgtid == 121){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
          else{
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useHessian;
          }
        }
      }
      else if (weightstype == madgraph_0offset && 973 <= wgtid && wgtid <= 1075){ // These are the NNPDF 3.0 NLO PDF and MC variations
        if (pdfChoice==tryNNPDF30 && orderChoice==tryNLO){
          if (wgtid == 973){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
          else{
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useMC;
          }
        }
      }
      else if (weightstype == madgraph_0offset && wgtid == 1076){ // This is NNPDF 3.0 NNLO with a_s=0.118, no variations
        if (pdfChoice==tryNNPDF30 && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (weightstype == madgraph_0offset && wgtid == 1078){ // This is NNPDF 3.1 LO with a_s=0.130, no variations
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && orderChoice==tryLO){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (weightstype == madgraph_0offset && wgtid == 1080){ // This is NNPDF 3.0 LO with a_s=0.130, no variations
        if (pdfChoice==tryNNPDF30 && orderChoice==tryLO){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (weightstype == madgraph_0offset && 224 <= wgtid && wgtid <= 1080){/*do nothing, these are other various weights*/ }

      //QCD variations for all the other weight types
      else if (specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && weightstype == unknown && 1001 <= wgtid && wgtid <= 1045){
        if (wgtid % 5 == 1) LHEWeight.push_back(wgtval);
      }
      else if (weightstype == unknown && 1001 <= wgtid && wgtid <= 1009){
        LHEWeight.push_back(wgtval);
      }

      //Madgraph 1000 offset
      else if (
        (!specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && wgtid == 1010 && weightstype == unknown)
        ||
        (specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && wgtid == 1046 && weightstype == unknown)
        ){ // Begins NNPDF 3.1 NNLO a_s=0.118
        weightstype = madgraph_1000offset;
        defaultMemberZeroWeight = wgtval; found_defaultMemberZeroWeight=true;
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      // specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 with madgraph_1000offset
      else if (specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && weightstype == madgraph_1000offset && 1047 <= wgtid && wgtid <= 1148){ // These are the NNPDF 3.1 NNLO PDF and a_s variations (Hessian-type)
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){
          LHEPDFVariationWgt.push_back(wgtval);
          if (pdfVarMode == useNone) pdfVarMode = useHessian;
        }
      }
      else if (specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && weightstype == madgraph_1000offset && (wgtid == 1153 || wgtid == 1154)){ /*do nothing, these are the NNLO a_s=0.117/0.119 variations from a slightly different pdf set*/ }
      else if (specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && weightstype == madgraph_1000offset && 1149 <= wgtid && wgtid <= 1156){ /*do nothing, these are other NNLO a_s value variations*/ }
      else if (specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && weightstype == madgraph_1000offset && 1157 <= wgtid && wgtid <= 1259){ // These are the NNPDF 3.1 NLO PDF and a_s variations (Hessian-type)
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && orderChoice==tryNLO){
          if (wgtid == 1157){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
          else{
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useHessian;
          }
        }
      }
      else if (specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && weightstype == madgraph_1000offset && 1260 <= wgtid && wgtid <= 2008){ /*do nothing, these are other various weights*/ }
      else if (specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && weightstype == madgraph_1000offset && 2009 <= wgtid && wgtid <= 2111){ // These are the NNPDF 3.0 NLO PDF and MC variations
        if (pdfChoice==tryNNPDF30 && orderChoice==tryNLO){
          if (wgtid == 2009){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
          else{
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useMC;
          }
        }
      }
      else if (specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && weightstype == madgraph_1000offset && wgtid == 2112){ // This is NNPDF 3.0 NNLO with a_s=0.118, no variations
        if (pdfChoice==tryNNPDF30 && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && weightstype == madgraph_1000offset && wgtid == 2114){ // This is NNPDF 3.1 LO with a_s=0.130, no variations
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && orderChoice==tryLO){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && weightstype == madgraph_1000offset && wgtid == 2116){ // This is NNPDF 3.0 LO with a_s=0.130, no variations
        if (pdfChoice==tryNNPDF30 && orderChoice==tryLO){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1 && weightstype == madgraph_1000offset && 2112 <= wgtid && wgtid <= 2116){/*do nothing, these are other various weights*/ }
      // specialPDF_NNPDF31_NNLO_as_0118_nf_4 with madgraph_1000offset
      else if (specialPDF_NNPDF31_NNLO_as_0118_nf_4 && weightstype == madgraph_1000offset && 1011 <= wgtid && wgtid <= 1110){ // These are the NNPDF 3.1 NNLO PDF variations (MC-type)
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){
          LHEPDFVariationWgt.push_back(wgtval);
          if (pdfVarMode == useNone) pdfVarMode = useMC;
        }
      }
      else if (specialPDF_NNPDF31_NNLO_as_0118_nf_4 && weightstype == madgraph_1000offset && 1575 <= wgtid && wgtid <= 1675){ // These are the NNPDF 3.1 NLO PDF variations (MC-type)
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && orderChoice==tryNLO){
          if (wgtid == 1575){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
          else{
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useMC;
          }
        }
      }
      else if (specialPDF_NNPDF31_NNLO_as_0118_nf_4 && weightstype == madgraph_1000offset && 1779 <= wgtid && wgtid <= 1881){ // These are the NNPDF 3.0 NLO PDF and MC variations
        if (pdfChoice==tryNNPDF30 && orderChoice==tryNLO){
          if (wgtid == 1779){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
          else{
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useMC;
          }
        }
      }
      else if (specialPDF_NNPDF31_NNLO_as_0118_nf_4 && weightstype == madgraph_1000offset && wgtid == 1778){ // This is NNPDF 3.0 LO with a_s=0.130, no variations
        if (pdfChoice==tryNNPDF30 && orderChoice==tryLO){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (specialPDF_NNPDF31_NNLO_as_0118_nf_4 && weightstype == madgraph_1000offset){ /*Ignore all other weights*/ }
      // Regular madgraph_1000offset
      else if (weightstype == madgraph_1000offset && 1011 <= wgtid && wgtid <= 1112){ // These are the NNPDF 3.1 NNLO PDF (<=1110) and a_s (1111, 1112) variations (Hessian-type)
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){
          LHEPDFVariationWgt.push_back(wgtval);
          if (pdfVarMode == useNone) pdfVarMode = useHessian;
        }
      }
      else if (weightstype == madgraph_1000offset && (wgtid == 1117 || wgtid == 1118)){ /*do nothing, these are the NNLO a_s=0.117/0.119 variations from a slightly different pdf set*/ }
      else if (weightstype == madgraph_1000offset && 1113 <= wgtid && wgtid <= 1120){ /*do nothing, these are other NNLO a_s value variations*/ }
      else if (weightstype == madgraph_1000offset && 1121 <= wgtid && wgtid <= 1223){ // These are the NNPDF 3.1 NLO PDF and a_s variations (Hessian-type)
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && orderChoice==tryNLO){
          if (wgtid == 1121){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
          else{
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useHessian;
          }
        }
      }
      else if (weightstype == madgraph_1000offset && 1224 <= wgtid && wgtid <= 1972){ /*do nothing, these are other various weights*/ }
      else if (weightstype == madgraph_1000offset && 1973 <= wgtid && wgtid <= 2075){ // These are the NNPDF 3.0 NLO PDF and MC variations
        if (pdfChoice==tryNNPDF30 && orderChoice==tryNLO){
          if (wgtid == 1973){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
          else{
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useMC;
          }
        }
      }
      else if (weightstype == madgraph_1000offset && wgtid == 2076){ // This is NNPDF 3.0 NNLO with a_s=0.118, no variations
        if (pdfChoice==tryNNPDF30 && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (weightstype == madgraph_1000offset && wgtid == 2078){ // This is NNPDF 3.1 LO with a_s=0.130, no variations
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && orderChoice==tryLO){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (weightstype == madgraph_1000offset && wgtid == 2080){ // This is NNPDF 3.0 LO with a_s=0.130, no variations
        if (pdfChoice==tryNNPDF30 && orderChoice==tryLO){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (weightstype == madgraph_1000offset && 2076 <= wgtid && wgtid <= 2080){/*do nothing, these are other various weights*/ }

      //powheg
      else if (weightstype == unknown && (wgtid == 2000 || wgtid == 2001)){ // This is the NNPDF 3.1 NNLO central value, but do nothing if wgtid == 2000
        weightstype = powheg;
        if (specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1){ // How fun is that?! This special case is actually an NNPDF 3.0 NLO type!
          if (wgtid == 2000){
            defaultMemberZeroWeight = wgtval; found_defaultMemberZeroWeight=true;
            if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF30) && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNLO)){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
          }
          else if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF30) && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNLO)){
            // Do not assign defaultMemberZeroWeight or defaultWeightScale. This will be handled correctly later on via LHEWeight[0].
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useMC;
          }
        }
        else{
          if (wgtid == 2000){
            defaultMemberZeroWeight = wgtval; found_defaultMemberZeroWeight=true;
            if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
          }
          else if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){
            // Do not assign defaultMemberZeroWeight or defaultWeightScale. This will be handled correctly later on via LHEWeight[0].
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useHessian;
          }
        }
      }
      else if (weightstype == powheg && 2001 <= wgtid && wgtid <= 2102){ // These are more NNLO variations
        if (specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1){ // How fun is that?! This special case is actually an NNPDF 3.0 NLO type!
          if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF30) && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNLO)){
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useMC;
          }
        }
        else{
          if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useHessian;
          }
        }
      }
      else if (weightstype == powheg && 2103 <= wgtid && wgtid <= 2111){ /*do nothing, these are other a_s variations*/ }
      else if (weightstype == powheg && 1500 <= wgtid && wgtid <= 1602){ // These are the NNPDF 3.0 NLO pdf and MC variations
        if (pdfChoice==tryNNPDF30 && orderChoice==tryNLO){
          if (wgtid == 1500){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
          else{
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useMC;
          }
        }
      }

      else if (weightstype == powheg && wgtid == 1700){ // This is the NNLO pdf for NNPDF 3.0, no variations
        if (pdfChoice==tryNNPDF30 && (orderChoice==keepDefaultQCDOrder || orderChoice==tryNNLO)){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (weightstype == powheg && wgtid == 1850){ // This is the LO pdf for NNPDF 3.1 with a_s=0.130, no variations
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && orderChoice==tryLO){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (weightstype == powheg && wgtid == 1950){ // This is the LO pdf for NNPDF 3.0 with a_s=0.130, no variations
        if (pdfChoice==tryNNPDF30 && orderChoice==tryLO){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
      }
      else if (weightstype == powheg && (wgtid == 1800 || wgtid == 1900)){/*do nothing, these are LO pdfs*/ }
      else if (weightstype == powheg && wgtid >= 4000)                  {/*do nothing, these are other various weights*/ }

      else if (weightstype == powheg && 3000 <= wgtid && wgtid <= 3102){ // These are the NNPDF3.1 NLO and variations
        if ((pdfChoice==keepDefaultPDF || pdfChoice==tryNNPDF31) && orderChoice==tryNLO){
          if (wgtid == 3000){ defaultWeightScale = wgtval; found_defaultWeightScale = true; }
          else{
            LHEPDFVariationWgt.push_back(wgtval);
            if (pdfVarMode == useNone) pdfVarMode = useHessian;
          }
        }
      }

      else throw cms::Exception("LHEWeights") << "Don't know what to do with alternate weight id = " << wgtid << "(weightstype == " << weightstype << ")";
    }
    else throw cms::Exception("LHEWeights") << "Unknown year " << year;
    it_weightid++;
  }
  // Handle LO samples with 101 or 103 alternative weights
  if (year == 2016 && (LHEPDFVariationWgt.size()==101 || LHEPDFVariationWgt.size()==103)) LHEPDFVariationWgt.erase(LHEPDFVariationWgt.begin(), LHEPDFVariationWgt.begin()+1);
  /*
  static bool firsttimeaccount=true;
  if (firsttimeaccount){
    edm::LogWarning warning("PDFWeightsInfo");
    warning << "There are 1 + " << LHEWeight.size() << " + " << LHEPDFVariationWgt.size() << " weights.";
    firsttimeaccount=false;
  }
  */

  /* NO LOOP OVER WEIGHTS BEYOND THIS POINT */

  // Check if member 0 weight exists. This is the value from the default PDF with which the MC is simulated.
  // It might not exist if the alternative for the default weights skip the central member (happens in some samples).
  // In that case, if muR1_muF1 weight exists, it should be the same as member 0, so pick that one.
  // Otherwise, most likely no alternative weights exist, so just pick the original weight.
  if (!found_defaultMemberZeroWeight){
    if (!LHEWeight.empty()) defaultMemberZeroWeight = LHEWeight.at(0);
    else defaultMemberZeroWeight = LHEOriginalWeight;
  }
  if (LHEOriginalWeight==1.f && defaultMemberZeroWeight!=1.f){ // This is a Pythia bug in 2016 samples
    LHEOriginalWeight = defaultMemberZeroWeight;
  }
  // defaultWeightScale is the central weight of the alternative PDFChoice or QCDOrderChoice weights.
  // It serves the same purpose as defaultMemberZeroWeight, but it belongs to the PDF variations that are actually picked by this program.
  // It might not exist if the deafult PDF and QCD order is chosen and member 0 of the variations for this set is skipped (i.e. found_defaultMemberZeroWeight==false).
  // Therefore, if it doesn't exist, assign defaultMemberZeroWeight value to it.
  if (!found_defaultWeightScale) defaultWeightScale = defaultMemberZeroWeight;
  if (!LHEWeight.empty()){ // Divide LHE weights by LHEWeight[0]
    float lheweight0 = LHEWeight.at(0);
    for (float& w:LHEWeight) w = safeDivide(w, lheweight0);
  }
  // Divide PDF variation weights by the member 0 weight that corresponds to thte variations
  for (float& w:LHEPDFVariationWgt) w = safeDivide(w, defaultWeightScale);
  // Divide defaultWeightScale by defaultMemberZeroWeight to get a true scale factor
  defaultWeightScale = safeDivide(defaultWeightScale, defaultMemberZeroWeight);

  static bool warning_LHEWeight_NoLHEPDFVars=false;
  if (!warning_LHEWeight_NoLHEPDFVars){
    if (year == 2016 && !weights.empty() && !(LHEWeight.size() == 9 && (LHEPDFVariationWgt.size()==100 || LHEPDFVariationWgt.size()==102))){
      edm::LogWarning warning("LHEWeights");
      warning << "For " << year << " MC, PDF choice " << pdfChoice << ", and QCD order " << orderChoice << ", expect to find either\n"
        << "\t- no alternative LHE weights, or\n"
        << "\t- all of the following:\n"
        << "\t\t- muR and muF variations (found " << LHEWeight.size() << " / 9 of them)\n"
        << "\t\t- PDF (+ alpha_s) weight variations (found " << LHEPDFVariationWgt.size() << " / 100 (or 102) of them)";
    }
    else if ((year == 2017 || year == 2018) && !weights.empty() && orderChoice!=tryLO && !(LHEWeight.size() == 9 && (LHEPDFVariationWgt.size()==100 || LHEPDFVariationWgt.size()==102))){
      edm::LogWarning warning("LHEWeights");
      warning << "For " << year << " MC, PDF choice " << pdfChoice << ", and QCD order " << orderChoice << ", expect to find either\n"
        << "\t- no alternative LHE weights, or\n"
        << "\t- all of the following:\n"
        << "\t\t- muR and muF variations (found " << LHEWeight.size() << " / 9 of them)\n"
        << "\t\t- PDF (+ alpha_s) weight variations (found " << LHEPDFVariationWgt.size() << " / 100 (or 102) of them)";
    }
    warning_LHEWeight_NoLHEPDFVars=true;
  }

  static bool warning_LHEWeight_transposed=false;
  if (weightstype == madgraph_1000offset && LHEWeight.size() == 9){ // but not 0offset!  0offset does it the same way as powheg
    if (year == 2016 || year == 2017 || year == 2018){
      if (!warning_LHEWeight_transposed){
        edm::LogWarning warning("LHEWeightTransposed");
        warning << "For " << year << " MC and weight type " << weightstype << ", the matrix of LHE muR/muF weights is transposed.";
        warning_LHEWeight_transposed=true;
      }
      LHEWeight={
        LHEWeight[0], LHEWeight[3], LHEWeight[6],
        LHEWeight[1], LHEWeight[4], LHEWeight[7],
        LHEWeight[2], LHEWeight[5], LHEWeight[8]
      };
    }
    else throw cms::Exception("LHEWeightTransposed") << "Unknown year " << year;
  }

  if (LHEPDFVariationWgt.size() > 101){
    std::vector<float>::iterator firstalphasweight = LHEPDFVariationWgt.begin() + 100; //= iterator to LHEPDFVariationWgt[100] = 101st entry
    LHEPDFAlphaSMZWgt.assign(firstalphasweight, LHEPDFVariationWgt.end());
    LHEPDFVariationWgt.erase(firstalphasweight, LHEPDFVariationWgt.end());
  }

  // Find the proper PDF and alphas(mZ) variations
  if (!LHEPDFVariationWgt.empty() && defaultWeightScale!=0.f){
    if (pdfVarMode == useMC){
      // LHE weight suppresssion not needed
      std::sort(LHEPDFVariationWgt.begin(), LHEPDFVariationWgt.end(), compareAbsIsLess);
      LHEWeight_PDFVariationUpDn.push_back(findNearestOneSigma(1., 1, LHEPDFVariationWgt));
      LHEWeight_PDFVariationUpDn.push_back(findNearestOneSigma(1., -1, LHEPDFVariationWgt));
      if (LHEPDFAlphaSMZWgt.size()>1){
        float const& asdn = LHEPDFAlphaSMZWgt.at(0);
        float const& asup = LHEPDFAlphaSMZWgt.at(1);
        // Rescale alphas(mZ) variations from 0.118+-0.001 to 0.118+-0.0015
        // 0.001 is valid both for the alternates used in 2016 and for NNPDF30_nlo_nf_5_pdfas in 2017/18
        LHEWeight_AsMZUpDn.push_back(1.f + (asup-1.f)*1.5f);
        LHEWeight_AsMZUpDn.push_back(1.f + (asdn-1.f)*1.5f);
      }
    }
    else if (pdfVarMode == useHessian){
      //https://arxiv.org/pdf/1706.00428v2.pdf page 88
      //we are working with NNPDF31_nlo_hessian_pdfas
      //see the routine starting on line 101 of PDFSet.cc in LHAPDF-6.2.1
      //based on NNPDF31_nlo_hessian_pdfas.info, which confirms that errorType() is symmhessian+as
      float error = 0.f;
      for (const auto& wt:LHEPDFVariationWgt) error += pow(wt - 1.f, 2);
      error = sqrt(error);
      
      LHEWeight_PDFVariationUpDn ={ (1.f + error), (1.f - error) };

      if (LHEPDFAlphaSMZWgt.size()>1){
        float asdn = LHEPDFAlphaSMZWgt.at(0);
        float asup = LHEPDFAlphaSMZWgt.at(1);
        // Rescale alphas(mZ) variations from 0.118+-0.002 to 0.118+-0.0015
        //                           Note this number (^) is different than 2016!
        LHEWeight_AsMZUpDn.push_back(1.f + (asup-1.f)*0.75);
        LHEWeight_AsMZUpDn.push_back(1.f + (asdn-1.f)*0.75);
      }
    }
  }
}

bool LHEWeightHandler::compareAbsIsLess(float val1, float val2){ return std::abs(val1) < std::abs(val2); }

float LHEWeightHandler::findNearestOneSigma(float ref, int lowhigh, std::vector<float> const& wgt_array){
  int nrep = wgt_array.size();
  int pos_low=-1, pos_high=nrep;
  for (int irep=0; irep<nrep; irep++){ // Assume ordered from low to high in absolute value
    float tmp = fabs(wgt_array.at(irep));
    if (fabs(tmp)<fabs(ref)) pos_low=irep;
    else if (fabs(tmp)>fabs(ref)){
      pos_high=irep;
      break;
    }
  }
  if (lowhigh==-1 && nrep>0){ // Low
    float ninst = pos_low+1;
    if (ninst==0) return ref;
    int nprogress = (ninst*0.6827+0.5); // truncate down
    int taken = ninst-nprogress;
    return wgt_array.at(taken);
  }
  else if (lowhigh==1 && nrep>0){
    float ninst = nrep-pos_high;
    if (ninst==0) return ref;
    int nprogress = (ninst*0.6827+0.5); // truncate down
    int taken = nprogress+pos_high-1;
    return wgt_array.at(taken);
  }
  else return ref;
}

void LHEWeightHandler::suppressLargeWeights(std::vector<float>& wgt_array){
  constexpr float c_wgt_thr=5.f;
  constexpr float c_trunc_min=0.35f;
  constexpr float c_trunc_max=0.7f;
  if (wgt_array.empty()) return;

  std::vector<float> weightListNP[2];
  for (float const& weight:wgt_array){
    if (weight<0.) weightListNP[0].push_back(fabs(weight));
    else if (weight>0.) weightListNP[1].push_back(fabs(weight));
  }

  float weight_meanlow=0;
  float weight_meanhigh=0;
  float weight_thrlow=0;
  float weight_thrhigh=0;
  for (unsigned int inp=0; inp<2; inp++){
    std::vector<float>& weightList = weightListNP[inp];
    float& weight_mean = (inp==0 ? weight_meanlow : weight_meanhigh);
    float& weight_thr = (inp==0 ? weight_thrlow : weight_thrhigh);
    const unsigned int wsize=weightList.size();
    if (wsize>0){
      std::sort(weightList.begin(), weightList.end());
      const unsigned int ibegin=float(wsize)*c_trunc_min;
      const unsigned int iend=float(wsize)*c_trunc_max+1;
      std::vector<float>::const_iterator it_begin=weightList.cbegin()+ibegin;
      std::vector<float>::const_iterator it_end=weightList.cbegin()+iend;
      // Find mean and rms
      for (std::vector<float>::const_iterator it=it_begin; it!=it_end; it++){
        float const& wgt = *it;
        weight_mean += wgt;
      }
      weight_mean /= float(iend-ibegin);
      unsigned int nhigh=0;
      for (std::vector<float>::const_iterator it=it_begin; it!=it_end; it++){
        float const& wgt = *it;
        if (wgt>=weight_mean){ weight_thr += pow(wgt, 2); nhigh++; }
      }
      if (nhigh>0){
        weight_thr /= float(nhigh);
        weight_thr = sqrt(fabs(weight_thr-pow(weight_mean, 2)));
        weight_thr = c_wgt_thr*weight_thr + weight_mean;
      }
      else weight_thr = c_wgt_thr*weight_mean;
      //cout << "Weight threshold [" << inp << "]: " << weight_mean << " +" << weight_thr << "/-0" << endl;
    }
  }
  weight_thrlow = -weight_thrlow;
  weight_meanlow = -weight_meanlow;

  {
    float const diffthrlow = weight_thrlow - weight_meanlow;
    float const diffthrhigh = weight_thrhigh - weight_meanhigh;
    for (float& weight:wgt_array){
      if (weight==0.) continue;
      if (weight<weight_thrlow){
        float weightdiff = weight - weight_meanlow;
        weightdiff = pow(diffthrlow, 2)/weightdiff;
        weight = weightdiff + weight_meanlow;
      }
      else if (weight>weight_thrhigh){
        float weightdiff = weight - weight_meanhigh;
        weightdiff = pow(diffthrhigh, 2)/weightdiff;
        weight = weightdiff + weight_meanhigh;
      }
    }
  }

}
