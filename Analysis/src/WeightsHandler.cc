#include <cassert>
#include <iterator>
#include <algorithm>
#include <utility>
#include "WeightsHandler.h"
#include "SampleHelpers.h"
#include "FrameworkVariables.hh"
#include "FrameworkSet.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


WeightsHandler::WeightsHandler() :
  IvyBase(),
  use2016Scheme(true),
  product(nullptr),
  weightHandler_DefaultPDF(nullptr),
  weightHandler_2016(nullptr)
{}

bool WeightsHandler::constructWeights(){
  static bool first_event=true;

  clear();
  if (!currentTree){
    if (verbosity>=TVar::ERROR) MELAerr << "WeightsHandler::constructWeights: Current tree is null!" << endl;
    return false;
  }
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(currentTree);
  if (!fwktree){
    if (verbosity>=TVar::ERROR) MELAerr << "WeightsHandler::constructWeights: Current tree is not derived from a FrameworkTree class!" << endl;
    return false;
  }

  std::vector<float>* const* weightsPtr = valVfloats[_genweights_];
  std::vector<float> const* weights = (weightsPtr ? *weightsPtr : nullptr);
  std::vector<std::string>* const* weightIdsPtr = valVstrings[_genweightIDs_];
  std::vector<std::string> const* weightIds = (weightIdsPtr ? *weightIdsPtr : nullptr);

  product = new ProductType_t;
  auto& productExtras = product->extras;
  bool hasNewWeights = currentTree->branchExists(_genHEPMCweight_); // Should have been already booked by now
  if (hasNewWeights){
    if (first_event && verbosity>=TVar::INFO) MELAout << "WeightsHandler::constructWeights: Case hasNewWeights=true" << endl;

    float const* genHEPMCweight = valfloats[_genHEPMCweight_];
    float const* LHEweight_QCDscale_muR1_muF2 = valfloats[_LHEweight_QCDscale_muR1_muF2_];
    float const* LHEweight_QCDscale_muR1_muF0p5 = valfloats[_LHEweight_QCDscale_muR1_muF0p5_];
    float const* LHEweight_QCDscale_muR2_muF1 = valfloats[_LHEweight_QCDscale_muR2_muF1_];
    float const* LHEweight_QCDscale_muR0p5_muF1 = valfloats[_LHEweight_QCDscale_muR0p5_muF1_];
    float const* LHEweight_PDFVariation_Up = valfloats[_LHEweight_PDFVariation_Up_];
    float const* LHEweight_PDFVariation_Dn = valfloats[_LHEweight_PDFVariation_Dn_];
    float const* LHEweight_AsMZ_Up = valfloats[_LHEweight_AsMZ_Up_];
    float const* LHEweight_AsMZ_Dn = valfloats[_LHEweight_AsMZ_Dn_];
    float const* genHEPMCweight_2016 = nullptr;
    float const* LHEweight_PDFVariation_Up_2016 = nullptr;
    float const* LHEweight_PDFVariation_Dn_2016 = nullptr;
    float const* LHEweight_AsMZ_Up_2016 = nullptr;
    float const* LHEweight_AsMZ_Dn_2016 = nullptr;
    if (SampleHelpers::theDataYear > 2016 && use2016Scheme){
      if (first_event && verbosity>=TVar::INFO) MELAout << "\t- Linking also 2016-like weights" << endl;

      genHEPMCweight_2016 = valfloats[_genHEPMCweight_2016_];
      LHEweight_PDFVariation_Up_2016 = valfloats[_LHEweight_PDFVariation_Up_2016_];
      LHEweight_PDFVariation_Dn_2016 = valfloats[_LHEweight_PDFVariation_Dn_2016_];
      LHEweight_AsMZ_Up_2016 = valfloats[_LHEweight_AsMZ_Up_2016_];
      LHEweight_AsMZ_Dn_2016 = valfloats[_LHEweight_AsMZ_Dn_2016_];
    }
    // If 2016-like weights exist and 2016 PDF variations are present...
    bool use2016 = (genHEPMCweight_2016!=nullptr);
    if (use2016) use2016 &= (((*genHEPMCweight_2016)==0.) || !((*LHEweight_PDFVariation_Up_2016)==1.f && (*LHEweight_PDFVariation_Dn_2016)==1.f));
    if (use2016){
      if (first_event && verbosity>=TVar::INFO) MELAout << "\t- Using 2016-like weights" << endl;

      productExtras.wgt_central = (*genHEPMCweight_2016);
      productExtras.wgt_PDFVariationUp = productExtras.wgt_central * (*LHEweight_PDFVariation_Up_2016);
      productExtras.wgt_PDFVariationDn = productExtras.wgt_central * (*LHEweight_PDFVariation_Dn_2016);
      productExtras.wgt_AsMZUp = productExtras.wgt_central * (*LHEweight_AsMZ_Up_2016);
      productExtras.wgt_AsMZDn = productExtras.wgt_central * (*LHEweight_AsMZ_Dn_2016);

      productExtras.wgt_central_default = (*genHEPMCweight);
      productExtras.wgt_PDFVariationUp_default = productExtras.wgt_central_default * (*LHEweight_PDFVariation_Up);
      productExtras.wgt_PDFVariationDn_default = productExtras.wgt_central_default * (*LHEweight_PDFVariation_Dn);
    }
    else{
      if (first_event && verbosity>=TVar::INFO) MELAout << "\t- Using default weights" << endl;

      productExtras.wgt_central_default = productExtras.wgt_central = (*genHEPMCweight);
      productExtras.wgt_PDFVariationUp_default = productExtras.wgt_PDFVariationUp = productExtras.wgt_central * (*LHEweight_PDFVariation_Up);
      productExtras.wgt_PDFVariationDn_default = productExtras.wgt_PDFVariationDn = productExtras.wgt_central * (*LHEweight_PDFVariation_Dn);
      productExtras.wgt_AsMZUp = productExtras.wgt_central * (*LHEweight_AsMZ_Up);
      productExtras.wgt_AsMZDn = productExtras.wgt_central * (*LHEweight_AsMZ_Dn);
    }
    productExtras.wgt_muF2 = productExtras.wgt_central * (*LHEweight_QCDscale_muR1_muF2);
    productExtras.wgt_muF0p5 = productExtras.wgt_central * (*LHEweight_QCDscale_muR1_muF0p5);
    productExtras.wgt_muR2 = productExtras.wgt_central * (*LHEweight_QCDscale_muR2_muF1);
    productExtras.wgt_muR0p5 = productExtras.wgt_central * (*LHEweight_QCDscale_muR0p5_muF1);
  }
  else if (weights && weightIds){
    if (first_event && verbosity>=TVar::INFO) MELAout << "WeightsHandler::constructWeights: Case hasNewWeights=false" << endl;

    int const& year = SampleHelpers::theDataYear;
    if (!weightHandler_DefaultPDF) weightHandler_DefaultPDF = new LHEWeightHandler(year, LHEWeightHandler::keepDefaultPDF, LHEWeightHandler::keepDefaultQCDOrder);
    if (year!=2016 && use2016Scheme && !weightHandler_2016) weightHandler_2016 = new LHEWeightHandler(year, LHEWeightHandler::tryNNPDF30, LHEWeightHandler::tryNLO);

    float const* genHEPMCweight_old = valfloats[_genHEPMCweight_old_];

    auto const& exceptionalCases = fwktree->getOptions().getLHEExceptionalCases();
    bool use2016 = false;
    if (weightHandler_2016){
      if (first_event && verbosity>=TVar::INFO) MELAout << "\t- Constructing 2016-like weights" << endl;

      weightHandler_2016->set_specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1(exceptionalCases.specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1);
      weightHandler_2016->set_specialPDF_NNPDF31_NNLO_as_0118_nf_4(exceptionalCases.specialPDF_NNPDF31_NNLO_as_0118_nf_4);
      weightHandler_2016->set_specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1(exceptionalCases.specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1);
      weightHandler_2016->extract(*genHEPMCweight_old, *weights, *weightIds);

      float pdfup = weightHandler_DefaultPDF->getLHEWeight_PDFVariationUpDn(1, 1.);
      float pdfdn = weightHandler_DefaultPDF->getLHEWeight_PDFVariationUpDn(-1, 1.);
      float cvwgt = weightHandler_2016->getLHEOriginalWeight() * weightHandler_2016->getWeightRescale();
      use2016 = ((cvwgt==0.) || !(pdfup==1.f && pdfdn==1.f));
    }
    if (weightHandler_DefaultPDF){
      if (first_event && verbosity>=TVar::INFO) MELAout << "\t- Constructing default weights" << endl;

      weightHandler_DefaultPDF->set_specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1(exceptionalCases.specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1);
      weightHandler_DefaultPDF->set_specialPDF_NNPDF31_NNLO_as_0118_nf_4(exceptionalCases.specialPDF_NNPDF31_NNLO_as_0118_nf_4);
      weightHandler_DefaultPDF->set_specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1(exceptionalCases.specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1);
      weightHandler_DefaultPDF->extract(*genHEPMCweight_old, *weights, *weightIds);
    }
    if (use2016){
      if (first_event && verbosity>=TVar::INFO) MELAout << "\t- Using 2016-like weights" << endl;

      productExtras.wgt_central = weightHandler_2016->getLHEOriginalWeight() * weightHandler_2016->getWeightRescale();
      productExtras.wgt_muF2 = productExtras.wgt_central * weightHandler_2016->getLHEWeight(1, 1.);
      productExtras.wgt_muF0p5 = productExtras.wgt_central * weightHandler_2016->getLHEWeight(2, 1.);
      productExtras.wgt_muR2 = productExtras.wgt_central * weightHandler_2016->getLHEWeight(3, 1.);
      productExtras.wgt_muR0p5 = productExtras.wgt_central * weightHandler_2016->getLHEWeight(6, 1.);
      productExtras.wgt_PDFVariationUp = productExtras.wgt_central * weightHandler_2016->getLHEWeight_PDFVariationUpDn(1, 1.);
      productExtras.wgt_PDFVariationDn = productExtras.wgt_central * weightHandler_2016->getLHEWeight_PDFVariationUpDn(-1, 1.);
      productExtras.wgt_AsMZUp = productExtras.wgt_central * weightHandler_2016->getLHEWeigh_AsMZUpDn(1, 1.);
      productExtras.wgt_AsMZDn = productExtras.wgt_central * weightHandler_2016->getLHEWeigh_AsMZUpDn(-1, 1.);

      if (weightHandler_DefaultPDF){
        productExtras.wgt_central_default = weightHandler_DefaultPDF->getLHEOriginalWeight() * weightHandler_DefaultPDF->getWeightRescale();
        productExtras.wgt_PDFVariationUp_default = productExtras.wgt_central_default * weightHandler_DefaultPDF->getLHEWeight_PDFVariationUpDn(1, 1.);
        productExtras.wgt_PDFVariationDn_default = productExtras.wgt_central_default * weightHandler_DefaultPDF->getLHEWeight_PDFVariationUpDn(-1, 1.);
      }
    }
    else if (weightHandler_DefaultPDF){
      if (first_event && verbosity>=TVar::INFO) MELAout << "\t- Using default weights" << endl;

      productExtras.wgt_central_default = productExtras.wgt_central = weightHandler_DefaultPDF->getLHEOriginalWeight() * weightHandler_DefaultPDF->getWeightRescale();
      productExtras.wgt_muF2 = productExtras.wgt_central * weightHandler_DefaultPDF->getLHEWeight(1, 1.);
      productExtras.wgt_muF0p5 = productExtras.wgt_central * weightHandler_DefaultPDF->getLHEWeight(2, 1.);
      productExtras.wgt_muR2 = productExtras.wgt_central * weightHandler_DefaultPDF->getLHEWeight(3, 1.);
      productExtras.wgt_muR0p5 = productExtras.wgt_central * weightHandler_DefaultPDF->getLHEWeight(6, 1.);
      productExtras.wgt_PDFVariationUp_default = productExtras.wgt_PDFVariationUp = productExtras.wgt_central * weightHandler_DefaultPDF->getLHEWeight_PDFVariationUpDn(1, 1.);
      productExtras.wgt_PDFVariationDn_default = productExtras.wgt_PDFVariationDn = productExtras.wgt_central * weightHandler_DefaultPDF->getLHEWeight_PDFVariationUpDn(-1, 1.);
      productExtras.wgt_AsMZUp = productExtras.wgt_central * weightHandler_DefaultPDF->getLHEWeigh_AsMZUpDn(1, 1.);
      productExtras.wgt_AsMZDn = productExtras.wgt_central * weightHandler_DefaultPDF->getLHEWeigh_AsMZUpDn(-1, 1.);
    }
  }
  else{
    if (!weights && verbosity>=TVar::ERROR) MELAerr << "WeightsHandler::constructWeights: The weight list is null!" << endl;
    if (!weightIds && verbosity>=TVar::ERROR) MELAerr << "WeightsHandler::constructWeights: The weight id list is null!" << endl;
    assert(0);
    return false;
  }

  // Get the parton showe variations (i.e. ISR*FSR with scale variation of x/ 4)
  if (weights && weightIds && !weights->empty()){
    std::vector<float> genwgtvars;
    std::vector<float>::const_iterator it_copy_begin = weights->cend();
    auto it_wgts = weights->cbegin();
    for (std::string const& wgtid:(*weightIds)){
      if (wgtid==""){ it_copy_begin=it_wgts; break; }
      it_wgts++;
    }
    if (it_copy_begin!=weights->cend()) std::copy(it_copy_begin, weights->cend(), std::back_inserter(genwgtvars));
    if (genwgtvars.size()>1){
      if (genwgtvars.size()<14 || genwgtvars.at(0) != genwgtvars.at(1)){
        if (verbosity>=TVar::ERROR) MELAerr << "WeightsHandler::constructWeights: Expected to find 1 gen weight, or >=14 with the first two the same, found " << genwgtvars.size() << ": " << genwgtvars << endl;
        assert(0);
        return false;
      }
      float const& nominal = genwgtvars.at(0);
      /*
      float PythiaWeight_isr_muRoneoversqrt2 = (nominal!=0.f ? genwgtvars[2]/nominal : 0.f);
      float PythiaWeight_fsr_muRoneoversqrt2 = (nominal!=0.f ? genwgtvars[3]/nominal : 0.f);
      float PythiaWeight_isr_muRsqrt2 = (nominal!=0.f ? genwgtvars[4]/nominal : 0.f);
      float PythiaWeight_fsr_muRsqrt2 = (nominal!=0.f ? genwgtvars[5]/nominal : 0.f);
      float PythiaWeight_isr_muR0p5 = (nominal!=0.f ? genwgtvars[6]/nominal : 0.f);
      float PythiaWeight_fsr_muR0p5 = (nominal!=0.f ? genwgtvars[7]/nominal : 0.f);
      float PythiaWeight_isr_muR2 = (nominal!=0.f ? genwgtvars[8]/nominal : 0.f);
      float PythiaWeight_fsr_muR2 = (nominal!=0.f ? genwgtvars[9]/nominal : 0.f);
      */
      float PythiaWeight_isr_muR0p25 = (nominal!=0.f ? genwgtvars[10]/nominal : 0.f);
      float PythiaWeight_fsr_muR0p25 = (nominal!=0.f ? genwgtvars[11]/nominal : 0.f);
      float PythiaWeight_isr_muR4 = (nominal!=0.f ? genwgtvars[12]/nominal : 0.f);
      float PythiaWeight_fsr_muR4 = (nominal!=0.f ? genwgtvars[13]/nominal : 0.f);
      productExtras.wgt_PSUp = productExtras.wgt_central * PythiaWeight_isr_muR4*PythiaWeight_fsr_muR4;
      productExtras.wgt_PSDn = productExtras.wgt_central * PythiaWeight_isr_muR0p25*PythiaWeight_fsr_muR0p25;
      productExtras.wgt_ISRUp = productExtras.wgt_central * PythiaWeight_isr_muR4;
      productExtras.wgt_ISRDn = productExtras.wgt_central * PythiaWeight_isr_muR0p25;
      productExtras.wgt_FSRUp = productExtras.wgt_central * PythiaWeight_fsr_muR4;
      productExtras.wgt_FSRDn = productExtras.wgt_central * PythiaWeight_fsr_muR0p25;
    }
  }

  if (first_event) first_event=false;

  return true;
}

void WeightsHandler::bookBranches(BaseTree* tree){
  if (!tree || !tree->isValid()) return;
  FrameworkTree* fwktree = dynamic_cast<FrameworkTree*>(tree);
  if (!fwktree) return;
  SampleHelpers::setupUsingOptions(fwktree->getOptions());

  bool hasNewWeights = (SampleHelpers::branchExists(fwktree->getSelectedTree(), _genHEPMCweight_) || SampleHelpers::aliasExists(fwktree->getSelectedTree(), _genHEPMCweight_));
  if (hasNewWeights){
    fwktree->bookEDMBranch<float>(_genHEPMCweight_, 0);
    fwktree->bookEDMBranch<float>(_LHEweight_QCDscale_muR1_muF2_, 0);
    fwktree->bookEDMBranch<float>(_LHEweight_QCDscale_muR1_muF0p5_, 0);
    fwktree->bookEDMBranch<float>(_LHEweight_QCDscale_muR2_muF1_, 0);
    fwktree->bookEDMBranch<float>(_LHEweight_QCDscale_muR0p5_muF1_, 0);
    fwktree->bookEDMBranch<float>(_LHEweight_PDFVariation_Up_, 0);
    fwktree->bookEDMBranch<float>(_LHEweight_PDFVariation_Dn_, 0);
    fwktree->bookEDMBranch<float>(_LHEweight_AsMZ_Up_, 0);
    fwktree->bookEDMBranch<float>(_LHEweight_AsMZ_Dn_, 0);

    this->addConsumed<float>(_genHEPMCweight_);
    this->addConsumed<float>(_LHEweight_QCDscale_muR1_muF2_);
    this->addConsumed<float>(_LHEweight_QCDscale_muR1_muF0p5_);
    this->addConsumed<float>(_LHEweight_QCDscale_muR2_muF1_);
    this->addConsumed<float>(_LHEweight_QCDscale_muR0p5_muF1_);
    this->addConsumed<float>(_LHEweight_PDFVariation_Up_);
    this->addConsumed<float>(_LHEweight_PDFVariation_Dn_);
    this->addConsumed<float>(_LHEweight_AsMZ_Up_);
    this->addConsumed<float>(_LHEweight_AsMZ_Dn_);

    this->defineConsumedSloppy(_genHEPMCweight_);
    this->defineConsumedSloppy(_LHEweight_QCDscale_muR1_muF2_);
    this->defineConsumedSloppy(_LHEweight_QCDscale_muR1_muF0p5_);
    this->defineConsumedSloppy(_LHEweight_QCDscale_muR2_muF1_);
    this->defineConsumedSloppy(_LHEweight_QCDscale_muR0p5_muF1_);
    this->defineConsumedSloppy(_LHEweight_PDFVariation_Up_);
    this->defineConsumedSloppy(_LHEweight_PDFVariation_Dn_);
    this->defineConsumedSloppy(_LHEweight_AsMZ_Up_);
    this->defineConsumedSloppy(_LHEweight_AsMZ_Dn_);

    if (SampleHelpers::theDataYear > 2016 && use2016Scheme){ // Not needed for 2016
      fwktree->bookEDMBranch<float>(_genHEPMCweight_2016_, 0);
      fwktree->bookEDMBranch<float>(_LHEweight_PDFVariation_Up_2016_, 0);
      fwktree->bookEDMBranch<float>(_LHEweight_PDFVariation_Dn_2016_, 0);
      fwktree->bookEDMBranch<float>(_LHEweight_AsMZ_Up_2016_, 0);
      fwktree->bookEDMBranch<float>(_LHEweight_AsMZ_Dn_2016_, 0);

      this->addConsumed<float>(_genHEPMCweight_2016_);
      this->addConsumed<float>(_LHEweight_PDFVariation_Up_2016_);
      this->addConsumed<float>(_LHEweight_PDFVariation_Dn_2016_);
      this->addConsumed<float>(_LHEweight_AsMZ_Up_2016_);
      this->addConsumed<float>(_LHEweight_AsMZ_Dn_2016_);

      this->defineConsumedSloppy(_genHEPMCweight_2016_);
      this->defineConsumedSloppy(_LHEweight_PDFVariation_Up_2016_);
      this->defineConsumedSloppy(_LHEweight_PDFVariation_Dn_2016_);
      this->defineConsumedSloppy(_LHEweight_AsMZ_Up_2016_);
      this->defineConsumedSloppy(_LHEweight_AsMZ_Dn_2016_);
    }
  }
  else{
    fwktree->bookEDMBranch<float>(_genHEPMCweight_old_, 0);
    this->addConsumed<float>(_genHEPMCweight_old_);
  }
  // These weights are also used for Pythia variations, so keep them as common
  fwktree->bookEDMBranch<std::vector<float>*>(_genweights_, nullptr);
  fwktree->bookEDMBranch<std::vector<std::string>*>(_genweightIDs_, nullptr);

  this->addConsumed<std::vector<float>*>(_genweights_);
  this->addConsumed<std::vector<std::string>*>(_genweightIDs_);
}

bool WeightsHandler::recordWeights(SimpleEntry& entry, float multiplier) const{
  if (product){
    for (int iw=WeightVariables::wCentral; iw<WeightVariables::nWeightTypes; iw++){
      TString s = WeightVariables::getWeightName((WeightVariables::WeightType)iw);
      float v = multiplier * product->extras.getWeight((WeightVariables::WeightType)iw);
      entry.setNamedVal<float>(s, v);
    }
    return true;
  }
  else{
    MELAerr << "WeightsHandler::recordWeights: Product is invalid! You should probably run WeightsHandler::constructWeights() first." << endl;
    for (int iw=WeightVariables::wCentral; iw<WeightVariables::nWeightTypes; iw++){
      TString s = WeightVariables::getWeightName((WeightVariables::WeightType)iw);
      entry.setNamedVal<float>(s, 0);
    }
    return false;
  }
}

