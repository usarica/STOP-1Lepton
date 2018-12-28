#ifndef LHEWEIGHTHADNLER_H
#define LHEWEIGHTHADNLER_H

#include <vector>
#include <string>


class LHEWeightHandler{
public:
  enum PDFChoice{
    keepDefaultPDF=0,
    tryNNPDF30,
    tryNNPDF31
  };
  enum QCDOrderChoice{
    keepDefaultQCDOrder=0,
    tryLO,
    tryNLO,
    tryNNLO
  };

  LHEWeightHandler(int year_, LHEWeightHandler::PDFChoice pdfChoice_, LHEWeightHandler::QCDOrderChoice orderChoice_);
  virtual ~LHEWeightHandler();

  void extract(float wgt_central, std::vector<float> const& weights, std::vector<std::string> const& weightIds);
  void clear();

  float const& getLHEOriginalWeight() const; // Weight written in the <event> block, supposed to = genhepmcweight if no Pythia reweighting is done
  float const& getMemberZeroWeight() const; // Weight from POWHEG before JHUGen reweighting, taken from alternate weight 1001.  If there are no alternate weights this is the same as the LHEOriginalWeight
  float const& getWeightRescale() const; // Nominal weight should be getLHEOriginalWeight() * getWeightRescale()
  float getLHEWeight(unsigned int whichWeight, float defaultValue=1) const; // = {Weights written in LHE weight variations} / getLHEOriginalWeight()
  float getLHEWeight_PDFVariationUpDn(int whichUpDn, float defaultValue=1) const; // = {Weights written in LHE weight variations} / getLHEOriginalWeight()
  float getLHEWeigh_AsMZUpDn(int whichUpDn, float defaultValue=1) const; // = {Weights written in LHE weight variations} / getLHEOriginalWeight()

  // Misc. functions needed for ordering the PDF weights
  static bool compareAbsIsLess(float val1, float val2);
  static void suppressLargeWeights(std::vector<float>& wgt_array);
  static float findNearestOneSigma(float ref, int lowhigh, std::vector<float> const& wgt_array);
  static float safeDivide(float numerator, float denominator);

  // Include functions for exceptional cases from a separate file here
  // so that these functions are members of the LHEWeightHandler class.
#include "LHEWeightHandler_ExceptionalCases.h"

protected:
  enum PDFVariationMode{
    useNone=0,
    useMC=1,
    useHessian=2,
    useAlternativePDFs=3
  };
  enum AlternateWeightsType{
    unknown,
    powheg,
    madgraph_0offset,
    madgraph_1000offset
  };

  // Year of the MC
  const int year;

  // These options influence which pdf sets/members are selected
  PDFChoice pdfChoice;
  QCDOrderChoice orderChoice;

  float defaultMemberZeroWeight;
  float defaultWeightScale;
  float LHEOriginalWeight;
  std::vector<float> LHEWeight;
  std::vector<float> LHEWeight_PDFVariationUpDn;
  std::vector<float> LHEWeight_AsMZUpDn;

};


#endif
