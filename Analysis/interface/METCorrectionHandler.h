#ifndef METCORRECTIONHANDLER_H
#define METCORRECTIONHANDLER_H

#include <vector>
#include <unordered_map>
#include <utility>
#include "Samples.h"
#include "HelperFunctions.h"
#include "METObject.h"
#include "SystematicVariations.h"
#include "ScaleFactorHandlerBase.h"


class METCorrectionHandler : public ScaleFactorHandlerBase{
protected:
  bool applyCorrection;
  std::vector<float> lumilist;
  std::vector<std::vector<std::pair<float, float>>> values_data_map;
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, std::vector<std::vector<std::pair<float, float>>>> values_MC_map;

  void readFile(TString const& strinput);

public:
  METCorrectionHandler();
  ~METCorrectionHandler();

  bool setup();
  void reset();

  bool hasMETCorrection() const{ return applyCorrection; }
  void correctMET(float const& genMET, float const& genMETPhi, METObject* obj, bool useFastSim) const;

  void printParameters() const;

};



#endif
