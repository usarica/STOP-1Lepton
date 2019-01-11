#ifndef JECSCALEFACTORHANDLER_H
#define JECSCALEFACTORHANDLER_H


#include <vector>
#include <utility>
#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "Samples.h"
#include "ScaleFactorHandlerBase.h"
#include "JECJERHelpers.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"
#include <CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h>
#include <CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h>


class JECScaleFactorHandler : public ScaleFactorHandlerBase{
public:
  JECJERHelpers::JECJERType const type;

protected:
  FactorizedJetCorrector* corrector_data;
  FactorizedJetCorrector* corrector_MC_noFS;
  FactorizedJetCorrector* corrector_MC_FS;
  JetCorrectionUncertainty* uncertaintyEstimator_MC_noFS;
  JetCorrectionUncertainty* uncertaintyEstimator_MC_FS;

  FactorizedJetCorrector* makeCorrector(std::vector<TString> const& fnames);
  JetCorrectionUncertainty* makeUncertaintyEstimator(TString const& fname);

public:
  JECScaleFactorHandler(JECJERHelpers::JECJERType type_);
  ~JECScaleFactorHandler();

  bool setup();
  void reset();

  void applyJEC(AK4JetObject* obj, bool isMC, bool isFastSim);
  void applyJEC(AK8JetObject* obj, bool isMC, bool isFastSim);

};


#endif
