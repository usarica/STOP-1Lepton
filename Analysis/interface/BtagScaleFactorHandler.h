#ifndef BTAGSCALEFACTORHANDLER_H
#define BTAGSCALEFACTORHANDLER_H

// Implementation is a mixture of
// CJLST/ZZAnalysis/AnalysisStep/plugins/JetFiller.cc
// and
// cmstas/StopAnalysis/StopBabyMaker/JetTree.cc

#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "Samples.h"
#include "ScaleFactorHandlerBase.h"
#include "BtagHelpers.h"
#include <cmstas/CORE/Tools/btagsf/BTagCalibrationStandalone.h>


class BtagScaleFactorHandler : public ScaleFactorHandlerBase{
public:
  BtagHelpers::BtagWPType const type;

protected:
  bool isFastSim;
  // related to scale factors
  BTagCalibration* m_calib;
  std::vector<BTagCalibrationReader*> m_readers; // Loads all of [central, up, down](b, c, udsg)

  // related to b tag efficiency
  TFile* m_fileEff;
  TH1* m_hEff[3]; // [0: b, 1: c 2: udsg]

public:
  float WPval;

  BtagScaleFactorHandler(BtagHelpers::BtagWPType type_, bool isFastSim_);
  ~BtagScaleFactorHandler();

  bool setup();
  void reset();

  float getSF(int syst, int jetFlavor, float pt, float eta);
  float getEff(int jetFlavor, float pt, float eta);

};



#endif
