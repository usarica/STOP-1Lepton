#ifndef JERSCALEFACTORHANDLER_H
#define JERSCALEFACTORHANDLER_H


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
#include <JetMETCorrections/Modules/interface/JetResolution.h>


class JERScaleFactorHandler : public ScaleFactorHandlerBase{
protected:
  JME::JetResolution resolution_pt_data;
  JME::JetResolution resolution_phi_data;
  JME::JetResolution resolution_pt_mc;
  JME::JetResolution resolution_phi_mc;
  JME::JetResolutionScaleFactor resolution_sf; // Only MC

public:
  JECJERHelpers::JECJERType const type;

  JERScaleFactorHandler(JECJERHelpers::JECJERType type_);
  ~JERScaleFactorHandler();

  bool setup();
  void reset();

  void smear(std::vector<AK4JetObject*>& jets, float rho, bool isMC);
  void smear(std::vector<AK8JetObject*>& jets, float rho, bool isMC);

};



#endif
