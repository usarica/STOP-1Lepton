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
  JME::JetResolution resolution_pt;
  JME::JetResolution resolution_phi;
  JME::JetResolutionScaleFactor resolution_sf;

public:
  JECJERHelpers::JECJERType const type;

  JERScaleFactorHandler(JECJERHelpers::JECJERType type_);
  ~JERScaleFactorHandler();

  bool setup();
  void reset();

  void smear(std::vector<std::pair<AK4JetObject*, CMSLorentzVector*>>& jet_genjet_pairs, float rho);
  void smear(std::vector<std::pair<AK8JetObject*, CMSLorentzVector*>>& jet_genjet_pairs, float rho);

};



#endif
