#ifndef AK4JETSELECTIONHELPERS_H
#define AK4JETSELECTIONHELPERS_H


#include "AK4JetObject.h"


namespace AK4JetSelectionHelpers{
  // Taken from StopBabyMaker/runBabyMaker.cc
  constexpr float ptThr_skim_preselection = 30.;
  constexpr float ptThr_analysis = 30.;
  constexpr float ptThr_analysis_btagged = 30.;
  constexpr float etaThr_skim_preselection = 2.4;

  // Veto, loose, medium, tight etc. selection bits
  enum SelectionBits{
    kLooseID,
    kTightID,
    kSkimPtEta,
    kPreselection
  };

  bool isLooseAK4JetPOG(AK4JetObject const& part);
  bool isTightAK4JetPOG(AK4JetObject const& part);

  bool testLooseId(AK4JetObject const& part);
  bool testLooseSelection(AK4JetObject const& part);

  bool testTightId(AK4JetObject const& part);
  bool testTightSelection(AK4JetObject const& part);

  void setSelectionBits(AK4JetObject& part);

}


#endif