#ifndef AK8JETSELECTIONHELPERS_H
#define AK8JETSELECTIONHELPERS_H


#include "AK8JetObject.h"
#include "METObject.h"


namespace AK8JetSelectionHelpers{
  // Taken from StopBabyMaker/runBabyMaker.cc
  constexpr float ptThr_skim_preselection = 200.;
  constexpr float ptThr_analysis = ptThr_skim_preselection;
  constexpr float ptThr_analysis_btagged = ptThr_skim_preselection;
  constexpr float etaThr_skim_preselection = 2.4;

  // Veto, loose, medium, tight etc. selection bits
  enum SelectionBits{
    kSkimPtEta,
    kPreselection,

    kSkimPtEta_JECUp,
    kPreselection_JECUp,

    kSkimPtEta_JECDn,
    kPreselection_JECDn,

    kSkimPtEta_JERUp,
    kPreselection_JERUp,

    kSkimPtEta_JERDn,
    kPreselection_JERDn
  };

  bool testSkimPtEta(AK8JetObject const& part, int icorr);
  bool testPreselection(AK8JetObject const& part, int icorr);

  void setSelectionBits(AK8JetObject& part);

}


#endif