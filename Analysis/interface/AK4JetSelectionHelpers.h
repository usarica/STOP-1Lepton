#ifndef AK4JETSELECTIONHELPERS_H
#define AK4JETSELECTIONHELPERS_H


#include "AK4JetObject.h"
#include "METObject.h"
#include "BtagHelpers.h"


namespace AK4JetSelectionHelpers{
  // Taken from StopBabyMaker/runBabyMaker.cc
  constexpr float ptThr_skim_preselection = 30.;
  constexpr float ptThr_analysis = 30.;
  constexpr float ptThr_analysis_btagged = 30.;
  constexpr float etaThr_skim_preselection = 2.4;
  constexpr BtagHelpers::BtagWPType AK4Jets_BTagWPType = BtagHelpers::kDeepCSV_Medium;

  // Veto, loose, medium, tight etc. selection bits
  enum SelectionBits{
    kLooseID,
    kTightID,
    kNotBadMuonJet,

    kSkimPtEta,
    kPreselection,

    kSkimPtEta_JECUp,
    kPreselection_JECUp,

    kSkimPtEta_JECDn,
    kPreselection_JECDn,

    kSkimPtEta_JERUp,
    kPreselection_JERUp,

    kSkimPtEta_JERDn,
    kPreselection_JERDn,

    kIsBTagged,
    kIsBTagged_SFUp,
    kIsBTagged_SFDn
  };

  bool isLooseAK4JetPOG(AK4JetObject const& part);
  bool isTightAK4JetPOG(AK4JetObject const& part);

  bool isBadMuonJet(AK4JetObject& part, METObject const& obj);

  bool testLooseId(AK4JetObject const& part);
  bool testLooseSelection(AK4JetObject const& part);

  bool testTightId(AK4JetObject const& part);
  bool testTightSelection(AK4JetObject const& part);

  bool testSkimPtEta(AK4JetObject const& part, int icorr);
  bool testPreselection(AK4JetObject const& part, int icorr);

  void setSelectionBits(AK4JetObject& part); // Does not set the btagging bits!

}


#endif