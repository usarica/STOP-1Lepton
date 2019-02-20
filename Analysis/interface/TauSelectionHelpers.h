#ifndef TAUSELECTIONHELPERS_H
#define TAUSELECTIONHELPERS_H

#include "TauObject.h"
#include "TString.h"


namespace TauSelectionHelpers{
  // Taken from StopBabyMaker/looper.cc
  constexpr float ptThr_skim_veto = 20.;
  constexpr float etaThr_skim_veto = 2.3; // Based on tracker acceptance
  constexpr float deltaR_veto_comparison = 0.4;
  const TString strPFIsoVar = "byMediumIsolationMVArun2v1DBoldDMwLT"; // From StopSelections.cc/isVetoTau_v2

  // Veto, loose, medium, tight etc. selection bits
  enum SelectionBits{
    kVetoID,
    kVetoIDIso,
    kSkimPtEta,

    nSelectionBits
  };
  const SelectionBits bit_preselection = kVetoIDIso;

  bool testVetoId(TauObject const& part);
  bool testVetoSelection(TauObject const& part);
  bool testPtEtaSkim(TauObject const& part);

  void setSelectionBits(TauObject& part);

}


#endif