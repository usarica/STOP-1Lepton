#ifndef ISOTRACKSELECTIONHELPERS_H
#define ISOTRACKSELECTIONHELPERS_H

#include "IsoTrackObject.h"


namespace IsoTrackSelectionHelpers{
  // Taken from StopBabyMaker/looper.cc
  constexpr float ptThr_skim_veto = 10.;
  constexpr float etaThr_skim_veto = 2.4; // Based on tracker acceptance
  constexpr float deltaR_veto_comparison = 0.4;

  // Veto, loose, medium, tight etc. selection bits
  enum SelectionBits{
    kVetoID,
    kVetoIDIso,
    kSkimPtEta,

    nSelectionBits
  };
  const SelectionBits bit_preselection = kVetoIDIso;

  bool testVetoId(IsoTrackObject const& part);
  bool testVetoSelection(IsoTrackObject const& part);
  bool testPtEtaSkim(IsoTrackObject const& part);

  void setSelectionBits(IsoTrackObject& part);

}


#endif