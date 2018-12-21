#ifndef ELECTRONSELECTIONHELPERS_H
#define ELECTRONSELECTIONHELPERS_H

#include "ElectronObject.h"


namespace ElectronSelectionHelpers{
  extern int eleEffAreaVersion;

  // Taken from StopBabyMaker/runBabyMaker.cc
  constexpr float ptThr_skim_veto = 5.;
  constexpr float ptThr_skim_loose = 10.;
  constexpr float ptThr_skim_preselection = 20.;
  constexpr float etaThr_skim_veto = 2.4;
  constexpr float etaThr_skim_loose = 2.4;
  constexpr float etaThr_skim_preselection = 1.4442;

  int setEleEffAreaVersion();
  float electronEffArea_DR0p3(ElectronObject const& part);

  bool testVetoCutBasedId(ElectronObject const& part);
  bool testLooseCutBasedId(ElectronObject const& part);
  bool testMediumCutBasedId(ElectronObject const& part);

  bool testVetoSelection(ElectronObject const& part);
  bool testLooseSelection(ElectronObject const& part);
  bool testMediumSelection(ElectronObject const& part);

}


#endif