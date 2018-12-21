#ifndef ELECTRONSELECTIONHELPERS_H
#define ELECTRONSELECTIONHELPERS_H

#include "ElectronObject.h"


namespace ElectronSelectionHelpers{
  extern int eleEffAreaVersion;

  // Taken from StopBabyMaker/runBabyMaker.cc
  constexpr float ptThr_skim_veto = 5.;
  constexpr float ptThr_skim_loose = 10.;
  constexpr float ptThr_skim_good = 20.;
  constexpr float etaThr_skim_veto = 2.4;
  constexpr float etaThr_skim_loose = 2.4;
  constexpr float etaThr_skim_good = 1.4442;

  int setEleEffAreaVersion();
  float electronEffArea_DR0p3(ElectronObject const& part);

  bool testVetoCutBasedId(ElectronObject const& part);
  bool testVetoSelection(ElectronObject const& part);

}


#endif