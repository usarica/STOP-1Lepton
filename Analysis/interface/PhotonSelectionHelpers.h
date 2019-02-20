#ifndef PHOTONSELECTIONHELPERS_H
#define PHOTONSELECTIONHELPERS_H

#include "PhotonObject.h"


namespace PhotonSelectionHelpers{
  // Taken from StopBabyMaker/runBabyMaker.cc
  constexpr float ptThr_gen = 5.;
  constexpr float ptThr_skim_loose = 25.;
  constexpr float ptThr_skim_medium = 50.;
  constexpr float ptThr_skim_tight = 50.;
  constexpr float etaThr_gen = 2.5;
  constexpr float etaThr_skim_loose = 2.5;
  constexpr float etaThr_skim_medium = 2.5; // 1.4442 last ECAL crystal in barrel
  constexpr float etaThr_skim_tight = 2.5; // 1.4442 last ECAL crystal in barrel

  // Veto, loose, medium, tight etc. selection bits
  enum SelectionBits{
    kGenPtEta,
    kLooseIDReco,
    kMediumIDReco,
    kTightIDReco,
    kSkimPtEta,
    kPreselection,

    nSelectionBits
  };
  const SelectionBits bit_preselection_idisoreco = kTightIDReco;


  float photonCHEffArea_DR0p3(PhotonObject const& part);
  float photonNHEffArea_DR0p3(PhotonObject const& part);
  float photonEMEffArea_DR0p3(PhotonObject const& part);

  float photonCHIsoCorr_DR0p3(PhotonObject const& part);
  float photonNHIsoCorr_DR0p3(PhotonObject const& part);
  float photonEMIsoCorr_DR0p3(PhotonObject const& part);

  bool testPtEtaGen(PhotonObject const& part);

  bool testLooseCutBasedId(PhotonObject const& part);
  bool testLooseSelection(PhotonObject const& part);

  bool testMediumCutBasedId(PhotonObject const& part);
  bool testMediumSelection(PhotonObject const& part);

  bool testTightCutBasedId(PhotonObject const& part);
  bool testTightSelection(PhotonObject const& part);

  bool testPtEtaSkim(PhotonObject const& part);
  bool testPreselection(PhotonObject const& part);

  void setSelectionBits(PhotonObject& part);

}


#endif
