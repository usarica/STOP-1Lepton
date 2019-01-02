#ifndef MUONSELECTIONHELPERS_H
#define MUONSELECTIONHELPERS_H


#include "MuonObject.h"


namespace MuonSelectionHelpers{
  // Taken from StopBabyMaker/runBabyMaker.cc
  constexpr float ptThr_gen = 5.;
  constexpr float ptThr_skim_veto = 5.;
  constexpr float ptThr_skim_loose = 10.;
  constexpr float ptThr_skim_medium = 20.;
  constexpr float etaThr_gen = 2.4;
  constexpr float etaThr_skim_veto = 2.4;
  constexpr float etaThr_skim_loose = 2.4;
  constexpr float etaThr_skim_medium = 2.4;

  // To keep updated with: https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonReco/interface/Muon.h#L189-L220
  enum POGSelectorBits{
    CutBasedIdLoose=0,
    CutBasedIdMedium,
    CutBasedIdMediumPrompt, // medium with IP cuts    
    CutBasedIdTight,
    CutBasedIdGlobalHighPt, // high pt muon for Z',W' (better momentum resolution)
    CutBasedIdTrkHighPt, // high pt muon for boosted Z (better efficiency)
    PFIsoVeryLoose, // reliso<0.40       
    PFIsoLoose, // reliso<0.25       
    PFIsoMedium, // reliso<0.20       
    PFIsoTight, // reliso<0.15       
    PFIsoVeryTight, // reliso<0.10       
    TkIsoLoose, // reliso<0.10       
    TkIsoTight, // reliso<0.05       
    SoftCutBasedId,
    SoftMvaId,
    MvaLoose,
    MvaMedium,
    MvaTight,
    MiniIsoLoose, // reliso<0.40       
    MiniIsoMedium, // reliso<0.20       
    MiniIsoTight, // reliso<0.10       
    MiniIsoVeryTight, // reliso<0.05       
    TriggerIdLoose, // robust selector for HLT    
    InTimeMuon,
    PFIsoVeryVeryTight, // reliso<0.05       
    MultiIsoLoose, // miniIso with ptRatio and ptRel   
    MultiIsoMedium // miniIso with ptRatio and ptRel   
  };

  // Veto, loose, medium, tight etc. selection bits
  enum SelectionBits{
    kGenPtEta,
    kVetoID,
    kVetoIDReco,
    kLooseID,
    kLooseIDReco,
    kMediumID,
    kMediumIDReco,
    kTightID,
    kTightIDReco,
    //kTrigger,
    kSkimPtEta,
    kPreselection
  };
  const SelectionBits bit_preselection_idiso = kMediumID;
  const SelectionBits bit_preselection_idisoreco = kMediumIDReco;

  bool checkPOGSelectorBit(MuonObject const& part, POGSelectorBits ibit);

  int setMuonEffAreaVersion();
  float muonEffArea(MuonObject const& part);

  float miniAbsIso_DR0p3(MuonObject const& part);
  float miniRelIso_DR0p3(MuonObject const& part);

  bool isLooseMuonPOG(MuonObject const& part);
  bool isMediumMuonPOG(MuonObject const& part);

  bool testPtEtaGen(MuonObject const& part);

  bool testVetoCutBasedId(MuonObject const& part);
  bool testVetoSelection(MuonObject const& part);

  bool testLooseCutBasedId(MuonObject const& part);
  bool testLooseSelection(MuonObject const& part);

  bool testMediumCutBasedId(MuonObject const& part);
  bool testMediumSelection(MuonObject const& part);

  bool testPtEtaSkim(MuonObject const& part);
  bool testPreselection(MuonObject const& part);

  void setSelectionBits(MuonObject& part);

}


#endif