#ifndef ELECTRONSCALEFACTORHANDLER_H
#define ELECTRONSCALEFACTORHANDLER_H

#include "TH2F.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "Samples.h"
#include "ScaleFactorHandlerBase.h"


class ElectronScaleFactorHandler : public ScaleFactorHandlerBase{
protected:
  TFile* finput_SF;
  TFile* finput_SF_tracking;
  TFile* finput_SF_veto_eff;

  TFile* finput_SF_FastSim_id;
  TFile* finput_SF_FastSim_iso;
  TFile* finput_SF_FastSim_veto_id;
  TFile* finput_SF_FastSim_veto_iso;

  // FullSim SFs
  TH2F* h_SF_id;
  TH2F* h_SF_iso;
  TH2F* h_SF_tracking;
  TH2F* h_SF_veto_id;
  TH2F* h_SF_veto_iso;

  // FastSim/FullSim SFs
  TH2F* h_SF_FastSim_id;
  TH2F* h_SF_FastSim_iso;
  TH2F* h_SF_FastSim_veto_id;
  TH2F* h_SF_FastSim_veto_iso;

  // FastSim/FullSim SFs
  TH2F* h_SF_veto_eff;

  void evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, float const& ele_pt, float const& ele_etasc, TH2F const* hist, bool etaOnY, bool useAbsEta) const;

public:
  ElectronScaleFactorHandler();
  ~ElectronScaleFactorHandler();

  void getIdIsoSFAndError(float& theSF, float& theSFRelErr, float const& ele_pt, float const& ele_etasc, bool isVeto, bool useFastSim) const;
  void getRecoSFAndError(float& theSF, float& theSFRelErr, float const& ele_pt, float const& ele_etasc) const;
  void getGenSFAndError(float& theSF, float& theSFRelErr, float const& ele_pt, float const& ele_eta, float const& theIdIsoSF, float const& theIdIsoSFRelErr) const;

};



#endif
