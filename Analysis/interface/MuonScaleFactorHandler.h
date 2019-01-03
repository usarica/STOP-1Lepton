#ifndef MUONSCALEFACTORHANDLER_H
#define MUONSCALEFACTORHANDLER_H

#include "TH2F.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "Samples.h"
#include "ScaleFactorHandlerBase.h"
#include "MuonObject.h"


class MuonScaleFactorHandler : public ScaleFactorHandlerBase{
protected:
  TFile* finput_SF_id;
  TFile* finput_SF_iso;
  TFile* finput_SF_ip;
  TFile* finput_SF_tracking;

  TFile* finput_SF_veto_eff;

  TFile* finput_SF_veto_id;
  TFile* finput_SF_veto_iso;
  TFile* finput_SF_veto_ip;

  TFile* finput_SF_FastSim_id;
  TFile* finput_SF_FastSim_iso;
  TFile* finput_SF_FastSim_ip;

  TFile* finput_SF_FastSim_veto_id;
  TFile* finput_SF_FastSim_veto_iso;
  TFile* finput_SF_FastSim_veto_ip;

  // FullSim SFs
  TH2F* h_SF_id;
  TH2F* h_SF_iso;
  TH2F* h_SF_ip;
  TH2F* h_SF_tracking;
  TH2F* h_SF_veto_id;
  TH2F* h_SF_veto_iso;
  TH2F* h_SF_veto_ip;

  // FastSim/FullSim SFs
  TH2F* h_SF_FastSim_id;
  TH2F* h_SF_FastSim_iso;
  TH2F* h_SF_FastSim_ip;
  TH2F* h_SF_FastSim_veto_id;
  TH2F* h_SF_FastSim_veto_iso;
  TH2F* h_SF_FastSim_veto_ip;

  // FastSim/FullSim SFs
  TH2F* h_SF_veto_eff;

  void evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, MuonObject const* obj, TH2F const* hist, bool etaOnY, bool useAbsEta) const;

public:
  MuonScaleFactorHandler();
  ~MuonScaleFactorHandler();

  bool setup();
  void reset();

  void getIdIsoSFAndError(float& theSF, float& theSFRelErr, MuonObject const* obj, bool useFastSim) const;
  void getRecoSFAndError(float& theSF, float& theSFRelErr, MuonObject const* obj) const;
  void getGenSFAndError(float& theSF, float& theSFRelErr, MuonObject const* obj, float const& theIdIsoSF, float const& theIdIsoSFRelErr) const;

};



#endif
