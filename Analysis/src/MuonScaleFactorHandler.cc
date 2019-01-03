#include "MuonScaleFactorHandler.h"
#include "MuonSelectionHelpers.h"
#include "TDirectory.h"


using namespace std;
using namespace SampleHelpers;
using namespace MuonSelectionHelpers;


MuonScaleFactorHandler::MuonScaleFactorHandler() :
  ScaleFactorHandlerBase(),
  finput_SF_id(nullptr),
  finput_SF_iso(nullptr),
  finput_SF_ip(nullptr),
  finput_SF_tracking(nullptr),
  finput_SF_veto_eff(nullptr),
  finput_SF_veto_id(nullptr),
  finput_SF_veto_iso(nullptr),
  finput_SF_veto_ip(nullptr),
  finput_SF_FastSim_id(nullptr),
  finput_SF_FastSim_iso(nullptr),
  finput_SF_FastSim_ip(nullptr),
  finput_SF_FastSim_veto_id(nullptr),
  finput_SF_FastSim_veto_iso(nullptr),
  finput_SF_FastSim_veto_ip(nullptr),
  h_SF_id(nullptr),
  h_SF_iso(nullptr),
  h_SF_ip(nullptr),
  h_SF_tracking(nullptr),
  h_SF_veto_id(nullptr),
  h_SF_veto_iso(nullptr),
  h_SF_veto_ip(nullptr),
  h_SF_FastSim_id(nullptr),
  h_SF_FastSim_iso(nullptr),
  h_SF_FastSim_ip(nullptr),
  h_SF_FastSim_veto_id(nullptr),
  h_SF_FastSim_veto_iso(nullptr),
  h_SF_FastSim_veto_ip(nullptr),
  h_SF_veto_eff(nullptr)
{
  setup();
}

MuonScaleFactorHandler::~MuonScaleFactorHandler(){ this->reset(); }


bool MuonScaleFactorHandler::setup(){
  bool res = true;
  this->reset();

  // Recipe: https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF
  // More info: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations
  if (theDataYear == 2016){
    // ID/Iso./IP and tracking SF files
    finput_SF_id = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root", "read");
    finput_SF_iso = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/TnP_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root", "read");
    finput_SF_ip = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root", "read");
    //finput_SF_tracking = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/Tracking_EfficienciesAndSF_BCDEFGH.root", "read");
    // Veto SFs
    finput_SF_veto_eff = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/lepeff__moriond17__ttbar_powheg_pythia8_25ns.root", "read");
    finput_SF_veto_id = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/TnP_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root", "read");
    finput_SF_veto_iso = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/TnP_NUM_MiniIsoTight_DENOM_LooseID_VAR_map_pt_eta.root", "read");
    finput_SF_veto_ip = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/TnP_NUM_MediumIP2D_DENOM_LooseID_VAR_map_pt_eta.root", "read");
    // Fastsim/Fullsim SF files
    finput_SF_FastSim_id = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/sf_mu_mediumID.root", "read");
    finput_SF_FastSim_iso = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/sf_mu_mediumID_mini02.root", "read");
    finput_SF_FastSim_ip = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/sf_mu_mediumID_tightIP2D.root", "read");
    finput_SF_FastSim_veto_id = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/sf_mu_looseID.root", "read");
    finput_SF_FastSim_veto_iso = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/sf_mu_looseID_mini02.root", "read");
    finput_SF_FastSim_veto_ip = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2016/sf_mu_mediumID_looseIP2D.root", "read");

    res = (
      // FullSim SFs
      getHistogram(h_SF_id, finput_SF_id, "SF")
      && getHistogram(h_SF_iso, finput_SF_iso, "SF")
      && getHistogram(h_SF_ip, finput_SF_ip, "SF")
      //&& getHistogram(h_SF_tracking, finput_SF_tracking, hname_SF_tracking)
      && getHistogram(h_SF_veto_eff, finput_SF_veto_eff, "h2_lepEff_vetoSel_Eff_mu")
      && getHistogram(h_SF_veto_id, finput_SF_veto_id, "SF")
      && getHistogram(h_SF_veto_iso, finput_SF_veto_iso, "SF")
      && getHistogram(h_SF_veto_ip, finput_SF_veto_ip, "SF")
      // FastSim/FullSim SFs
      && getHistogram(h_SF_FastSim_id, finput_SF_FastSim_id, "histo2D")
      && getHistogram(h_SF_FastSim_iso, finput_SF_FastSim_iso, "histo2D")
      && getHistogram(h_SF_FastSim_ip, finput_SF_FastSim_ip, "histo2D")
      && getHistogram(h_SF_FastSim_veto_id, finput_SF_FastSim_veto_id, "histo2D")
      && getHistogram(h_SF_FastSim_veto_iso, finput_SF_FastSim_veto_iso, "histo2D")
      && getHistogram(h_SF_FastSim_veto_ip, finput_SF_FastSim_veto_ip, "histo2D")
      );
  }
  else if (theDataYear == 2017 || theDataYear == 2018){
    // ID/Iso./IP and tracking SF files
    finput_SF_id = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/RunBCDEF_SF_ID.root", "read");
    finput_SF_iso = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/sf_mu_iso_susy_2017.root", "read");
    //finput_SF_ip = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root", "read");
    //finput_SF_tracking = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/Tracking_EfficienciesAndSF_BCDEFGH.root", "read");
    // Veto SFs
    finput_SF_veto_eff = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/lepeffs_tt2l_madgraph_mc2017.root", "read");
    //finput_SF_veto_id = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/RunBCDEF_SF_ID.root", "read");
    //finput_SF_veto_iso = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/sf_mu_iso_susy_2017.root", "read");
    //finput_SF_veto_ip = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/TnP_NUM_MediumIP2D_DENOM_LooseID_VAR_map_pt_eta.root", "read");
    // Fastsim/Fullsim SF files
    finput_SF_FastSim_id = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/sf_mu_mediumID.root", "read");
    finput_SF_FastSim_iso = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/sf_mu_mediumID_mini02.root", "read");
    finput_SF_FastSim_ip = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/sf_mu_mediumID_tightIP2D.root", "read");
    finput_SF_FastSim_veto_id = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/sf_mu_looseID.root", "read");
    finput_SF_FastSim_veto_iso = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/sf_mu_looseID_mini02.root", "read");
    finput_SF_FastSim_veto_ip = TFile::Open(STOP1LPKGDATAPATH+"MuSFs/2017/sf_mu_mediumID_looseIP2D.root", "read");

    res = (
      // FullSim SFs
      getHistogram(h_SF_id, finput_SF_id, "NUM_MediumPromptID_DEN_genTracks_pt_abseta")
      && getHistogram(h_SF_iso, finput_SF_iso, "TnP_MC_NUM_MiniIso02Cut_DEN_MediumID_PAR_pt_eta")
      //&& getHistogram(h_SF_ip, finput_SF_ip, "SF")
      //&& getHistogram(h_SF_tracking, finput_SF_tracking, hname_SF_tracking)
      && getHistogram(h_SF_veto_eff, finput_SF_veto_eff, "heff17_lepeff_veto_mu")
      && getHistogram(h_SF_veto_id, finput_SF_id/*finput_SF_veto_id*/, "NUM_LooseID_DEN_genTracks_pt_abseta")
      && getHistogram(h_SF_veto_iso, finput_SF_iso/*finput_SF_veto_iso*/, "TnP_MC_NUM_MiniIso02Cut_DEN_LooseID_PAR_pt_eta")
      //&& getHistogram(h_SF_veto_ip, finput_SF_veto_ip, "SF")
      // FastSim/FullSim SFs
      && getHistogram(h_SF_FastSim_id, finput_SF_FastSim_id, "histo2D")
      && getHistogram(h_SF_FastSim_iso, finput_SF_FastSim_iso, "histo2D")
      && getHistogram(h_SF_FastSim_ip, finput_SF_FastSim_ip, "histo2D")
      && getHistogram(h_SF_FastSim_veto_id, finput_SF_FastSim_veto_id, "histo2D")
      && getHistogram(h_SF_FastSim_veto_iso, finput_SF_FastSim_veto_iso, "histo2D")
      && getHistogram(h_SF_FastSim_veto_ip, finput_SF_FastSim_veto_ip, "histo2D")
      );
  }
  //else if (theDataPeriod == "2018"){
  // FIXME: To be implemented
  //}

  return res;
}
void MuonScaleFactorHandler::reset(){
  ScaleFactorHandlerBase::closeFile(finput_SF_id);
  ScaleFactorHandlerBase::closeFile(finput_SF_iso);
  ScaleFactorHandlerBase::closeFile(finput_SF_ip);
  ScaleFactorHandlerBase::closeFile(finput_SF_tracking);
  ScaleFactorHandlerBase::closeFile(finput_SF_veto_eff);
  ScaleFactorHandlerBase::closeFile(finput_SF_FastSim_id);
  ScaleFactorHandlerBase::closeFile(finput_SF_FastSim_iso);
  ScaleFactorHandlerBase::closeFile(finput_SF_FastSim_ip);
  ScaleFactorHandlerBase::closeFile(finput_SF_FastSim_veto_id);
  ScaleFactorHandlerBase::closeFile(finput_SF_FastSim_veto_iso);
  ScaleFactorHandlerBase::closeFile(finput_SF_FastSim_veto_ip);

  h_SF_id = nullptr;
  h_SF_iso = nullptr;
  h_SF_ip = nullptr;
  h_SF_tracking = nullptr;
  h_SF_veto_id = nullptr;
  h_SF_veto_iso = nullptr;
  h_SF_veto_ip = nullptr;
  h_SF_FastSim_id = nullptr;
  h_SF_FastSim_iso = nullptr;
  h_SF_FastSim_ip = nullptr;
  h_SF_FastSim_veto_id = nullptr;
  h_SF_FastSim_veto_iso = nullptr;
  h_SF_FastSim_veto_ip = nullptr;
  h_SF_veto_eff = nullptr;
}

void MuonScaleFactorHandler::evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, MuonObject const* obj, TH2F const* hist, bool etaOnY, bool useAbsEta) const{
  if (!hist) return;
  if (!obj) return;

  float const pt = obj->pt();
  float const eta = obj->eta();

  int ix, iy;
  int nbinsx = hist->GetNbinsX();
  int nbinsy = hist->GetNbinsY();
  if (!etaOnY){
    ix = hist->GetXaxis()->FindBin((!useAbsEta ? eta : fabs(eta)));
    iy = hist->GetYaxis()->FindBin(pt);
  }
  else{
    ix = hist->GetXaxis()->FindBin(pt);
    iy = hist->GetYaxis()->FindBin((!useAbsEta ? eta : fabs(eta)));
  }
  if (ix==0) ix=1;
  else if (ix==nbinsx+1) ix=nbinsx;
  if (iy==0) iy=1;
  else if (iy==nbinsy+1) iy=nbinsy;

  float bc = hist->GetBinContent(ix, iy);
  float be = hist->GetBinError(ix, iy);
  if (bc!=0.f) be /= bc;
  if (be<0.f) be=0;

  theSF *= bc; theSFRelErr = sqrt(pow(theSFRelErr, 2)+pow(be, 2));
}

void MuonScaleFactorHandler::getIdIsoSFAndError(float& theSF, float& theSFRelErr, MuonObject const* obj, bool useFastSim) const{
  theSF=1; theSFRelErr=0;

  if (!obj) return;
  bool passSel = obj->testSelection(bit_preselection_idiso);
  bool passVeto = obj->testSelection(kVetoID) && !passSel;
  if (!passSel && !passVeto) return;

  if (!useFastSim){
    if (passSel){
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_iso, true, true);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_ip, true, true);
      theSFRelErr *= 2.f; // Double iso and ip errors
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_id, true, true);
    }
    else{
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_veto_iso, true, true);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_veto_ip, true, true);
      theSFRelErr *= 2.f; // Double iso and ip errors
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_veto_id, true, true);
    }
  }
  else{
    if (passSel){
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_FastSim_iso, true, true);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_FastSim_ip, true, true);
      theSFRelErr *= 2.f; // Double iso and ip errors
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_FastSim_id, true, true);
    }
    else{
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_FastSim_veto_iso, true, true);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_FastSim_veto_ip, true, true);
      theSFRelErr *= 2.f; // Double iso and ip errors
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_FastSim_veto_id, true, true);
    }
  }
}
void MuonScaleFactorHandler::getRecoSFAndError(float& theSF, float& theSFRelErr, MuonObject const* obj) const{
  theSF=1; theSFRelErr=0;

  if (!obj) return;
  bool passSel = obj->testSelection(bit_preselection_idisoreco);
  //bool passVeto = obj->testSelection(kVetoIDReco) && !passSel;
  //if (!passSel && !passVeto) return;
  if (!passSel) return;

  evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_tracking, false, false);
}

void MuonScaleFactorHandler::getGenSFAndError(float& theSF, float& theSFRelErr, MuonObject const* obj, float const& theIdIsoSF, float const& theIdIsoSFRelErr) const{
  theSF=1; theSFRelErr=0;

  if (!obj) return;
  if (!obj->testSelection(kGenPtEta)) return;

  // FIXME: Histogram for 2017 might not use abs(eta) or be swapped.
  evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_veto_eff, true/*(theDataPeriod=="2016")*/, true/*(theDataPeriod=="2016")*/);

  if (theSF<1.){
    const float theSFtmp = theSF;
    theSF = (1.-(theSFtmp*theIdIsoSF))/(1.-theSFtmp);
    if (theSF<0.) theSF=0.;
    theSFRelErr = fabs(theSFtmp*theIdIsoSFRelErr*theIdIsoSF/(1.-theSFtmp));
    if (theSF>0.) theSFRelErr /= theSF;
    else theSFRelErr=0;
  }

  // NOTE: THIS WEIGHT SHOULD NOT BE MULTIPLIED WITH THE RECO (i.e. TRACKING) WEIGHT
}

