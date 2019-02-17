#include "ElectronScaleFactorHandler.h"
#include "ElectronSelectionHelpers.h"
#include "TDirectory.h"


using namespace std;
using namespace SampleHelpers;
using namespace ElectronSelectionHelpers;


ElectronScaleFactorHandler::ElectronScaleFactorHandler() :
  ScaleFactorHandlerBase(),
  finput_SF(nullptr),
  finput_SF_tracking(nullptr),
  finput_SF_tracking_lowpt(nullptr),
  finput_SF_veto_eff(nullptr),
  finput_SF_FastSim_id(nullptr),
  finput_SF_FastSim_iso(nullptr),
  finput_SF_FastSim_veto_id(nullptr),
  finput_SF_FastSim_veto_iso(nullptr),
  h_SF_id(nullptr),
  h_SF_iso(nullptr),
  h_SF_tracking(nullptr),
  h_SF_tracking_lowpt(nullptr),
  h_SF_veto_id(nullptr),
  h_SF_veto_iso(nullptr),
  h_SF_FastSim_id(nullptr),
  h_SF_FastSim_iso(nullptr),
  h_SF_FastSim_veto_id(nullptr),
  h_SF_FastSim_veto_iso(nullptr),
  h_SF_veto_eff(nullptr)
{
  setup();
}

ElectronScaleFactorHandler::~ElectronScaleFactorHandler(){ this->reset(); }


bool ElectronScaleFactorHandler::setup(){
  bool res = true;
  TDirectory* curdir = gDirectory;

  this->reset();

  // Recipe: https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF
  // More info: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations
  if (theDataYear == 2016){
    // ID/Iso. and tracking SF files
    finput_SF = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2016/ElectronScaleFactors_Run2016.root", "read");
    finput_SF_tracking = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2016/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root", "read");
    finput_SF_tracking_lowpt = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2016/EGM2D_BtoH_low_RecoSF_Legacy2016.root", "read");
    // FIXME: Only exists for 2016 at the moment
    finput_SF_veto_eff = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2016/lepeff__moriond17__ttbar_powheg_pythia8_25ns.root", "read");
    // Fastsim/Fullsim SF files
    finput_SF_FastSim_id = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2016/sf_el_mediumCB.root", "read");
    finput_SF_FastSim_iso = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2016/sf_el_mini01.root", "read");
    finput_SF_FastSim_veto_id = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2016/sf_el_vetoCB.root", "read");
    finput_SF_FastSim_veto_iso = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2016/sf_el_mini02.root", "read");

    res = (
      // FullSim SFs
      getHistogram(h_SF_id, finput_SF, "Run2016_CutBasedMediumNoIso94XV2")
      && getHistogram(h_SF_iso, finput_SF, "Run2016_Mini")
      && getHistogram(h_SF_tracking, finput_SF_tracking, "EGamma_SF2D")
      && getHistogram(h_SF_tracking_lowpt, finput_SF_tracking_lowpt, "EGamma_SF2D")
      && getHistogram(h_SF_veto_eff, finput_SF_veto_eff, "h2_lepEff_vetoSel_Eff_el")
      && getHistogram(h_SF_veto_id, finput_SF, "Run2016_CutBasedVetoNoIso94XV2")
      && getHistogram(h_SF_veto_iso, finput_SF, "Run2016_Mini2")
      // FastSim/FullSim SFs
      && getHistogram(h_SF_FastSim_id, finput_SF_FastSim_id, "histo2D")
      && getHistogram(h_SF_FastSim_iso, finput_SF_FastSim_iso, "histo2D")
      && getHistogram(h_SF_FastSim_veto_id, finput_SF_FastSim_veto_id, "histo2D")
      && getHistogram(h_SF_FastSim_veto_iso, finput_SF_FastSim_veto_iso, "histo2D")
      );
  }
  else if (theDataYear == 2017){
    // ID/Iso. and tracking SF files
    finput_SF = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2017/ElectronScaleFactors_Run2017.root", "read");
    finput_SF_tracking = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root", "read");
    finput_SF_tracking_lowpt = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root", "read");
    finput_SF_veto_eff = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2017/lepeffs_tt2l_madgraph_mc2017.root", "read");
    // Fastsim/Fullsim SF files
    // FIXME: Update for 2017: No fastsim sample for 2017, use 2016 files at the moment
    finput_SF_FastSim_id = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2017/sf_el_mediumCB.root", "read");
    finput_SF_FastSim_iso = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2017/sf_el_mini01.root", "read");
    finput_SF_FastSim_veto_id = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2017/sf_el_vetoCB.root", "read");
    finput_SF_FastSim_veto_iso = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2017/sf_el_mini02.root", "read");

    res = (
      // FullSim SFs
      getHistogram(h_SF_id, finput_SF, "Run2017_CutBasedMediumNoIso94XV2")
      && getHistogram(h_SF_iso, finput_SF, "Run2017_MVAVLooseTightIP2DMini")
      && getHistogram(h_SF_tracking, finput_SF_tracking, "EGamma_SF2D")
      && getHistogram(h_SF_tracking_lowpt, finput_SF_tracking_lowpt, "EGamma_SF2D")
      && getHistogram(h_SF_veto_eff, finput_SF_veto_eff, "heff17_lepeff_veto_el")
      && getHistogram(h_SF_veto_id, finput_SF, "Run2017_CutBasedVetoNoIso94XV2")
      && getHistogram(h_SF_veto_iso, finput_SF, "Run2017_MVAVLooseTightIP2DMini2")
      // FastSim/FullSim SFs
      && getHistogram(h_SF_FastSim_id, finput_SF_FastSim_id, "histo2D")
      && getHistogram(h_SF_FastSim_iso, finput_SF_FastSim_iso, "histo2D")
      && getHistogram(h_SF_FastSim_veto_id, finput_SF_FastSim_veto_id, "histo2D")
      && getHistogram(h_SF_FastSim_veto_iso, finput_SF_FastSim_veto_iso, "histo2D")
      );
  }
  else if (theDataYear == 2018){
    // ID/Iso. and tracking SF files
    finput_SF = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2018/ElectronScaleFactors_Run2018.root", "read");
    finput_SF_tracking = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2018/egammaEffi.txt_EGM2D.root", "read");
    finput_SF_tracking_lowpt = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2018/egammaEffi.txt_EGM2D_low.root", "read");
    finput_SF_veto_eff = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2018/lepeffs_tt2l_madgraph_mc2017.root", "read"); // FIXME: Needs to be updated
    // Fastsim/Fullsim SF files
    // FIXME: Update for 2018: No fastsim sample for 2018, use 2016 files at the moment
    finput_SF_FastSim_id = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2018/sf_el_mediumCB.root", "read");
    finput_SF_FastSim_iso = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2018/sf_el_mini01.root", "read");
    finput_SF_FastSim_veto_id = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2018/sf_el_vetoCB.root", "read");
    finput_SF_FastSim_veto_iso = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/2018/sf_el_mini02.root", "read");

    res = (
      // FullSim SFs
      getHistogram(h_SF_id, finput_SF, "Run2018_CutBasedMediumNoIso94XV2")
      && getHistogram(h_SF_iso, finput_SF, "Run2018_Mini")
      && getHistogram(h_SF_tracking, finput_SF_tracking, "EGamma_SF2D")
      && getHistogram(h_SF_tracking_lowpt, finput_SF_tracking_lowpt, "EGamma_SF2D")
      && getHistogram(h_SF_veto_eff, finput_SF_veto_eff, "heff17_lepeff_veto_el")
      && getHistogram(h_SF_veto_id, finput_SF, "Run2018_CutBasedVetoNoIso94XV2")
      && getHistogram(h_SF_veto_iso, finput_SF, "Run2018_Mini2")
      // FastSim/FullSim SFs
      && getHistogram(h_SF_FastSim_id, finput_SF_FastSim_id, "histo2D")
      && getHistogram(h_SF_FastSim_iso, finput_SF_FastSim_iso, "histo2D")
      && getHistogram(h_SF_FastSim_veto_id, finput_SF_FastSim_veto_id, "histo2D")
      && getHistogram(h_SF_FastSim_veto_iso, finput_SF_FastSim_veto_iso, "histo2D")
      );
  }

  curdir->cd();

  return res;
}
void ElectronScaleFactorHandler::reset(){
  ScaleFactorHandlerBase::closeFile(finput_SF);
  ScaleFactorHandlerBase::closeFile(finput_SF_tracking);
  ScaleFactorHandlerBase::closeFile(finput_SF_tracking_lowpt);
  ScaleFactorHandlerBase::closeFile(finput_SF_veto_eff);
  ScaleFactorHandlerBase::closeFile(finput_SF_FastSim_id);
  ScaleFactorHandlerBase::closeFile(finput_SF_FastSim_iso);
  ScaleFactorHandlerBase::closeFile(finput_SF_FastSim_veto_id);
  ScaleFactorHandlerBase::closeFile(finput_SF_FastSim_veto_iso);

  h_SF_id = nullptr;
  h_SF_iso = nullptr;
  h_SF_tracking = nullptr;
  h_SF_tracking_lowpt = nullptr;
  h_SF_veto_id = nullptr;
  h_SF_veto_iso = nullptr;
  h_SF_FastSim_id = nullptr;
  h_SF_FastSim_iso = nullptr;
  h_SF_FastSim_veto_id = nullptr;
  h_SF_FastSim_veto_iso = nullptr;
  h_SF_veto_eff = nullptr;
}

void ElectronScaleFactorHandler::evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, ElectronObject const* obj, TH2F const* hist, bool etaOnY, bool useAbsEta) const{
  if (!hist) return;
  if (!obj) return;

  float const ele_pt = obj->pt();
  float const& ele_etasc = obj->extras.etaSC;

  int ix, iy;
  int nbinsx = hist->GetNbinsX();
  int nbinsy = hist->GetNbinsY();
  if (!etaOnY){
    ix = hist->GetXaxis()->FindBin((!useAbsEta ? ele_etasc : fabs(ele_etasc)));
    iy = hist->GetYaxis()->FindBin(ele_pt);
  }
  else{
    ix = hist->GetXaxis()->FindBin(ele_pt);
    iy = hist->GetYaxis()->FindBin((!useAbsEta ? ele_etasc : fabs(ele_etasc)));
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

void ElectronScaleFactorHandler::getIdIsoSFAndError(float& theSF, float& theSFRelErr, ElectronObject const* obj, bool useFastSim) const{
  theSF=1; theSFRelErr=0;

  if (!obj) return;
  bool passSel = obj->testSelection(bit_preselection_idiso);
  bool passVeto = obj->testSelection(kVetoID) && !passSel;
  if (!passSel && !passVeto) return;

  if (!useFastSim){
    if (passSel){
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_id, false, false);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_iso, false, false);
    }
    else{
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_veto_id, false, false);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_veto_iso, false, false);
    }
  }
  else{
    // FIXME: Need to check axis inversion after 2017 hists are obtained
    if (passSel){
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_FastSim_id, true, true);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_FastSim_iso, true, true);
    }
    else{
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_FastSim_veto_id, true, true);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_FastSim_veto_iso, true, true);
    }
  }

  // Uncertainty in veto id+iso eff. is doubled as pT thr. is lower in this analysis
  if (passVeto && !useFastSim) theSFRelErr *= 2.f;
}
void ElectronScaleFactorHandler::getRecoSFAndError(float& theSF, float& theSFRelErr, ElectronObject const* obj) const{
  theSF=1; theSFRelErr=0;

  if (!obj) return;
  bool passSel = obj->testSelection(bit_preselection_idisoreco);
  //bool passVeto = obj->testSelection(kVetoIDReco) && !passSel;
  //if (!passSel && !passVeto) return;
  if (!passSel) return;

  evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, (obj->pt()>=20.f ? h_SF_tracking : h_SF_tracking_lowpt), false, false);
}

void ElectronScaleFactorHandler::getGenSFAndError(float& theSF, float& theSFRelErr, ElectronObject const* obj, float const& theIdIsoSF, float const& theIdIsoSFRelErr) const{
  theSF=1; theSFRelErr=0;

  if (!obj) return;
  if (!obj->testSelection(kGenPtEta)) return;

  // FIXME: Histogram for 2017 might not use abs(eta) or be swapped.
  evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_veto_eff, true, true);

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

