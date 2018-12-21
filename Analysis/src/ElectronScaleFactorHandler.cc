#include "ElectronScaleFactorHandler.h"
#include "TDirectory.h"


using namespace std;
using namespace SampleHelpers;


ElectronScaleFactorHandler::ElectronScaleFactorHandler() :
  finput_SF(nullptr),
  finput_SF_tracking(nullptr),
  finput_SF_veto_eff(nullptr),
  finput_SF_FastSim_id(nullptr),
  finput_SF_FastSim_iso(nullptr),
  finput_SF_FastSim_veto_id(nullptr),
  finput_SF_FastSim_veto_iso(nullptr),
  h_SF_id(nullptr),
  h_SF_iso(nullptr),
  h_SF_tracking(nullptr),
  h_SF_veto_id(nullptr),
  h_SF_veto_iso(nullptr),
  h_SF_FastSim_id(nullptr),
  h_SF_FastSim_iso(nullptr),
  h_SF_FastSim_veto_id(nullptr),
  h_SF_FastSim_veto_iso(nullptr),
  h_SF_veto_eff(nullptr)
{
  TDirectory* curdir = gDirectory;
  // Recipe: https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF
  // More info: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations
  if (theDataPeriod == "2016"){
    // ID/Iso. and tracking SF files
    finput_SF = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/scaleFactors.root", "read");
    finput_SF_tracking = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/egammaEffi.txt_EGM2D.root", "read");
    // FIXME: Only exists for 2016 at the moment
    finput_SF_veto_eff = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/lepeff__moriond17__ttbar_powheg_pythia8_25ns.root", "read");
    // Fastsim/Fullsim SF files
    finput_SF_FastSim_id = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/sf_el_mediumCB.root", "read");
    finput_SF_FastSim_iso = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/sf_el_mini01.root", "read");
    finput_SF_FastSim_veto_id = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/sf_el_vetoCB.root", "read");
    finput_SF_FastSim_veto_iso = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/sf_el_mini02.root", "read");

    // FullSim SFs
    finput_SF->cd();
    h_SF_id = (TH2F*) finput_SF->Get("GsfElectronToCutBasedSpring15M");
    h_SF_iso = (TH2F*) finput_SF->Get("MVAVLooseElectronToMini");
    h_SF_veto_id = (TH2F*) finput_SF->Get("GsfElectronToCutBasedSpring15V");
    h_SF_veto_iso = (TH2F*) finput_SF->Get("MVAVLooseElectronToMini2");
    finput_SF_tracking->cd();
    h_SF_tracking = (TH2F*) finput_SF_tracking->Get("EGamma_SF2D");
    finput_SF_veto_eff->cd();
    h_SF_veto_eff = (TH2F*) finput_SF_veto_eff->Get("h2_lepEff_vetoSel_Eff_el");

    // FastSim/FullSim SFs
    finput_SF_FastSim_id->cd(); h_SF_FastSim_id = (TH2F*) finput_SF_FastSim_id->Get("histo2D");
    finput_SF_FastSim_iso->cd(); h_SF_FastSim_iso = (TH2F*) finput_SF_FastSim_iso->Get("histo2D");
    finput_SF_FastSim_veto_id->cd(); h_SF_FastSim_veto_id = (TH2F*) finput_SF_FastSim_veto_id->Get("histo2D");
    finput_SF_FastSim_veto_iso->cd(); h_SF_FastSim_veto_iso = (TH2F*) finput_SF_FastSim_veto_iso->Get("histo2D");

  }
  else if (theDataPeriod == "2017"){
    // ID/Iso. and tracking SF files
    finput_SF = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/ElectronScaleFactors_Run2017.root", "read");
    finput_SF_tracking = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root", "read");
    // FIXME: Only exists for 2016 at the moment
    //finput_SF_veto_eff = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/lepeff__moriond17__ttbar_powheg_pythia8_25ns.root", "read");
    // Fastsim/Fullsim SF files
    // FIXME: Update for 2017: No fastsim sample for 2017, use 2016 files at the moment
    finput_SF_FastSim_id = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/sf_el_mediumCB.root", "read");
    finput_SF_FastSim_iso = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/sf_el_mini01.root", "read");
    finput_SF_FastSim_veto_id = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/sf_el_vetoCB.root", "read");
    finput_SF_FastSim_veto_iso = TFile::Open(STOP1LPKGDATAPATH+"EleSFs/sf_el_mini02.root", "read");

    // FullSim SFs
    finput_SF->cd();
    h_SF_id = (TH2F*) finput_SF->Get("Run2017_CutBasedMediumNoIso94XV2"); // Cut-based medium ID
    h_SF_iso = (TH2F*) finput_SF->Get("Run2017_MVAVLooseTightIP2DMini"); // MiniIso < 0.1
    h_SF_veto_id = (TH2F*) finput_SF->Get("Run2017_CutBasedVetoNoIso94XV2"); // Cut-based veto ID
    h_SF_veto_iso = (TH2F*) finput_SF->Get("Run2017_MVAVLooseTightIP2DMini2"); // MiniIso<0.2
    finput_SF_tracking->cd();
    h_SF_tracking = (TH2F*) finput_SF_tracking->Get("EGamma_SF2D");
    // FIXME: Get veto eff. for 2017
    //finput_SF_veto_eff->cd();
    //h_SF_veto_eff = (TH2F*) finput_SF_veto_eff->Get("h2_lepEff_vetoSel_Eff_el");

    // FastSim/FullSim SFs
    finput_SF_FastSim_id->cd(); h_SF_FastSim_id = (TH2F*) finput_SF_FastSim_id->Get("histo2D");
    finput_SF_FastSim_iso->cd(); h_SF_FastSim_iso = (TH2F*) finput_SF_FastSim_iso->Get("histo2D");
    finput_SF_FastSim_veto_id->cd(); h_SF_FastSim_veto_id = (TH2F*) finput_SF_FastSim_veto_id->Get("histo2D");
    finput_SF_FastSim_veto_iso->cd(); h_SF_FastSim_veto_iso = (TH2F*) finput_SF_FastSim_veto_iso->Get("histo2D");
  }
  else if (theDataPeriod == "2018"){
    // FIXME: To be implemented
  }

  curdir->cd();
}

ElectronScaleFactorHandler::~ElectronScaleFactorHandler(){
  this->closeFile(finput_SF);
  this->closeFile(finput_SF_tracking);
  this->closeFile(finput_SF_veto_eff);
  this->closeFile(finput_SF_FastSim_id);
  this->closeFile(finput_SF_FastSim_iso);
  this->closeFile(finput_SF_FastSim_veto_id);
  this->closeFile(finput_SF_FastSim_veto_iso);
}



void ElectronScaleFactorHandler::evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, float const& ele_pt, float const& ele_etasc, TH2F const* hist, bool etaOnY, bool useAbsEta) const{
  if (!hist) return;

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
  if (bc!=0.) be /= bc;
  if (be<0.) be=0;

  theSF *= bc; theSFRelErr = sqrt(pow(theSFRelErr, 2)+pow(be, 2));
}

void ElectronScaleFactorHandler::getIdIsoSFAndError(float& theSF, float& theSFRelErr, float const& ele_pt, float const& ele_etasc, bool isVeto, bool useFastSim) const{
  theSF=1; theSFRelErr=0;

  if (!useFastSim){
    if (!isVeto){
      evalScaleFactorFromHistogram(theSF, theSFRelErr, ele_pt, ele_etasc, h_SF_id, (theDataPeriod=="2016"), false);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, ele_pt, ele_etasc, h_SF_iso, (theDataPeriod=="2016"), false);
    }
    else{
      evalScaleFactorFromHistogram(theSF, theSFRelErr, ele_pt, ele_etasc, h_SF_veto_id, (theDataPeriod=="2016"), false);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, ele_pt, ele_etasc, h_SF_veto_iso, (theDataPeriod=="2016"), false);
    }
  }
  else{
    // FIXME: Need to check axis inversion after 2017 hists are obtained
    if (!isVeto){
      evalScaleFactorFromHistogram(theSF, theSFRelErr, ele_pt, ele_etasc, h_SF_FastSim_id, true/*(theDataPeriod=="2016")*/, true/*(theDataPeriod=="2016")*/);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, ele_pt, ele_etasc, h_SF_FastSim_iso, true/*(theDataPeriod=="2016")*/, true/*(theDataPeriod=="2016")*/);
    }
    else{
      evalScaleFactorFromHistogram(theSF, theSFRelErr, ele_pt, ele_etasc, h_SF_FastSim_veto_id, true/*(theDataPeriod=="2016")*/, true/*(theDataPeriod=="2016")*/);
      evalScaleFactorFromHistogram(theSF, theSFRelErr, ele_pt, ele_etasc, h_SF_FastSim_veto_iso, true/*(theDataPeriod=="2016")*/, true/*(theDataPeriod=="2016")*/);
    }
  }

  // Uncertainty in veto id+iso eff. is doubled as pT thr. is lower in this analysis
  if (isVeto && !useFastSim) theSFRelErr *= 2.f;
}
void ElectronScaleFactorHandler::getRecoSFAndError(float& theSF, float& theSFRelErr, float const& ele_pt, float const& ele_etasc) const{
  theSF=1; theSFRelErr=0;
  evalScaleFactorFromHistogram(theSF, theSFRelErr, ele_pt, ele_etasc, h_SF_tracking, false, false);
}

void ElectronScaleFactorHandler::getGenSFAndError(float& theSF, float& theSFRelErr, float const& ele_pt, float const& ele_eta, float const& theIdIsoSF, float const& theIdIsoSFRelErr) const{
  theSF=1; theSFRelErr=0;
  if (ele_pt<5. || fabs(ele_eta)>2.4) return; // FIXME: MAKE 2.4 AND 5. ADJUSTABLE VALUES (ANALYSIS CUTS)

  // FIXME: Histogram for 2017 might not use abs(eta) or be swapped.
  evalScaleFactorFromHistogram(theSF, theSFRelErr, ele_pt, ele_eta, h_SF_veto_eff, true/*(theDataPeriod=="2016")*/, true/*(theDataPeriod=="2016")*/);

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

