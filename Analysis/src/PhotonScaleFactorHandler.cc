#include "PhotonScaleFactorHandler.h"
#include "PhotonSelectionHelpers.h"
#include "TDirectory.h"


using namespace std;
using namespace SampleHelpers;
using namespace PhotonSelectionHelpers;


PhotonScaleFactorHandler::PhotonScaleFactorHandler() :
  ScaleFactorHandlerBase(),
  finput_SF(nullptr),
  h_SF_id(nullptr)
{
  setup();
}

PhotonScaleFactorHandler::~PhotonScaleFactorHandler(){ this->reset(); }


bool PhotonScaleFactorHandler::setup(){
  bool res = true;
  TDirectory* curdir = gDirectory;

  this->reset();

  // More info: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2#Photon_efficiencies_and_scale_fa
  if (theDataYear == 2016){
    // ID/Iso. SF files
    finput_SF = TFile::Open(STOP1LPKGDATAPATH+"PhotonSFs/Fall17V2_2016_Tight_photons.root", "read");

    res = (
      getHistogram(h_SF_id, finput_SF, "EGamma_SF2D")
      );
  }
  else if (theDataYear == 2017){
    // ID/Iso. SF files
    finput_SF = TFile::Open(STOP1LPKGDATAPATH+"PhotonSFs/2017_PhotonsTight.root", "read");

    res = (
      getHistogram(h_SF_id, finput_SF, "EGamma_SF2D")
      );
  }
  else if (theDataYear == 2018){
    // ID/Iso. SF files
    finput_SF = TFile::Open(STOP1LPKGDATAPATH+"PhotonSFs/2018_PhotonsTight.root", "read");

    res = (
      getHistogram(h_SF_id, finput_SF, "EGamma_SF2D")
      );
  }

  curdir->cd();

  return res;
}
void PhotonScaleFactorHandler::reset(){
  ScaleFactorHandlerBase::closeFile(finput_SF);

  h_SF_id = nullptr;
}

void PhotonScaleFactorHandler::evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, PhotonObject const* obj, TH2F const* hist, bool etaOnY, bool useAbsEta) const{
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

void PhotonScaleFactorHandler::getIdIsoSFAndError(float& theSF, float& theSFRelErr, PhotonObject const* obj, bool /*useFastSim*/) const{
  theSF=1; theSFRelErr=0;

  if (!obj) return;
  bool passSel = obj->testSelection(bit_preselection_idisoreco);
  if (!passSel) return;

  evalScaleFactorFromHistogram(theSF, theSFRelErr, obj, h_SF_id, false, false);
}
