#include "Samples.h"
#include "PUScaleFactorHandler.h"


using namespace std;
using namespace SampleHelpers;


PUScaleFactorHandler::PUScaleFactorHandler() :
  ScaleFactorHandlerBase(),

  finput(nullptr),

  h_nominal(nullptr),
  h_up(nullptr),
  h_dn(nullptr)
{
  this->setup();
}
PUScaleFactorHandler::~PUScaleFactorHandler(){
  this->reset();
}

bool PUScaleFactorHandler::setup(){
  bool res = true;
  this->reset();

  finput = TFile::Open(STOP1LPKGDATAPATH+"PileUpSFs/puWeights_Run2.root", "read");
  res = (
    finput && finput->IsOpen()
    && getHistogram(h_nominal, finput, Form("puWeight%i", SampleHelpers::theDataYear))
    && getHistogram(h_up, finput, Form("puWeight%iUp", SampleHelpers::theDataYear))
    && getHistogram(h_dn, finput, Form("puWeight%iDown", SampleHelpers::theDataYear))
    );

  return res;
}
void PUScaleFactorHandler::reset(){
  ScaleFactorHandlerBase::closeFile(finput);

  h_nominal = nullptr;
  h_up = nullptr;
  h_dn = nullptr;
}
void PUScaleFactorHandler::getPileUpWeight(float& theSF, float& theSFUp, float& theSFDn, int ntruevtxs) const{
  theSF = theSFUp = theSFDn = 1.f;

  int ibin = h_nominal->GetXaxis()->FindBin(ntruevtxs);
  if (ibin==0) ibin++;
  if (ibin==h_nominal->GetNbinsX()+1) ibin--;

  theSF = h_nominal->GetBinContent(ibin);
  theSFUp = h_up->GetBinContent(ibin);
  theSFDn = h_dn->GetBinContent(ibin);
}
