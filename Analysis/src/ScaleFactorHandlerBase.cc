#include "ScaleFactorHandlerBase.h"
#include "TDirectory.h"


void ScaleFactorHandlerBase::closeFile(TFile*& f){
  if (f){
    if (f->IsOpen()) f->Close();
    else delete f;
  }
  f = nullptr;
}

bool ScaleFactorHandlerBase::getHistogram(TH2F*& h, TFile*& f, TString s){
  TDirectory* curdir = gDirectory;
  if (s=="") return false;
  if (!f) return false;
  if (!f->IsOpen()) return false;
  if (f->IsZombie()) return false;

  f->cd();
  h = (TH2F*) f->Get(s);
  curdir->cd();
  return (h!=nullptr);
}
