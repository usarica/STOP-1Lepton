#include "ScaleFactorHandlerBase.h"


void ScaleFactorHandlerBase::closeFile(TFile*& f){
  if (f){
    if (f->IsOpen()) f->Close();
    else if (f->IsZombie()) delete f;
  }
  f = nullptr;
}
