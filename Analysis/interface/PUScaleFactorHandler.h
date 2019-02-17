#ifndef PUSCALEFACTORHANDLER_H
#define PUSCALEFACTORHANDLER_H

#include "TH1D.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "ScaleFactorHandlerBase.h"


class PUScaleFactorHandler : public ScaleFactorHandlerBase{
protected:
  TFile* finput;

  TH1D* h_nominal;
  TH1D* h_up;
  TH1D* h_dn;

public:
  PUScaleFactorHandler();
  ~PUScaleFactorHandler();

  bool setup();
  void reset();

  void getPileUpWeight(float& theSF, float& theSFUp, float& theSFDn, int ntruevtxs) const;

};



#endif
