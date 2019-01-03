#ifndef SCALEFACTORHANDLERBASE_H
#define SCALEFACTORHANDLERBASE_H

#include "TFile.h"
#include "TH2F.h"
#include "TString.h"


class ScaleFactorHandlerBase{
public:
  ScaleFactorHandlerBase(){}
  virtual ~ScaleFactorHandlerBase(){}

  static void closeFile(TFile*& f);
  static bool getHistogram(TH2F*& h, TFile*& f, TString s);

  virtual bool setup() = 0;
  virtual void reset() = 0;

};


#endif
