#ifndef SCALEFACTORHANDLERBASE_H
#define SCALEFACTORHANDLERBASE_H

#include "TFile.h"


class ScaleFactorHandlerBase{
public:
  ScaleFactorHandlerBase(){}
  virtual ~ScaleFactorHandlerBase(){}

  static void closeFile(TFile*& f);

  virtual bool setup() = 0;
  virtual void reset() = 0;

};


#endif
