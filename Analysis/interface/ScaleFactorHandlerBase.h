#ifndef SCALEFACTORHANDLERBASE_H
#define SCALEFACTORHANDLERBASE_H


class ScaleFactorHandlerBase{
public:
  ScaleFactorHandlerBase(){}
  virtual ~ScaleFactorHandlerBase(){}

  static void closeFile(TFile* f){ if (f && f->IsOpen()) f->Close(); }

};



#endif
