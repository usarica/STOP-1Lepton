#ifndef PHOTONSCALEFACTORHANDLER_H
#define PHOTONSCALEFACTORHANDLER_H

#include "TH2F.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "Samples.h"
#include "ScaleFactorHandlerBase.h"
#include "PhotonObject.h"


class PhotonScaleFactorHandler : public ScaleFactorHandlerBase{
protected:
  TFile* finput_SF;

  TH2F* h_SF_id;

  void evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, PhotonObject const* obj, TH2F const* hist, bool etaOnY, bool useAbsEta) const;

public:
  PhotonScaleFactorHandler();
  ~PhotonScaleFactorHandler();

  bool setup();
  void reset();

  void getIdIsoSFAndError(float& theSF, float& theSFRelErr, PhotonObject const* obj, bool useFastSim) const;

};



#endif
