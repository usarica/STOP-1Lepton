#ifndef ROORELBW2PRONGPDF_H
#define ROORELBW2PRONGPDF_H

#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"


class RooRelBW2ProngPdf : public RooAbsPdf{
public:
  static constexpr Double_t GeVunit=1e-2;

  struct modelParameters{
    RooAbsReal* mX;
    RooAbsReal* gamX;
  };
  struct modelMeasurables{
    RooAbsReal* pT1;
    RooAbsReal* eta1;
    RooAbsReal* phi1;
    RooAbsReal* mass1;

    RooAbsReal* pT2;
    RooAbsReal* eta2;
    RooAbsReal* phi2;
    RooAbsReal* mass2;
  };

  RooRelBW2ProngPdf() : RooAbsPdf(){}
  RooRelBW2ProngPdf(
    const char* name, const char* title,
    modelParameters const& parameters,
    modelMeasurables const& measurables
  );
  RooRelBW2ProngPdf(const RooRelBW2ProngPdf& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooRelBW2ProngPdf(*this, newname); }
  inline virtual ~RooRelBW2ProngPdf(){}

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const { return 0; }
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const { return 1.; }

  Float_t getQSq() const;

protected:
  RooRealProxy mX;
  RooRealProxy gamX;

  RooRealProxy pT1;
  RooRealProxy eta1;
  RooRealProxy phi1;
  RooRealProxy mass1;

  RooRealProxy pT2;
  RooRealProxy eta2;
  RooRealProxy phi2;
  RooRealProxy mass2;

  Double_t evaluate() const;

  virtual void setProxies(modelParameters const& parameters, modelMeasurables const& measurables);
  virtual void setProxy(RooRealProxy& proxy, RooAbsReal* objectPtr);

};

#endif