#ifndef ROORELBW3PRONGPDF_H
#define ROORELBW3PRONGPDF_H

#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "TLorentzVector.h"


class RooRelBW3ProngPdf : public RooAbsPdf{
public:
  enum VDecayType{
    kOnlyBW=0,
    kWany,
    kWplus, // Top, p2 is neutrino, p3 is l+
    kWminus // ATop, p2 is l-, p3 is anti-neutrino
  };

protected:
  VDecayType kinematicsFlag;

public:
  static constexpr Double_t GeVunit=1e-2;

  struct modelParameters{
    RooAbsReal* mX;
    RooAbsReal* gamX;

    RooAbsReal* mV;
    RooAbsReal* gamV;
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

    RooAbsReal* pT3;
    RooAbsReal* eta3;
    RooAbsReal* phi3;
    RooAbsReal* mass3;
  };

  RooRelBW3ProngPdf() : RooAbsPdf(){}
  RooRelBW3ProngPdf(
    const char* name, const char* title,
    modelParameters const& parameters,
    modelMeasurables const& measurables,
    RooRelBW3ProngPdf::VDecayType dtype = RooRelBW3ProngPdf::kOnlyBW
  );
  RooRelBW3ProngPdf(const RooRelBW3ProngPdf& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooRelBW3ProngPdf(*this, newname); }
  inline virtual ~RooRelBW3ProngPdf(){}

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const { return 0; }
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const { return 1.; }

  Float_t getQXSq(bool unscaled=false) const;
  Float_t getQVSq() const;

  void setKinematicsFlag(VDecayType flag){ kinematicsFlag=flag; }

  static Float_t getEuclideanProduct(TLorentzVector const& p1, TLorentzVector const& p2);

protected:
  RooRealProxy mX;
  RooRealProxy gamX;

  RooRealProxy mV;
  RooRealProxy gamV;

  // Bottom/anti-bottom
  RooRealProxy pT1;
  RooRealProxy eta1;
  RooRealProxy phi1;
  RooRealProxy mass1;

  // Jet or Lepton
  RooRealProxy pT2;
  RooRealProxy eta2;
  RooRealProxy phi2;
  RooRealProxy mass2;

  // Jet or neutrino
  RooRealProxy pT3;
  RooRealProxy eta3;
  RooRealProxy phi3;
  RooRealProxy mass3;

  Double_t evaluate() const;

  virtual void setProxies(modelParameters const& parameters, modelMeasurables const& measurables);
  virtual void setProxy(RooRealProxy& proxy, RooAbsReal* objectPtr);

};

#endif