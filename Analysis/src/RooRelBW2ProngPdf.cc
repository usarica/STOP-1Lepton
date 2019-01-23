#include "RooRelBW2ProngPdf.h"
#include "MELAStreamHelpers.hh"
#include <cmath>


using namespace std;
using namespace MELAStreamHelpers;


RooRelBW2ProngPdf::RooRelBW2ProngPdf(
  const char* name, const char* title,
  modelParameters const& parameters,
  modelMeasurables const& measurables
) :
  RooAbsPdf(name, title),

  mX("mX", "mX", this),
  gamX("gamX", "gamX", this),

  pT1("pT1", "pT1", this),
  eta1("eta1", "eta1", this),
  phi1("phi1", "phi1", this),
  mass1("mass1", "mass1", this),

  pT2("pT2", "pT2", this),
  eta2("eta2", "eta2", this),
  phi2("phi2", "phi2", this),
  mass2("mass2", "mass2", this)
{
  setProxies(parameters, measurables);
}
RooRelBW2ProngPdf::RooRelBW2ProngPdf(const RooRelBW2ProngPdf& other, const char* name) :
  RooAbsPdf(other, name),

  mX("mX", this, other.mX),
  gamX("gamX", this, other.gamX),

  pT1("pT1", this, other.pT1),
  eta1("eta1", this, other.eta1),
  phi1("phi1", this, other.phi1),
  mass1("mass1", this, other.mass1),

  pT2("pT2", this, other.pT2),
  eta2("eta2", this, other.eta2),
  phi2("phi2", this, other.phi2),
  mass2("mass2", this, other.mass2)
{}


void RooRelBW2ProngPdf::setProxies(modelParameters const& parameters, modelMeasurables const& measurables){
  setProxy(mX, (RooAbsReal*) parameters.mX);
  setProxy(gamX, (RooAbsReal*) parameters.gamX);

  setProxy(pT1, (RooAbsReal*) measurables.pT1);
  setProxy(eta1, (RooAbsReal*) measurables.eta1);
  setProxy(phi1, (RooAbsReal*) measurables.phi1);
  setProxy(mass1, (RooAbsReal*) measurables.mass1);

  setProxy(pT2, (RooAbsReal*) measurables.pT2);
  setProxy(eta2, (RooAbsReal*) measurables.eta2);
  setProxy(phi2, (RooAbsReal*) measurables.phi2);
  setProxy(mass2, (RooAbsReal*) measurables.mass2);
}
void RooRelBW2ProngPdf::setProxy(RooRealProxy& proxy, RooAbsReal* objectPtr){
  if (objectPtr) proxy.setArg((RooAbsReal&) *objectPtr);
}

Double_t RooRelBW2ProngPdf::evaluate() const{
  Double_t res=0;
  const float mx = mX * GeVunit;
  const float gamx = gamX * GeVunit;

  const float px1 = pT1 * cos(phi1) * GeVunit;
  const float py1 = pT1 * sin(phi1) * GeVunit;
  const float pz1 = pT1 * sinh(eta1) * GeVunit;
  const float E1 = sqrt(pow(pT1 * cosh(eta1) * GeVunit, 2) + pow(mass1 * GeVunit, 2)*(mass1>=0. ? +1.f : -1.f));

  const float px2 = pT2 * cos(phi2) * GeVunit;
  const float py2 = pT2 * sin(phi2) * GeVunit;
  const float pz2 = pT2 * sinh(eta2) * GeVunit;
  const float E2 = sqrt(pow(pT2 * cosh(eta2) * GeVunit, 2) + pow(mass2 * GeVunit, 2)*(mass2>=0. ? +1.f : -1.f));

  const float px = px1 + px2;
  const float py = py1 + py2;
  const float pz = pz1 + pz2;
  const float E = E1 + E2;

  const float qsq = pow(E, 2) - pow(px, 2) - pow(py, 2) - pow(pz, 2);
  const float denom = pow(qsq - mx*mx, 2) + pow(mx*gamx, 2);
  if (denom>0.f) res = 1.f/denom;

  return res;
}


