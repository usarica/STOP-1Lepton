#include "RooRelBW3ProngPdf.h"
#include "MELAStreamHelpers.hh"
#include <cmath>
#include "TUtil.hh"


using namespace std;
using namespace MELAStreamHelpers;


RooRelBW3ProngPdf::RooRelBW3ProngPdf(
  const char* name, const char* title,
  modelParameters const& parameters,
  modelMeasurables const& measurables,
  RooRelBW3ProngPdf::VDecayType dtype
) :
  RooAbsPdf(name, title),

  kinematicsFlag(dtype),

  mX("mX", "mX", this),
  gamX("gamX", "gamX", this),

  mV("mV", "mV", this),
  gamV("gamV", "gamV", this),

  pT1("pT1", "pT1", this),
  eta1("eta1", "eta1", this),
  phi1("phi1", "phi1", this),
  mass1("mass1", "mass1", this),

  pT2("pT2", "pT2", this),
  eta2("eta2", "eta2", this),
  phi2("phi2", "phi2", this),
  mass2("mass2", "mass2", this),

  pT3("pT3", "pT3", this),
  eta3("eta3", "eta3", this),
  phi3("phi3", "phi3", this),
  mass3("mass3", "mass3", this)
{
  setProxies(parameters, measurables);
}
RooRelBW3ProngPdf::RooRelBW3ProngPdf(const RooRelBW3ProngPdf& other, const char* name) :
  RooAbsPdf(other, name),

  kinematicsFlag(other.kinematicsFlag),

  mX("mX", this, other.mX),
  gamX("gamX", this, other.gamX),

  mV("mV", this, other.mV),
  gamV("gamV", this, other.gamV),

  pT1("pT1", this, other.pT1),
  eta1("eta1", this, other.eta1),
  phi1("phi1", this, other.phi1),
  mass1("mass1", this, other.mass1),

  pT2("pT2", this, other.pT2),
  eta2("eta2", this, other.eta2),
  phi2("phi2", this, other.phi2),
  mass2("mass2", this, other.mass2),

  pT3("pT3", this, other.pT3),
  eta3("eta3", this, other.eta3),
  phi3("phi3", this, other.phi3),
  mass3("mass3", this, other.mass3)
{}


void RooRelBW3ProngPdf::setProxies(modelParameters const& parameters, modelMeasurables const& measurables){
  setProxy(mX, (RooAbsReal*) parameters.mX);
  setProxy(gamX, (RooAbsReal*) parameters.gamX);

  setProxy(mV, (RooAbsReal*) parameters.mV);
  setProxy(gamV, (RooAbsReal*) parameters.gamV);

  setProxy(pT1, (RooAbsReal*) measurables.pT1);
  setProxy(eta1, (RooAbsReal*) measurables.eta1);
  setProxy(phi1, (RooAbsReal*) measurables.phi1);
  setProxy(mass1, (RooAbsReal*) measurables.mass1);

  setProxy(pT2, (RooAbsReal*) measurables.pT2);
  setProxy(eta2, (RooAbsReal*) measurables.eta2);
  setProxy(phi2, (RooAbsReal*) measurables.phi2);
  setProxy(mass2, (RooAbsReal*) measurables.mass2);

  setProxy(pT3, (RooAbsReal*) measurables.pT3);
  setProxy(eta3, (RooAbsReal*) measurables.eta3);
  setProxy(phi3, (RooAbsReal*) measurables.phi3);
  setProxy(mass3, (RooAbsReal*) measurables.mass3);
}
void RooRelBW3ProngPdf::setProxy(RooRealProxy& proxy, RooAbsReal* objectPtr){
  if (objectPtr) proxy.setArg((RooAbsReal&) *objectPtr);
}

Double_t RooRelBW3ProngPdf::evaluate() const{
  Double_t res=0;
  Float_t betaFactor=1;
  const Float_t mx = mX * GeVunit;
  const Float_t gamx = gamX * GeVunit;
  const Float_t qXsq = this->getQXSq(false);
  Float_t denom = pow(qXsq - mx*mx, 2) + pow(mx*gamx, 2);
  if (mV>=0. && gamV>=0.){
    const Float_t mv = mV * GeVunit;
    const Float_t gamv = gamV * GeVunit;
    const Float_t qVsq = this->getQVSq();
    denom *= pow(qVsq - mv*mv, 2) + pow(mv*gamv, 2);
    if (qXsq>qVsq && qVsq>0.f && qXsq>0.f) betaFactor = 1.-qVsq/qXsq;
    else betaFactor = 0;
  }
  if (denom>0.f) res = betaFactor/denom;
  if (kinematicsFlag > kOnlyBW){
    using namespace TUtil;

    TLorentzVector p2, p2_massive; p2_massive.SetPtEtaPhiM(pT1, eta1, phi1, mass1); // Bottom or anti-bottom
    scaleMomentumToEnergy(p2_massive, p2, 0);

    TLorentzVector p3; p3.SetPtEtaPhiM(pT2, eta2, phi2, mass2); // Lepton or anti-lepton
    TLorentzVector p4; p4.SetPtEtaPhiM(pT3, eta3, phi3, mass3); // Neutrino or anti-neutrino
    constrainedRemovePairMass(p3, p4, 0, 0);

    // Scale everything by the GeV unit for computational stability
    p2 *= GeVunit; p3 *= GeVunit; p4 *= GeVunit;

    // Top vector = bottom + f1 + f2
    TLorentzVector p1 = p2 + p3 + p4;
    // In the amplitude below, every vector is outgoing, so top -> -top
    p1 = -p1;

    // Swap for the correct fermion currents
    if (kinematicsFlag == kWplus) std::swap(p3, p4);
    if (kinematicsFlag == kWminus) std::swap(p1, p2);

    // Note that this amplitude omits the longitudinal polarization of the W, which appears as (g_munu - k_mu k_nu/k^2) instead of g_munu in the unitary gauge in order to require no ghost fields.
    // It corresponds to omitting the dampening of the decay amplitude by these longitudinal polarizations.
    // It should be a ~q_W^2/q_t^2 - dependent factor, so not as worrisome as spin correlations.
    const Float_t gV1=0.5, gA1=-0.5, gV2=0.5, gA2=-0.5;
    Float_t kinFactor = (pow(gV1, 2)+pow(gA1, 2))*(pow(gV2, 2)+pow(gA2, 2))*(
      p2.Dot(p4)*p1.Dot(p3) + p2.Dot(p3)*p1.Dot(p4)
      //+ 2*m1*m2*m3*m4 + m1*m2*p3.Dot(p4) + p1.Dot(p2)*m3*m4 
      );
    if (kinematicsFlag != kWany){
      Float_t asymFactor = RooRelBW3ProngPdf::getEuclideanProduct(p2, p4)*RooRelBW3ProngPdf::getEuclideanProduct(p1, p3) - RooRelBW3ProngPdf::getEuclideanProduct(p2, p3)*RooRelBW3ProngPdf::getEuclideanProduct(p1, p4);
      kinFactor += 4.f*gV1*gA1*gV2*gA2*asymFactor;
    }
    if (kinFactor>0.f) res *= kinFactor;
  }
  return res;
}

Float_t RooRelBW3ProngPdf::getQXSq(bool unscaled) const{
  using namespace TUtil;

  // Important to scale everything for m3f
  TLorentzVector p2, p2_massive; p2_massive.SetPtEtaPhiM(pT1, eta1, phi1, mass1); p2=p2_massive;
  if (!unscaled) scaleMomentumToEnergy(p2_massive, p2, 0);

  TLorentzVector p3; p3.SetPtEtaPhiM(pT2, eta2, phi2, mass2);
  TLorentzVector p4; p4.SetPtEtaPhiM(pT3, eta3, phi3, mass3);
  if (!unscaled) constrainedRemovePairMass(p3, p4, 0, 0);

  TLorentzVector p1 = p2 + p3 + p4;
  p1 = p1 * GeVunit;

  const Float_t qsq = p1.M2();
  return qsq;
}
Float_t RooRelBW3ProngPdf::getQVSq() const{
  // Not important to scale everything for m23; the scaling function constrainedRemovePairMass is supposed to preserve qsq.
  const Float_t px2 = pT2 * cos(phi2) * GeVunit;
  const Float_t py2 = pT2 * sin(phi2) * GeVunit;
  const Float_t pz2 = pT2 * sinh(eta2) * GeVunit;
  const Float_t E2 = sqrt(pow(pT2 * cosh(eta2) * GeVunit, 2) + pow(mass2 * GeVunit, 2)*(mass2>=0. ? +1.f : -1.f));

  const Float_t px3 = pT3 * cos(phi3) * GeVunit;
  const Float_t py3 = pT3 * sin(phi3) * GeVunit;
  const Float_t pz3 = pT3 * sinh(eta3) * GeVunit;
  const Float_t E3 = sqrt(pow(pT3 * cosh(eta3) * GeVunit, 2) + pow(mass3 * GeVunit, 2)*(mass3>=0. ? +1.f : -1.f));

  const Float_t px = px2 + px3;
  const Float_t py = py2 + py3;
  const Float_t pz = pz2 + pz3;
  const Float_t E = E2 + E3;

  const Float_t qsq = pow(E, 2) - pow(px, 2) - pow(py, 2) - pow(pz, 2);
  return qsq;
}
Float_t RooRelBW3ProngPdf::getEuclideanProduct(TLorentzVector const& p1, TLorentzVector const& p2){ return (2.*p1.T()*p2.T() - p1.Dot(p2)); }
