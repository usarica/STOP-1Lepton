#include <cassert>
#include "common_includes.h"
#include "PDGHelpers.h"
#include "MELAParticle.h"
#include "RooMsgService.h"
#include "RooRelBW2ProngPdf.h"
#include "RooRelBW3ProngPdf.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "TNumericUtil.hh"
#include <HiggsAnalysis/CombinedLimit/interface/AsymPow.h>
#include <limits>
#include "TStyle.h"
#include "TLegend.h"


#define DATA_CLASS(t, name, ...) \
  t* name=nullptr; \
  tree->SetBranchStatus(#name, 1); \
  tree->SetBranchAddress(#name, &name);
#define DATA_SIMPLE(t, name, ...) \
  t name=0; \
  tree->SetBranchStatus(#name, 1); \
  tree->SetBranchAddress(#name, &name);
#define BOOK_BRANCH(t, name, outtree, ...) \
  t name; (outtree).Branch(#name, &name);


const bool doSingleTheta=true;


int getHistogramColorByWeightType(WeightVariables::WeightType type){
  switch (type){
  case WeightVariables::wCentral:
  case WeightVariables::wCentral_Default:
    return kBlack;

  case WeightVariables::wFacScaleUp:
  case WeightVariables::wFacScaleDn:
    return kRed;

  case WeightVariables::wRenScaleUp:
  case WeightVariables::wRenScaleDn:
    return kBlue;

  case WeightVariables::wPDFUp:
  case WeightVariables::wPDFDn:
    return kSpring+9;

  case WeightVariables::wAsMZUp:
  case WeightVariables::wAsMZDn:
    return kViolet;

  case WeightVariables::wPSUp:
  case WeightVariables::wPSDn:
    return kPink-7;

  case WeightVariables::wPDFUp_Default:
  case WeightVariables::wPDFDn_Default:
    return kCyan-3;
  case WeightVariables::wISRUp:
  case WeightVariables::wISRDn:
    return kBlue;
  case WeightVariables::wFSRUp:
  case WeightVariables::wFSRDn:
    return kViolet;
  default:
    return kBlack;
  };
}
int getHistogramDashByWeightType(WeightVariables::WeightType type){
  switch (type){
  case WeightVariables::wCentral:
    return 1;

  case WeightVariables::wCentral_Default:
    return 2;

  case WeightVariables::wFacScaleUp:
  case WeightVariables::wRenScaleUp:
  case WeightVariables::wPDFUp:
  case WeightVariables::wAsMZUp:
  case WeightVariables::wPSUp:
  case WeightVariables::wPDFUp_Default:
  case WeightVariables::wISRUp:
  case WeightVariables::wFSRUp:
    return 1;

  case WeightVariables::wFacScaleDn:
  case WeightVariables::wRenScaleDn:
  case WeightVariables::wPDFDn:
  case WeightVariables::wAsMZDn:
  case WeightVariables::wPSDn:
  case WeightVariables::wPDFDn_Default:
  case WeightVariables::wISRDn:
  case WeightVariables::wFSRDn:
    return 7;

  default:
    return 1;
  };
}

bool getBestZCandidate(std::vector<MELAParticle>& leptons, MELAParticle& bestZ){
  std::vector<MELAParticle*> leps, aleps;
  for (auto& l:leptons){
    if (l.id>0) leps.push_back(&l);
    else if (l.id<0) aleps.push_back(&l);
  }
  std::pair<MELAParticle*, MELAParticle*> Zpair=std::pair<MELAParticle*, MELAParticle*>(nullptr, nullptr);
  for (auto& lep:leps){
    for (auto& alep:aleps){
      if (PDGHelpers::getCoupledVertex(lep->id, alep->id)!=23) continue;
      if (
        (!Zpair.first || !Zpair.second)
        ||
        fabs((lep->p4 + alep->p4).M()-PDGHelpers::Zmass)<fabs((Zpair.first->p4 + Zpair.second->p4).M()-PDGHelpers::Zmass)
        ){ Zpair.first = lep; Zpair.second = alep; }
    }
  }
  if (Zpair.first && Zpair.second){
    MELAParticle tmpZ(23, Zpair.first->p4 + Zpair.second->p4);
    tmpZ.addDaughter(Zpair.first);
    tmpZ.addDaughter(Zpair.second);
    bestZ.swap(tmpZ);
    return true;
  }
  else return false;
}
void doMinimization(RooAbsReal& fcn, std::vector<RooRealVar>& thetas, bool& isFitOK, int& fitstatus){
  RooFitResult* fitResult=nullptr;
  for (auto& var:thetas) var.setVal(0);

  RooMinuit minimizer(fcn);
  minimizer.setPrintLevel(-1);
  minimizer.setNoWarn();
  minimizer.setStrategy(2);
  minimizer.migrad();
  fitResult = minimizer.save();
  if (fitResult){
    fitstatus = fitResult->status();
    isFitOK = (fitstatus==0);
    delete fitResult; fitResult=nullptr;
  }
  if (fitstatus==4){
    for (unsigned int itry=0; itry<5; itry++){
      minimizer.migrad();
      fitResult = minimizer.save();
      if (fitResult){
        fitstatus=fitResult->status();
        isFitOK=(fitstatus==0);
        delete fitResult; fitResult=nullptr;
      }
      if (isFitOK) break;
    }
  }
  else{
    minimizer.setStrategy(0);
    for (unsigned int itry=0; itry<5; itry++){
      minimizer.migrad();
      fitResult = minimizer.save();
      if (fitResult){
        fitstatus=fitResult->status();
        isFitOK=(fitstatus==0);
        delete fitResult; fitResult=nullptr;
      }
      if (isFitOK) break;
    }
  }
}

struct FitResultSummary{
  int fitstatus;
  int objid;
  float pdfval;

  float constraintval_Xpdf;
  float Xmass;
  float Xmass_prefit;

  float constraintval_Vpdf;
  float Vmass;
  float Vmass_prefit;

  std::vector<unsigned int> jetindices;
  std::vector<float> jetnuisances;


  FitResultSummary();
  FitResultSummary& operator=(FitResultSummary const& other);

  bool operator==(FitResultSummary const& other) const;
  bool operator<(FitResultSummary const& other) const;
  bool operator>(FitResultSummary const& other) const;
  bool operator<=(FitResultSummary const& other) const;
  bool operator>=(FitResultSummary const& other) const;
};
FitResultSummary::FitResultSummary() :
  fitstatus(-1),
  objid(-9000),
  pdfval(-1),

  constraintval_Xpdf(-1),
  Xmass(-1),
  Xmass_prefit(-1),

  constraintval_Vpdf(-1),
  Vmass(-1),
  Vmass_prefit(-1)
{}
FitResultSummary& FitResultSummary::operator=(FitResultSummary const& other){
  this->fitstatus = other.fitstatus;
  this->objid = other.objid;
  this->pdfval = other.pdfval;

  this->constraintval_Xpdf = other.constraintval_Xpdf;
  this->Xmass = other.Xmass;
  this->Xmass_prefit = other.Xmass_prefit;

  this->constraintval_Vpdf = other.constraintval_Vpdf;
  this->Vmass = other.Vmass;
  this->Vmass_prefit = other.Vmass_prefit;

  this->jetindices = other.jetindices;
  this->jetnuisances = other.jetnuisances;

  return (*this);
}
bool FitResultSummary::operator==(FitResultSummary const& other) const{
  if (this->objid!=other.objid) return false;
  if (this->fitstatus<0 && other.fitstatus<0) return true;
  else if (this->fitstatus<0 && other.fitstatus>=0) return false;
  else if (this->fitstatus>=0 && other.fitstatus<0) return false;
  else if (this->fitstatus>0 && other.fitstatus==0) return false;
  else if (this->fitstatus==0 && other.fitstatus>0) return false;
  else if (other.fitstatus==this->fitstatus) return (this->pdfval==other.pdfval);
  return false;
}
bool FitResultSummary::operator<(FitResultSummary const& other) const{
  if (this->objid!=other.objid) return false;
  if (this->fitstatus<0 && other.fitstatus<0) return false;
  else if (this->fitstatus<0 && other.fitstatus>=0) return true;
  else if (this->fitstatus>=0 && other.fitstatus<0) return false;
  else if (this->fitstatus>0 && other.fitstatus==0) return true;
  else if (this->fitstatus==0 && other.fitstatus>0) return false;
  else if (other.fitstatus==this->fitstatus) return (this->pdfval<other.pdfval);
  return false;
}
bool FitResultSummary::operator>(FitResultSummary const& other) const{
  if (this->objid!=other.objid) return false;
  if (this->fitstatus<0 && other.fitstatus<0) return false;
  else if (this->fitstatus<0 && other.fitstatus>=0) return false;
  else if (this->fitstatus>=0 && other.fitstatus<0) return true;
  else if (this->fitstatus>0 && other.fitstatus==0) return false;
  else if (this->fitstatus==0 && other.fitstatus>0) return true;
  else if (other.fitstatus==this->fitstatus) return (this->pdfval>other.pdfval);
  return false;
}
bool FitResultSummary::operator<=(FitResultSummary const& other) const{ return ((*this)<other || (*this)==other); }
bool FitResultSummary::operator>=(FitResultSummary const& other) const{ return ((*this)>other || (*this)==other); }


void getAllHadTopFits(std::vector<std::pair<MELAParticle, float>> const& ak4jets, std::vector<FitResultSummary>& res){
  constexpr size_t nprongs = 3;
  unsigned int const nak4jets_preselected = ak4jets.size();
  if (nak4jets_preselected<nprongs) return;

  RooConstVar one_rcv("one_rcv", "one_rcv", 1);
  RooConstVar zero_rcv("zero_rcv", "zero_rcv", 0);

  RooConstVar mT_3prong("mT_3prong", "mT", PDGHelpers::Topmass);
  RooConstVar gamT_3prong("gamT_3prong", "gamT", PDGHelpers::Topwidth);
  RooConstVar mW_3prong("mW_3prong", "mW", PDGHelpers::Wmass);
  RooConstVar gamW_3prong("gamW_3prong", "gamW", PDGHelpers::Wwidth);

  RooRelBW3ProngPdf::modelParameters Tpars;
  Tpars.mX = &mT_3prong;
  Tpars.gamX = &gamT_3prong;
  Tpars.mV = &mW_3prong;
  Tpars.gamV = &gamW_3prong;

  std::vector<std::vector<int>> comblist;
  TNumericUtil::PermutationGenerator(nak4jets_preselected, nprongs, comblist, 0);
  for (auto const& comb:comblist){
    if (comb.at(1)>comb.at(2)) continue;

    RooRelBW3ProngPdf::modelMeasurables vars;

    std::vector<RooConstVar> ptraw, eta, phi, massraw, kappaLow, kappaHigh;
    ptraw.reserve(nprongs); eta.reserve(nprongs); phi.reserve(nprongs); massraw.reserve(nprongs); kappaLow.reserve(nprongs); kappaHigh.reserve(nprongs);
    std::vector<RooRealVar> theta_ptunc; theta_ptunc.reserve(nprongs); if (doSingleTheta) theta_ptunc.emplace_back("theta_ptunc_alljets", "theta_ptunc_alljets", 0., -5., 5.);
    std::vector<AsymPow> asympow_mult; asympow_mult.reserve(nprongs);
    std::vector<RooFormulaVar> pt_rfv, mass_rfv;
    pt_rfv.reserve(nprongs); mass_rfv.reserve(nprongs);
    std::vector<RooGaussian> pdf_ptunc; pdf_ptunc.reserve(nprongs);
    for (unsigned int ijet=0; ijet<nprongs; ijet++){
      const int& ipart = comb.at(ijet);
      std::pair<MELAParticle, float> const& thePair = ak4jets.at(ipart);
      MELAParticle const& part = thePair.first;
      float const& ptRes = thePair.second;
      float const relPtRes = (part.pt()>0.f ? ptRes/part.pt() : 0.f);

      kappaLow.emplace_back(Form("kappaLow_jet%i", ijet+1), Form("kappaLow_jet%i", ijet+1), 1.f/(1.f+relPtRes));
      kappaHigh.emplace_back(Form("kappaHigh_jet%i", ijet+1), Form("kappaHigh_jet%i", ijet+1), 1.f+relPtRes);
      if (!doSingleTheta) theta_ptunc.emplace_back(Form("theta_ptunc_jet%i", ijet+1), Form("theta_ptunc_jet%i", ijet+1), 0., -5., 5.);
      asympow_mult.emplace_back(Form("asympow_mult_jet%i", ijet+1), Form("asympow_mult_jet%i", ijet+1), kappaLow.back(), kappaHigh.back(), theta_ptunc.back());
      pdf_ptunc.emplace_back(Form("pdf_ptunc_jet%i", ijet+1), Form("pdf_ptunc_jet%i", ijet+1), theta_ptunc.back(), zero_rcv, one_rcv);

      eta.emplace_back(Form("eta_jet%i", ijet+1), Form("eta_jet%i", ijet+1), ak4jets.at(ipart).first.eta());
      phi.emplace_back(Form("phi_jet%i", ijet+1), Form("phi_jet%i", ijet+1), ak4jets.at(ipart).first.phi());

      ptraw.emplace_back(Form("ptraw_jet%i", ijet+1), Form("ptraw_jet%i", ijet+1), ak4jets.at(ipart).first.pt());
      pt_rfv.emplace_back(Form("pt_rfv_jet%i", ijet+1), "(@0*@1)", RooArgList(ptraw.back(), asympow_mult.back()));

      massraw.emplace_back(Form("massraw_jet%i", ijet+1), Form("massraw_jet%i", ijet+1), ak4jets.at(ipart).first.m());
      mass_rfv.emplace_back(Form("mass_rfv_jet%i", ijet+1), "(@0*@1)", RooArgList(massraw.back(), asympow_mult.back()));
    }

    vars.pT1 = &(pt_rfv.at(0));
    vars.eta1 = &(eta.at(0));
    vars.phi1 = &(phi.at(0));
    vars.mass1 = &(mass_rfv.at(0));

    vars.pT2 = &(pt_rfv.at(1));
    vars.eta2 = &(eta.at(1));
    vars.phi2 = &(phi.at(1));
    vars.mass2 = &(mass_rfv.at(1));

    vars.pT3 = &(pt_rfv.at(2));
    vars.eta3 = &(eta.at(2));
    vars.phi3 = &(phi.at(2));
    vars.mass3 = &(mass_rfv.at(2));

    RooProdPdf constraintspdf("constraintspdf", "constraintspdf", RooArgList(pdf_ptunc.at(0), pdf_ptunc.at(1), pdf_ptunc.at(2)));
    RooProdPdf constraintspdf_onlyV("constraintspdf_onlyV", "constraintspdf_onlyV", RooArgList(pdf_ptunc.at(1), pdf_ptunc.at(2))); // Constraints only for the V->ff side of the pdf

    RooRelBW3ProngPdf topWpdf("topWpdf", "topWpdf", Tpars, vars, RooRelBW3ProngPdf::kWany);
    RooProdPdf totalpdf("totalpdf", "totalpdf", RooArgList(topWpdf, constraintspdf));
    RooFormulaVar logL_rfv("logLikelihood_rfv", "-log(@0)", RooArgList(totalpdf));

    float topmass_prefit = topWpdf.getQXSq(); // Not important from which pdf you read it
    float Wmass_prefit = topWpdf.getQVSq(); // Not important from which pdf you read it

    bool const doWfits = (Wmass_prefit>0.f && sqrt(Wmass_prefit)/RooRelBW3ProngPdf::GeVunit>=40.f && sqrt(Wmass_prefit)/RooRelBW3ProngPdf::GeVunit<=210.f);
    bool const doTWfits = doWfits && (topmass_prefit>0.f && sqrt(topmass_prefit)/RooRelBW3ProngPdf::GeVunit>=100.f && sqrt(topmass_prefit)/RooRelBW3ProngPdf::GeVunit<=350.f);

    bool isFitOK=false;
    int fitstatus=-999;
    // Do T+W minimization
    if (doTWfits) doMinimization(logL_rfv, theta_ptunc, isFitOK, fitstatus);
    if (isFitOK){
      res.push_back(FitResultSummary());
      FitResultSummary& lastFit = res.back();

      lastFit.fitstatus = fitstatus;
      lastFit.pdfval = totalpdf.getVal();
      lastFit.objid = 6;

      lastFit.constraintval_Xpdf = constraintspdf.getVal();
      lastFit.Xmass = sqrt(topWpdf.getQXSq())*100.f;
      lastFit.Xmass_prefit = sqrt(topmass_prefit)*100.f;

      lastFit.constraintval_Vpdf = constraintspdf_onlyV.getVal();
      lastFit.Vmass = sqrt(topWpdf.getQVSq())*100.f;
      lastFit.Vmass_prefit = sqrt(Wmass_prefit)*100.f;

      for (size_t ijet=0; ijet<nprongs; ijet++){
        lastFit.jetnuisances.push_back((!doSingleTheta ? theta_ptunc.at(ijet) : theta_ptunc.back()).getVal());
        lastFit.jetindices.push_back(comb.at(ijet));
      }
    }
  } // End loop over permutations
  std::sort(res.begin(), res.end(), [](const FitResultSummary& a, const FitResultSummary& b){return a > b; });
}
void getAllHadVFits(int Vid, std::vector<std::pair<MELAParticle, float>> const& ak4jets, std::vector<FitResultSummary>& res){
  constexpr size_t nprongs = 2;
  unsigned int const nak4jets_preselected = ak4jets.size();
  if (nak4jets_preselected<nprongs) return;

  RooConstVar one_rcv("one_rcv", "one_rcv", 1);
  RooConstVar zero_rcv("zero_rcv", "zero_rcv", 0);

  RooConstVar mZ("mZ", "mZ", PDGHelpers::Zmass);
  RooConstVar gamZ("gamZ", "gamZ", PDGHelpers::Zwidth);
  RooConstVar mW("mW", "mW", PDGHelpers::Wmass);
  RooConstVar gamW("gamW", "gamW", PDGHelpers::Wwidth);

  RooRelBW2ProngPdf::modelParameters Vpars;
  if (Vid==23){
    Vpars.mX = &mZ;
    Vpars.gamX = &gamZ;
  }
  else if (abs(Vid)==24){
    Vpars.mX = &mW;
    Vpars.gamX = &gamW;
  }
  else assert(0);

  std::vector<std::vector<int>> comblist;
  TNumericUtil::CombinationGenerator(nak4jets_preselected, nprongs, comblist, 0);
  for (auto const& comb:comblist){
    RooRelBW2ProngPdf::modelMeasurables vars;

    std::vector<RooConstVar> ptraw, eta, phi, massraw, kappaLow, kappaHigh;
    ptraw.reserve(nprongs); eta.reserve(nprongs); phi.reserve(nprongs); massraw.reserve(nprongs); kappaLow.reserve(nprongs); kappaHigh.reserve(nprongs);
    std::vector<RooRealVar> theta_ptunc; theta_ptunc.reserve(nprongs); if (doSingleTheta) theta_ptunc.emplace_back("theta_ptunc_alljets", "theta_ptunc_alljets", 0., -5., 5.);
    std::vector<AsymPow> asympow_mult; asympow_mult.reserve(nprongs);
    std::vector<RooFormulaVar> pt_rfv, mass_rfv;
    pt_rfv.reserve(nprongs); mass_rfv.reserve(nprongs);
    std::vector<RooGaussian> pdf_ptunc; pdf_ptunc.reserve(nprongs);
    for (unsigned int ijet=0; ijet<nprongs; ijet++){
      const int& ipart = comb.at(ijet);
      std::pair<MELAParticle, float> const& thePair = ak4jets.at(ipart);
      MELAParticle const& part = thePair.first;
      float const& ptRes = thePair.second;
      float const relPtRes = (part.pt()>0.f ? ptRes/part.pt() : 0.f);

      kappaLow.emplace_back(Form("kappaLow_jet%i", ijet+1), Form("kappaLow_jet%i", ijet+1), 1.f/(1.f+relPtRes));
      kappaHigh.emplace_back(Form("kappaHigh_jet%i", ijet+1), Form("kappaHigh_jet%i", ijet+1), 1.f+relPtRes);
      if (!doSingleTheta) theta_ptunc.emplace_back(Form("theta_ptunc_jet%i", ijet+1), Form("theta_ptunc_jet%i", ijet+1), 0., -5., 5.);
      asympow_mult.emplace_back(Form("asympow_mult_jet%i", ijet+1), Form("asympow_mult_jet%i", ijet+1), kappaLow.back(), kappaHigh.back(), theta_ptunc.back());
      pdf_ptunc.emplace_back(Form("pdf_ptunc_jet%i", ijet+1), Form("pdf_ptunc_jet%i", ijet+1), theta_ptunc.back(), zero_rcv, one_rcv);

      eta.emplace_back(Form("eta_jet%i", ijet+1), Form("eta_jet%i", ijet+1), ak4jets.at(ipart).first.eta());
      phi.emplace_back(Form("phi_jet%i", ijet+1), Form("phi_jet%i", ijet+1), ak4jets.at(ipart).first.phi());

      ptraw.emplace_back(Form("ptraw_jet%i", ijet+1), Form("ptraw_jet%i", ijet+1), ak4jets.at(ipart).first.pt());
      pt_rfv.emplace_back(Form("pt_rfv_jet%i", ijet+1), "(@0*@1)", RooArgList(ptraw.back(), asympow_mult.back()));

      massraw.emplace_back(Form("massraw_jet%i", ijet+1), Form("massraw_jet%i", ijet+1), ak4jets.at(ipart).first.m());
      mass_rfv.emplace_back(Form("mass_rfv_jet%i", ijet+1), "(@0*@1)", RooArgList(massraw.back(), asympow_mult.back()));
    }

    vars.pT1 = &(pt_rfv.at(0));
    vars.eta1 = &(eta.at(0));
    vars.phi1 = &(phi.at(0));
    vars.mass1 = &(mass_rfv.at(0));

    vars.pT2 = &(pt_rfv.at(1));
    vars.eta2 = &(eta.at(1));
    vars.phi2 = &(phi.at(1));
    vars.mass2 = &(mass_rfv.at(1));

    RooProdPdf constraintspdf("constraintspdf", "constraintspdf", RooArgList(pdf_ptunc.at(0), pdf_ptunc.at(1)));

    RooRelBW2ProngPdf Vpdf("Vpdf", "Vpdf", Vpars, vars);
    RooProdPdf totalpdf("totalpdf", "totalpdf", RooArgList(Vpdf, constraintspdf));
    RooFormulaVar logL_rfv("logLikelihood_rfv", "-log(@0)", RooArgList(totalpdf));

    float mass_prefit = Vpdf.getQSq(); // Not important from which pdf you read it

    bool const dofits = (mass_prefit>0.f && sqrt(mass_prefit)/RooRelBW2ProngPdf::GeVunit>=40.f && sqrt(mass_prefit)/RooRelBW2ProngPdf::GeVunit<=210.f);

    bool isFitOK=false;
    int fitstatus=-999;
    if (dofits) doMinimization(logL_rfv, theta_ptunc, isFitOK, fitstatus);
    if (isFitOK){
      res.push_back(FitResultSummary());
      FitResultSummary& lastFit = res.back();

      lastFit.fitstatus = fitstatus;
      lastFit.pdfval = totalpdf.getVal();
      lastFit.objid = Vid;

      lastFit.constraintval_Vpdf = lastFit.constraintval_Xpdf = constraintspdf.getVal();
      lastFit.Vmass = lastFit.Xmass = sqrt(Vpdf.getQSq())*100.f;
      lastFit.Vmass_prefit = lastFit.Xmass_prefit = sqrt(mass_prefit)*100.f;

      for (size_t ijet=0; ijet<nprongs; ijet++){
        lastFit.jetnuisances.push_back((!doSingleTheta ? theta_ptunc.at(ijet) : theta_ptunc.back()).getVal());
        lastFit.jetindices.push_back(comb.at(ijet));
      }
      if (lastFit.jetindices.at(0)>lastFit.jetindices.at(1)){
        std::swap(lastFit.jetnuisances.at(0), lastFit.jetnuisances.at(1));
        std::swap(lastFit.jetindices.at(0), lastFit.jetindices.at(1));
      }
    }
  } // End loop over permutations
  std::sort(res.begin(), res.end(), [] (const FitResultSummary& a, const FitResultSummary& b){return a > b; });
}


void doFits_N_3_4(
  // Object inputs
  std::vector<std::pair<MELAParticle, float>> const& ak4jets, // The second member is the abs. pT resolution
  MELAParticle const& bestLeptonicZcand,

  // Outputs
  FitResultSummary& bestZfit,
  FitResultSummary& bestT1fit
){
  constexpr size_t nprongs = 3;
  unsigned int const nak4jets_preselected = ak4jets.size();
  if (nak4jets_preselected<3 || nak4jets_preselected>4) return;

  const bool doZWcomparison = (fabs(bestLeptonicZcand.m() - PDGHelpers::Zmass)>=5.f*PDGHelpers::Zwidth);

  std::vector<FitResultSummary> hadTopFits;
  getAllHadTopFits(ak4jets, hadTopFits);
  FitResultSummary const* preferredT1Fit = nullptr;
  if (!hadTopFits.empty() && hadTopFits.front().fitstatus>=0) preferredT1Fit = &(hadTopFits.front());

  // hadZFits and hadWFits are aligned because they require the same prefit cuts and the permutation lists are the same.
  std::vector<FitResultSummary> hadZFits;
  std::vector<FitResultSummary> hadWFits;
  if (doZWcomparison){
    getAllHadVFits(23, ak4jets, hadZFits);
    getAllHadVFits(24, ak4jets, hadWFits);
  }

  FitResultSummary const* preferredZFit = nullptr;
  for (FitResultSummary const& hadZFit:hadZFits){
    if (!preferredT1Fit) break;
    if (hadZFit.fitstatus<0) continue;

    bool isWLike=false;
    for (FitResultSummary const& hadWFit:hadWFits){
      bool allJetsSame=(hadZFit.jetindices.size()==hadWFit.jetindices.size());
      if (!allJetsSame) continue;
      for (unsigned int ijet=0; ijet<hadZFit.jetindices.size(); ijet++) allJetsSame &= hadZFit.jetindices.at(ijet)==hadWFit.jetindices.at(ijet);
      if (allJetsSame){
        isWLike = ((hadWFit.fitstatus==0 || hadWFit.fitstatus==4) && hadWFit.constraintval_Vpdf>hadZFit.constraintval_Vpdf);
        break;
      }
    }
    if (isWLike) continue;

    //if (fabs(hadZFit.Vmass-PDGHelpers::Zmass)>PDGHelpers::Zwidth) continue;
    if (
      !(hadZFit.jetindices.at(0)==preferredT1Fit->jetindices.at(1) && hadZFit.jetindices.at(1)==preferredT1Fit->jetindices.at(2))
      ||
      (
        hadZFit.jetindices.at(0)==preferredT1Fit->jetindices.at(1) && hadZFit.jetindices.at(1)==preferredT1Fit->jetindices.at(2)
        && (
          hadZFit.constraintval_Vpdf>preferredT1Fit->constraintval_Vpdf
          ||
          (hadZFit.constraintval_Vpdf==preferredT1Fit->constraintval_Vpdf && fabs(preferredT1Fit->Vmass-PDGHelpers::Wmass)>fabs(hadZFit.Vmass-PDGHelpers::Zmass))
          )
        )
      ){
      preferredZFit = &hadZFit;
      break;
    }
  }

  if (preferredT1Fit) bestT1fit = *preferredT1Fit;
  if (preferredZFit) bestZfit = *preferredZFit;
}

void matchRecoToGen(std::vector<MELAParticle> const& recolist, std::vector<MELAParticle*> const& genlist, std::vector<MELAParticle*>& reco_matchedGenParts){
  if (recolist.empty()) return;
  reco_matchedGenParts.clear();
  reco_matchedGenParts.assign(recolist.size(), nullptr);
  if (genlist.empty()) return;

  std::vector<unsigned int> remaining_recoparts; if (recolist.size()>0) remaining_recoparts.reserve(recolist.size());
  for (unsigned int ipart=0; ipart<recolist.size(); ipart++) remaining_recoparts.push_back(ipart);
  std::vector<unsigned int> remaining_genparts; remaining_genparts.reserve(genlist.size());
  for (unsigned int ipart=0; ipart<genlist.size(); ipart++) remaining_genparts.push_back(ipart);
  while (!remaining_recoparts.empty() && !remaining_genparts.empty()){
    int chosenRecoPart=-1;
    int chosenGenPart=-1;
    float minDeltaR=-1;
    for (unsigned int const& rpart:remaining_recoparts){
      TLorentzVector const& pReco = recolist.at(rpart).p4;
      for (unsigned int const& gpart:remaining_genparts){
        TLorentzVector const& pGen = genlist.at(gpart)->p4;
        float deltaR = pGen.DeltaR(pReco);
        if (minDeltaR==-1. || deltaR<minDeltaR){
          minDeltaR=deltaR;
          chosenRecoPart=rpart;
          chosenGenPart=gpart;
        }
      }
    }

    if (chosenRecoPart>=0 && chosenGenPart>=0) reco_matchedGenParts.at(chosenRecoPart) = genlist.at(chosenGenPart);
    for (auto it=remaining_recoparts.begin(); it!=remaining_recoparts.end(); it++){ if ((int) *it == chosenRecoPart){ remaining_recoparts.erase(it); break; } }
    for (auto it=remaining_genparts.begin(); it!=remaining_genparts.end(); it++){ if ((int) *it == chosenGenPart){ remaining_genparts.erase(it); break; } }
  }
}

void testAllFits(TString fname){
  gStyle->SetOptStat(0);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  SampleHelpers::theDataYear=2018;
  SampleHelpers::theDataPeriod="2018";
  SampleHelpers::theDataVersion=SampleHelpers::kCMSSW_10_X;
  float const btagWP = BtagHelpers::getBtagWP(BtagHelpers::kDeepCSV_Medium);
  float const totalLumi=150.f*1000.f;

  bool isWZ = (fname.Contains("WZTo"));
  bool isWZZ = (fname.Contains("WZZ"));
  bool isWWW = (fname.Contains("WWW"));
  bool isDYLLJets = (fname.Contains("DYJetsToLL"));
  bool isTTWLNu = (fname.Contains("TTWJetsToLNu"));
  bool isTTZLLNuNu = (fname.Contains("TTZToLLNuNu"));
  bool isTT2L2Nu = (fname.Contains("TTTo2L2Nu"));
  bool isTT2L2Nu_0j_special = (fname.Contains("TT_DiLept"));
  bool isTT2L2Nu_1j_special = (fname.Contains("TTPlus1Jet_DiLept"));
  bool isZZ4L = (fname.Contains("ZZTo4L"));
  if (!fname.Contains(".root")) fname += ".root";
  const TString strindir="output/test_GenWeights/2018/";
  TString stroutputcore=strindir+"/fitterTests/AllFits/";
  TString stroutput=stroutputcore;
  TString sampleLabel;
  if (isWZ) sampleLabel="WZTo3LNu_comparison";
  if (isWZZ) sampleLabel="WZZ_comparison";
  if (isWWW) sampleLabel="WWW_comparison";
  if (isDYLLJets) sampleLabel="DYJetsToLL_comparison";
  if (isTTWLNu) sampleLabel="TTWLNu_comparison";
  if (isTTZLLNuNu) sampleLabel="TTZLLNuNu_comparison";
  if (isTT2L2Nu) sampleLabel="TTTo2L2Nu_comparison";
  if (isTT2L2Nu_0j_special) sampleLabel="TTTo2L2Nu_0j_special_comparison";
  if (isTT2L2Nu_1j_special) sampleLabel="TTTo2L2Nu_1j_special_comparison";
  if (isZZ4L) sampleLabel="ZZTo4L_comparison";
  stroutput += sampleLabel + ".root";

  TFile* finput = TFile::Open(strindir+fname, "read");
  TTree* tree = (TTree*) finput->Get("test");
  tree->SetBranchStatus("*", 0);

  std::vector<float> genweights(WeightVariables::nWeightTypes, 0.f);
  std::vector<float> oneWeight(WeightVariables::nWeightTypes, 1.f);
  for (int iwt=WeightVariables::wCentral; iwt<WeightVariables::nWeightTypes; iwt++){
    TString bname = WeightVariables::getWeightName((WeightVariables::WeightType) iwt);
    tree->SetBranchStatus(bname, 1);
    tree->SetBranchAddress(bname, &(genweights.at(iwt)));
  }

  // XSEC
  float theXSec=-1;
  DATA_SIMPLE(float, xsec)

  // MET
  DATA_SIMPLE(float, pfmet)
  DATA_SIMPLE(float, pfmetPhi)
  // Electrons
  DATA_CLASS(std::vector<int>, electrons_id)
  DATA_CLASS(std::vector<float>, electrons_pt)
  DATA_CLASS(std::vector<float>, electrons_eta)
  DATA_CLASS(std::vector<float>, electrons_phi)
  DATA_CLASS(std::vector<float>, electrons_mass)
  DATA_CLASS(std::vector<long long>, electrons_selectionBits)
  // Muons
  DATA_CLASS(std::vector<int>, muons_id)
  DATA_CLASS(std::vector<float>, muons_pt)
  DATA_CLASS(std::vector<float>, muons_eta)
  DATA_CLASS(std::vector<float>, muons_phi)
  DATA_CLASS(std::vector<float>, muons_mass)
  DATA_CLASS(std::vector<long long>, muons_selectionBits)
  // AK4Jets
  DATA_CLASS(std::vector<float>, ak4jets_pt)
  DATA_CLASS(std::vector<float>, ak4jets_eta)
  DATA_CLASS(std::vector<float>, ak4jets_phi)
  DATA_CLASS(std::vector<float>, ak4jets_mass)
  DATA_CLASS(std::vector<long long>, ak4jets_selectionBits)
  DATA_CLASS(std::vector<float>, ak4jets_estimatedPtResolution)
  DATA_CLASS(std::vector<float>, ak4jets_deepCSVc)
  DATA_CLASS(std::vector<float>, ak4jets_deepCSVb)
  DATA_CLASS(std::vector<float>, ak4jets_deepCSVbb)
  // AK8Jets
  DATA_CLASS(std::vector<float>, ak8jets_pt)
  DATA_CLASS(std::vector<float>, ak8jets_eta)
  DATA_CLASS(std::vector<float>, ak8jets_phi)
  DATA_CLASS(std::vector<float>, ak8jets_mass)
  DATA_CLASS(std::vector<long long>, ak8jets_selectionBits)
  DATA_CLASS(std::vector<float>, ak8jets_estimatedPtResolution)
  DATA_CLASS(std::vector<float>, ak8jets_deepdisc_qcd)
  DATA_CLASS(std::vector<float>, ak8jets_deepdisc_top)
  DATA_CLASS(std::vector<float>, ak8jets_deepdisc_w)
  DATA_CLASS(std::vector<float>, ak8jets_deepdisc_z)
  DATA_CLASS(std::vector<float>, ak8jets_deepdisc_zbb)
  // TFTops
  DATA_CLASS(std::vector<float>, tftops_pt)
  DATA_CLASS(std::vector<float>, tftops_eta)
  DATA_CLASS(std::vector<float>, tftops_phi)
  DATA_CLASS(std::vector<float>, tftops_mass)
  DATA_CLASS(std::vector<long long>, tftops_selectionBits)
  DATA_CLASS(std::vector<int>, tftops_nSubjets)
  DATA_CLASS(std::vector<float>, tftops_disc)

  // Gen. particles
  //DATA_CLASS(std::vector<bool>, genparticles_isPromptFinalState)
  //DATA_CLASS(std::vector<bool>, genparticles_isPromptDecayed)
  //DATA_CLASS(std::vector<bool>, genparticles_isDirectPromptTauDecayProductFinalState)
  DATA_CLASS(std::vector<bool>, genparticles_isHardProcess)
  //DATA_CLASS(std::vector<bool>, genparticles_fromHardProcessFinalState)
  //DATA_CLASS(std::vector<bool>, genparticles_fromHardProcessDecayed)
  //DATA_CLASS(std::vector<bool>, genparticles_isDirectHardProcessTauDecayProductFinalState)
  //DATA_CLASS(std::vector<bool>, genparticles_fromHardProcessBeforeFSR)
  //DATA_CLASS(std::vector<bool>, genparticles_isLastCopy)
  //DATA_CLASS(std::vector<bool>, genparticles_isLastCopyBeforeFSR)
  DATA_CLASS(std::vector<int>, genparticles_id)
  DATA_CLASS(std::vector<int>, genparticles_status)
  DATA_CLASS(std::vector<float>, genparticles_px)
  DATA_CLASS(std::vector<float>, genparticles_py)
  DATA_CLASS(std::vector<float>, genparticles_pz)
  DATA_CLASS(std::vector<float>, genparticles_mass)


  TFile* foutput = TFile::Open(stroutput, "recreate");
  TTree tout("ZFitTest", "");
  tout.Branch("pfmet", &pfmet);
  tout.Branch("pfmetPhi", &pfmetPhi);
  BOOK_BRANCH(float, transversemass_lep_MET, tout);
  BOOK_BRANCH(unsigned int, nak4jets_preselected, tout);
  BOOK_BRANCH(std::vector<int>, ak4jets_preselected_progenitorV_Id, tout);
  BOOK_BRANCH(std::vector<int>, ak4jets_preselected_matchedParton_Id, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_matchedParton_DeltaR, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_pt, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_eta, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_phi, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_mass, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_relPtResolution, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_deepCSVc, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_deepCSVb, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_deepCSVbb, tout);
  BOOK_BRANCH(float, bestTFTop_mass, tout);
  BOOK_BRANCH(float, bestTFTop_disc, tout);
  BOOK_BRANCH(int, bestTFTop_bestGenTopParticle_id, tout);
  BOOK_BRANCH(float, bestTFTop_bestGenTopParticle_mass, tout);
  BOOK_BRANCH(float, bestTFTop_bestGenTopParticle_dR, tout);
  BOOK_BRANCH(unsigned int, nelectrons_preselected, tout);
  BOOK_BRANCH(std::vector<float>, electrons_preselected_pt, tout);
  BOOK_BRANCH(std::vector<float>, electrons_preselected_eta, tout);
  BOOK_BRANCH(std::vector<float>, electrons_preselected_phi, tout);
  BOOK_BRANCH(std::vector<float>, electrons_preselected_mass, tout);
  BOOK_BRANCH(unsigned int, nmuons_preselected, tout);
  BOOK_BRANCH(std::vector<float>, muons_preselected_pt, tout);
  BOOK_BRANCH(std::vector<float>, muons_preselected_eta, tout);
  BOOK_BRANCH(std::vector<float>, muons_preselected_phi, tout);
  BOOK_BRANCH(std::vector<float>, muons_preselected_mass, tout);

  BOOK_BRANCH(float, bestLeptonicZ_mass, tout);
  BOOK_BRANCH(float, bestLeptonicZ_id, tout); // 11 or 13
  BOOK_BRANCH(std::vector<unsigned int>, bestLeptonicZ_leptonIndices, tout);
  BOOK_BRANCH(float, bestLeptonicZ_bestGenZParticle_mass, tout);
  BOOK_BRANCH(float, bestLeptonicZ_bestGenZParticle_dR, tout);

  BOOK_BRANCH(std::vector<unsigned int>, hadVfit_jetIndices, tout);
  BOOK_BRANCH(std::vector<float>, hadVfit_jetnuisances, tout);
  BOOK_BRANCH(int, hadVfit_status, tout);
  BOOK_BRANCH(int, hadVfit_VId, tout);
  BOOK_BRANCH(float, hadVfit_mV_prefit, tout);
  BOOK_BRANCH(float, hadVfit_mV_postfit, tout);
  BOOK_BRANCH(int, hadVfit_bestGenZParticle_id, tout);
  BOOK_BRANCH(float, hadVfit_bestGenZParticle_mass, tout);
  BOOK_BRANCH(float, hadVfit_bestGenZParticle_dR, tout);
  BOOK_BRANCH(int, hadVfit_postfit_bestGenZParticle_id, tout);
  BOOK_BRANCH(float, hadVfit_postfit_bestGenZParticle_mass, tout);
  BOOK_BRANCH(float, hadVfit_postfit_bestGenZParticle_dR, tout);

  BOOK_BRANCH(std::vector<unsigned int>, hadT1fit_jetIndices, tout);
  BOOK_BRANCH(std::vector<float>, hadT1fit_jetnuisances, tout);
  BOOK_BRANCH(int, hadT1fit_status, tout);
  BOOK_BRANCH(int, hadT1fit_VId, tout);
  BOOK_BRANCH(float, hadT1fit_mW_prefit, tout);
  BOOK_BRANCH(float, hadT1fit_mW_postfit, tout);
  BOOK_BRANCH(float, hadT1fit_mtop_prefit, tout);
  BOOK_BRANCH(float, hadT1fit_mtop_postfit, tout);
  BOOK_BRANCH(float, bestHadT1fitMatchTFTop_mass, tout);
  BOOK_BRANCH(float, bestHadT1fitMatchTFTop_disc, tout);
  BOOK_BRANCH(float, bestHadT1fitMatchTFTop_dR, tout);
  BOOK_BRANCH(int, hadT1fit_bestGenTopParticle_id, tout);
  BOOK_BRANCH(float, hadT1fit_bestGenTopParticle_mass, tout);
  BOOK_BRANCH(float, hadT1fit_bestGenTopParticle_dR, tout);
  BOOK_BRANCH(int, hadT1fit_postfit_bestGenTopParticle_id, tout);
  BOOK_BRANCH(float, hadT1fit_postfit_bestGenTopParticle_mass, tout);
  BOOK_BRANCH(float, hadT1fit_postfit_bestGenTopParticle_dR, tout);

  BOOK_BRANCH(std::vector<unsigned int>, hadT2fit_jetIndices, tout);
  BOOK_BRANCH(std::vector<float>, hadT2fit_jetnuisances, tout);
  BOOK_BRANCH(int, hadT2fit_status, tout);
  BOOK_BRANCH(int, hadT2fit_VId, tout);
  BOOK_BRANCH(float, hadT2fit_mW_prefit, tout);
  BOOK_BRANCH(float, hadT2fit_mW_postfit, tout);
  BOOK_BRANCH(float, hadT2fit_mtop_prefit, tout);
  BOOK_BRANCH(float, hadT2fit_mtop_postfit, tout);
  BOOK_BRANCH(float, bestHadT2fitMatchTFTop_mass, tout);
  BOOK_BRANCH(float, bestHadT2fitMatchTFTop_disc, tout);
  BOOK_BRANCH(float, bestHadT2fitMatchTFTop_dR, tout);
  BOOK_BRANCH(int, hadT2fit_bestGenTopParticle_id, tout);
  BOOK_BRANCH(float, hadT2fit_bestGenTopParticle_mass, tout);
  BOOK_BRANCH(float, hadT2fit_bestGenTopParticle_dR, tout);
  BOOK_BRANCH(int, hadT2fit_postfit_bestGenTopParticle_id, tout);
  BOOK_BRANCH(float, hadT2fit_postfit_bestGenTopParticle_mass, tout);
  BOOK_BRANCH(float, hadT2fit_postfit_bestGenTopParticle_dR, tout);


  // Plotted variables
  struct Variable{
    TString name;
    TString title;
    unsigned int nbins;
    float min;
    float max;
    float val;
    std::vector<float> weights;

    Variable() : nbins(0), min(0), max(0), val(0), weights(WeightVariables::nWeightTypes, 0.f){}
    Variable(Variable const& other) : name(other.name), title(other.title), nbins(other.nbins), min(other.min), max(other.max), val(other.val), weights(WeightVariables::nWeightTypes, 0.f){}
    Variable(TString name_, TString title_) : name(name_), title(title_), nbins(0), min(0), max(0), val(0), weights(WeightVariables::nWeightTypes, 0.f){}
    Variable(TString name_, TString title_, unsigned int nbins_, float min_, float max_) : name(name_), title(title_), nbins(nbins_), min(min_), max(max_), val(0), weights(WeightVariables::nWeightTypes, 0.f){}

    void setVal(float v, std::vector<float> const& wgts){
      float binwidth = (max-min)/float(nbins);
      val=v;
      if (val>=max) val=max-binwidth/2.;
      if (val<=min) val=min+binwidth/2.;
      weights=wgts;
    }
    void reset(){ val=min - fabs(min); weights=std::vector<float>(WeightVariables::nWeightTypes, 0.f); } // So that filling is done on the 0th bin
  };

  std::vector<Variable> plotvars;
  //plotvars.emplace_back("mll", "m_{ll} (GeV)", 40, 80, 100);
  plotvars.emplace_back("hadT1fit_mW_prefit", "m_{qq} (pre-fit) (GeV)", 60, 70, 100);
  plotvars.emplace_back("hadT1fit_mW_postfit", "m_{qq} (post-fit) (GeV)", 60, 70, 100);
  size_t const nvars = plotvars.size();

  std::vector<std::vector<TH1F*>> histColl(plotvars.size(), std::vector<TH1F*>(WeightVariables::nWeightTypes, nullptr));
  for (size_t icol=0; icol<histColl.size(); icol++){
    auto& var=plotvars.at(icol);
    auto& col=histColl.at(icol);
    for (int iwt=0; iwt<WeightVariables::nWeightTypes; iwt++){
      TH1F*& h=col.at(iwt);
      TString hname = var.name+"_"+WeightVariables::getWeightName((WeightVariables::WeightType) iwt);
      HelperFunctions::replaceString<TString, const TString>(hname, TString("_")+WeightVariables::getWeightName(WeightVariables::wCentral), "");
      h = new TH1F(hname, "", var.nbins, var.min, var.max);
      h->GetXaxis()->SetTitle(var.title);
      h->GetYaxis()->SetTitle("a.u.");
      h->Sumw2();
      h->SetLineColor(getHistogramColorByWeightType((WeightVariables::WeightType) iwt));
      h->SetLineStyle(getHistogramDashByWeightType((WeightVariables::WeightType) iwt));
      h->SetLineWidth(2);
    }
  }

  float sumCtrWgt[WeightVariables::nWeightTypes]={ 0 };
  int const nEntries = tree->GetEntries();
  for (int ev=0; ev<nEntries; ev++){
    tree->GetEntry(ev);
    HelperFunctions::progressbar(ev, nEntries);

    if (ev==0){ theXSec=xsec; tree->SetBranchStatus("xsec", 0); }

    for (int iwt=0; iwt<WeightVariables::nWeightTypes; iwt++) sumCtrWgt[iwt] += genweights[iwt];

    transversemass_lep_MET=-1;
    nak4jets_preselected=0;
    ak4jets_preselected_progenitorV_Id.clear();
    ak4jets_preselected_matchedParton_Id.clear();
    ak4jets_preselected_matchedParton_DeltaR.clear();
    ak4jets_preselected_pt.clear();
    ak4jets_preselected_eta.clear();
    ak4jets_preselected_phi.clear();
    ak4jets_preselected_mass.clear();
    ak4jets_preselected_relPtResolution.clear();
    ak4jets_preselected_deepCSVb.clear();
    ak4jets_preselected_deepCSVbb.clear();

    bestTFTop_mass=-1;
    bestTFTop_disc=-1;
    bestTFTop_bestGenTopParticle_id=-9000;
    bestTFTop_bestGenTopParticle_mass=-1;
    bestTFTop_bestGenTopParticle_dR=-1;
    nelectrons_preselected=0;
    electrons_preselected_pt.clear();
    electrons_preselected_eta.clear();
    electrons_preselected_phi.clear();
    electrons_preselected_mass.clear();
    nmuons_preselected=0;
    muons_preselected_pt.clear();
    muons_preselected_eta.clear();
    muons_preselected_phi.clear();
    muons_preselected_mass.clear();
    bestLeptonicZ_mass=-1;
    bestLeptonicZ_id=-9000; // 11 or 13
    bestLeptonicZ_leptonIndices.clear();
    bestLeptonicZ_bestGenZParticle_mass=-1;
    bestLeptonicZ_bestGenZParticle_dR=-1;

    hadVfit_jetIndices.clear();
    hadVfit_jetnuisances.clear();
    hadVfit_status=-1;
    hadVfit_VId=-9000;
    hadVfit_mV_prefit=-1;
    hadVfit_mV_postfit=-1;
    hadVfit_bestGenZParticle_id=-9000;
    hadVfit_bestGenZParticle_mass=-1;
    hadVfit_bestGenZParticle_dR=-1;
    hadVfit_postfit_bestGenZParticle_id=-9000;
    hadVfit_postfit_bestGenZParticle_mass=-1;
    hadVfit_postfit_bestGenZParticle_dR=-1;

    hadT1fit_jetIndices.clear();
    hadT1fit_jetnuisances.clear();
    hadT1fit_status=-1;
    hadT1fit_VId=-9000;
    hadT1fit_mW_prefit=-1;
    hadT1fit_mW_postfit=-1;
    hadT1fit_mtop_prefit=-1;
    hadT1fit_mtop_postfit=-1;
    bestHadT1fitMatchTFTop_mass=-1;
    bestHadT1fitMatchTFTop_disc=-1;
    bestHadT1fitMatchTFTop_dR=-1;
    hadT1fit_bestGenTopParticle_id=-9000;
    hadT1fit_bestGenTopParticle_mass=-1;
    hadT1fit_bestGenTopParticle_dR=-1;
    hadT1fit_postfit_bestGenTopParticle_id=-9000;
    hadT1fit_postfit_bestGenTopParticle_mass=-1;
    hadT1fit_postfit_bestGenTopParticle_dR=-1;

    hadT2fit_jetIndices.clear();
    hadT2fit_jetnuisances.clear();
    hadT2fit_status=-1;
    hadT2fit_VId=-9000;
    hadT2fit_mW_prefit=-1;
    hadT2fit_mW_postfit=-1;
    hadT2fit_mtop_prefit=-1;
    hadT2fit_mtop_postfit=-1;
    bestHadT2fitMatchTFTop_mass=-1;
    bestHadT2fitMatchTFTop_disc=-1;
    bestHadT2fitMatchTFTop_dR=-1;
    hadT2fit_bestGenTopParticle_id=-9000;
    hadT2fit_bestGenTopParticle_mass=-1;
    hadT2fit_bestGenTopParticle_dR=-1;
    hadT2fit_postfit_bestGenTopParticle_id=-9000;
    hadT2fit_postfit_bestGenTopParticle_mass=-1;
    hadT2fit_postfit_bestGenTopParticle_dR=-1;


    // Gen particles
    size_t const ngenparticles = genparticles_px->size();
    std::vector<MELAParticle> genparticles; genparticles.reserve(ngenparticles);
    std::vector<MELAParticle*> genZlhe, genWlhe, genToplhe, genHardParticles;
    for (size_t ip=0; ip<ngenparticles; ip++){
      const int& id = genparticles_id->at(ip);
      double E = pow(genparticles_px->at(ip), 2)+pow(genparticles_py->at(ip), 2)+pow(genparticles_pz->at(ip), 2);
      E += (genparticles_mass->at(ip)>=0.f ? 1.f : -1.f) * pow(genparticles_mass->at(ip), 2);
      E = sqrt(E);
      TLorentzVector v(genparticles_px->at(ip), genparticles_py->at(ip), genparticles_pz->at(ip), E);
      genparticles.emplace_back(id, v);
      genparticles.back().setGenStatus(genparticles_status->at(ip));
      if (genparticles_isHardProcess->at(ip)){
        genHardParticles.push_back(&(genparticles.back()));
        if (std::abs(id)==6) genToplhe.push_back(&(genparticles.back()));
        else if (std::abs(id)==23) genZlhe.push_back(&(genparticles.back()));
        else if (std::abs(id)==24) genWlhe.push_back(&(genparticles.back()));
      }
    }
    // Find the mother/daughter relationship for the Z, W and tops
    for (auto Wpart:genWlhe){
      int const& idW = Wpart->id;
      TLorentzVector const& pW = Wpart->p4;
      const float WmassSq = pW.M2();
      MELAParticle* pDau1=nullptr;
      MELAParticle* pDau2=nullptr;
      float bestMassDiff=-1;
      for (auto ipart=genHardParticles.begin(); ipart!=genHardParticles.end(); ipart++){
        if (!((*ipart)->genStatus==23 || (*ipart)->genStatus==1 || (*ipart)->genStatus==2)) continue;
        for (auto jpart=ipart+1; jpart!=genHardParticles.end(); jpart++){
          if (!((*jpart)->genStatus==23 || (*jpart)->genStatus==1 || (*jpart)->genStatus==2)) continue;

          if (idW!=PDGHelpers::getCoupledVertex((*ipart)->id, (*jpart)->id)) continue;
          TLorentzVector pTotal = (*ipart)->p4 + (*jpart)->p4;
          float const tmpMassDiff = fabs(pTotal.M2() - WmassSq);
          if (bestMassDiff<0.f || bestMassDiff>tmpMassDiff){
            bestMassDiff = tmpMassDiff;
            pDau1 = *ipart;
            pDau2 = *jpart;
          }
        }
      }
      if (pDau1 && pDau2){
        if (pDau1->id<0 && pDau2->id>0) std::swap(pDau1, pDau2);
        Wpart->addDaughter(pDau1);
        Wpart->addDaughter(pDau2);
        pDau1->addMother(Wpart);
        pDau2->addMother(Wpart);
        /*
        MELAout << "\t- W particle:" << endl;
        MELAout << *Wpart << endl;
        MELAout << "\t\t- Daughters:" << endl;
        MELAout << *pDau1 << endl;
        MELAout << *pDau2 << endl;
        */
      }
    }
    for (auto Zpart:genZlhe){
      int const& idZ = Zpart->id;
      TLorentzVector const& pZ = Zpart->p4;
      const float ZmassSq = pZ.M2();
      MELAParticle* pDau1=nullptr;
      MELAParticle* pDau2=nullptr;
      float bestMassDiff=-1;
      for (auto ipart=genHardParticles.begin(); ipart!=genHardParticles.end(); ipart++){
        if (!((*ipart)->genStatus==23 || (*ipart)->genStatus==1 || (*ipart)->genStatus==2)) continue;
        for (auto jpart=ipart+1; jpart!=genHardParticles.end(); jpart++){
          if (!((*jpart)->genStatus==23 || (*jpart)->genStatus==1 || (*jpart)->genStatus==2)) continue;

          if (idZ!=PDGHelpers::getCoupledVertex((*ipart)->id, (*jpart)->id)) continue;
          TLorentzVector pTotal = (*ipart)->p4 + (*jpart)->p4;
          float const tmpMassDiff = fabs(pTotal.M2() - ZmassSq);
          if (bestMassDiff<0.f || bestMassDiff>tmpMassDiff){
            bestMassDiff = tmpMassDiff;
            pDau1 = *ipart;
            pDau2 = *jpart;
          }
        }
      }
      if (pDau1 && pDau2){
        if (pDau1->id<0 && pDau2->id>0) std::swap(pDau1, pDau2);
        Zpart->addDaughter(pDau1);
        Zpart->addDaughter(pDau2);
        pDau1->addMother(Zpart);
        pDau2->addMother(Zpart);
        /*
        MELAout << "\t- Z particle:" << endl;
        MELAout << *Zpart << endl;
        MELAout << "\t\t- Daughters:" << endl;
        MELAout << *pDau1 << endl;
        MELAout << *pDau2 << endl;
        */
      }
    }
    for (auto Toppart:genToplhe){
      int const& idTop = Toppart->id;
      int const bsign = (idTop>0 ? 1 : -1);
      int const Widmatch = 24*bsign;

      TLorentzVector const& pTop = Toppart->p4;
      const float TopmassSq = pTop.M2();
      MELAParticle* pDau1=nullptr;
      MELAParticle* pDau2=nullptr;
      MELAParticle* pDau3=nullptr;
      float bestMassDiff=-1;
      for (auto ipart=genHardParticles.begin(); ipart!=genHardParticles.end(); ipart++){
        if (!((*ipart)->genStatus==23 || (*ipart)->genStatus==1 || (*ipart)->genStatus==2)) continue;
        for (auto jpart=ipart+1; jpart!=genHardParticles.end(); jpart++){
          if (!((*jpart)->genStatus==23 || (*jpart)->genStatus==1 || (*jpart)->genStatus==2)) continue;
          for (auto kpart=jpart+1; kpart!=genHardParticles.end(); kpart++){
            if (!((*kpart)->genStatus==23 || (*kpart)->genStatus==1 || (*kpart)->genStatus==2)) continue;

            int const id12 = PDGHelpers::getCoupledVertex((*ipart)->id, (*jpart)->id);
            int const id13 = PDGHelpers::getCoupledVertex((*ipart)->id, (*kpart)->id);
            int const id23 = PDGHelpers::getCoupledVertex((*jpart)->id, (*kpart)->id);
            if (!((id12==Widmatch && (*kpart)->id*bsign>0) || (id13==Widmatch && (*jpart)->id*bsign>0) || (id23==Widmatch && (*ipart)->id*bsign>0))) continue;
            TLorentzVector pTotal = (*ipart)->p4 + (*jpart)->p4 + (*kpart)->p4;
            float const tmpMassDiff = fabs(pTotal.M2() - TopmassSq);
            if (bestMassDiff<0.f || bestMassDiff>tmpMassDiff){
              bestMassDiff = tmpMassDiff;
              pDau1 = *ipart;
              pDau2 = *jpart;
              pDau3 = *kpart;
            }
          }
        }
      }
      if (pDau1 && pDau2 && pDau3){
        /*
        MELAout << "\t- Top particle:" << endl;
        MELAout << *Toppart << endl;
        MELAout << "\t\t- Potential daughters:" << endl;
        MELAout << *pDau1 << endl;
        MELAout << *pDau2 << endl;
        MELAout << *pDau3 << endl;
        */
        int const id12 = PDGHelpers::getCoupledVertex(pDau1->id, pDau2->id);
        int const id13 = PDGHelpers::getCoupledVertex(pDau1->id, pDau3->id);
        int const id23 = PDGHelpers::getCoupledVertex(pDau2->id, pDau3->id);
        std::pair<MELAParticle*, MELAParticle*> bestWpair(nullptr, nullptr);
        MELAParticle* bestB=nullptr;
        float bestWpairMassDiff=-1;
        if (id12==Widmatch && PDGHelpers::isAQuark(pDau3->id) && pDau3->id*bsign>0 && (bestWpairMassDiff<0.f || bestWpairMassDiff>fabs((pDau1->p4+pDau2->p4).M()-PDGHelpers::Wmass))){
          bestWpairMassDiff = fabs((pDau1->p4+pDau2->p4).M()-PDGHelpers::Wmass);
          bestB=pDau3;
          bestWpair.first=pDau1;
          bestWpair.second=pDau2;
        }
        if (id13==Widmatch && PDGHelpers::isAQuark(pDau2->id) && pDau2->id*bsign>0 && (bestWpairMassDiff<0.f || bestWpairMassDiff>fabs((pDau1->p4+pDau3->p4).M()-PDGHelpers::Wmass))){
          bestWpairMassDiff = fabs((pDau1->p4+pDau3->p4).M()-PDGHelpers::Wmass);
          bestB=pDau2;
          bestWpair.first=pDau1;
          bestWpair.second=pDau3;
        }
        if (id23==Widmatch && PDGHelpers::isAQuark(pDau1->id) && pDau1->id*bsign>0 && (bestWpairMassDiff<0.f || bestWpairMassDiff>fabs((pDau2->p4+pDau3->p4).M()-PDGHelpers::Wmass))){
          bestWpairMassDiff = fabs((pDau2->p4+pDau3->p4).M()-PDGHelpers::Wmass);
          bestB=pDau1;
          bestWpair.first=pDau2;
          bestWpair.second=pDau3;
        }
        if (bestWpair.first->id<0 && bestWpair.second->id>0) std::swap(bestWpair.first, bestWpair.second);
        Toppart->addDaughter(bestB);
        Toppart->addDaughter(bestWpair.first);
        Toppart->addDaughter(bestWpair.second);
        if (!genWlhe.empty()){
          for (auto& genW:genWlhe){
            if (genW->hasDaughter(bestWpair.first) && genW->hasDaughter(bestWpair.second)){
              genW->addMother(Toppart);
              break;
            }
          }
        }
      }
      if (Toppart->getNDaughters()!=3){
        MELAerr << "LHE top particle does not have daughters! Summary of the hard process particles:" << endl;
        for (auto part:genHardParticles) MELAerr << *part << " | status: " << part->genStatus << endl;
        exit(1);
      }
    }

    // Reco particles
    MELAParticle bestLeptonicZcand;
    bool hasLeptonicZ=false;

    size_t const nelectrons = electrons_pt->size();
    std::vector<MELAParticle> electrons; electrons.reserve(nelectrons);
    for (size_t ip=0; ip<nelectrons; ip++){
      if (!HelperFunctions::test_bit(electrons_selectionBits->at(ip), ElectronSelectionHelpers::kMediumIDReco)) continue;
      TLorentzVector v; v.SetPtEtaPhiM(electrons_pt->at(ip), electrons_eta->at(ip), electrons_phi->at(ip), electrons_mass->at(ip));
      if (v.Pt()<ElectronSelectionHelpers::ptThr_skim_medium || fabs(v.Eta())>=ElectronSelectionHelpers::etaThr_skim_loose) continue; // Preselection requires some crazy tight eta cut, so apply manually
      electrons.emplace_back(electrons_id->at(ip), v);
    }
    nelectrons_preselected = electrons.size();
    if (getBestZCandidate(electrons, bestLeptonicZcand)){
      bestLeptonicZ_leptonIndices.clear();
      for (unsigned int ip=0; ip<electrons.size();ip++){
        for (unsigned int jp=0; jp<electrons.size(); jp++){
          if (ip==jp) continue;
          if (bestLeptonicZcand.getDaughter(0)==&(electrons.at(ip)) && bestLeptonicZcand.getDaughter(1)==&(electrons.at(jp))){
            bestLeptonicZ_leptonIndices.push_back(ip);
            bestLeptonicZ_leptonIndices.push_back(jp);
            bestLeptonicZ_id = 11;
            bestLeptonicZ_mass = bestLeptonicZcand.m();
          }
        }
      }
      hasLeptonicZ = true;
    }

    size_t const nmuons = muons_pt->size();
    std::vector<MELAParticle> muons; muons.reserve(nmuons);
    for (size_t ip=0; ip<nmuons; ip++){
      if (!HelperFunctions::test_bit(muons_selectionBits->at(ip), MuonSelectionHelpers::kPreselection)) continue;
      TLorentzVector v; v.SetPtEtaPhiM(muons_pt->at(ip), muons_eta->at(ip), muons_phi->at(ip), muons_mass->at(ip));
      muons.emplace_back(muons_id->at(ip), v);
    }
    nmuons_preselected = muons.size();
    if (getBestZCandidate(muons, bestLeptonicZcand)){
      bestLeptonicZ_leptonIndices.clear();
      for (unsigned int ip=0; ip<muons.size(); ip++){
        for (unsigned int jp=0; jp<muons.size(); jp++){
          if (ip==jp) continue;
          if (bestLeptonicZcand.getDaughter(0)==&(muons.at(ip)) && bestLeptonicZcand.getDaughter(1)==&(muons.at(jp))){
            bestLeptonicZ_leptonIndices.push_back(ip);
            bestLeptonicZ_leptonIndices.push_back(jp);
            bestLeptonicZ_id = 13;
            bestLeptonicZ_mass = bestLeptonicZcand.m();
          }
        }
      }
      hasLeptonicZ = true;
    }
    // Match leptonic Z to gen. Z
    if (hasLeptonicZ){
      float mindr=-1;
      int igen=-1;
      for (size_t ip=0; ip<ngenparticles; ip++){
        auto const& genparticle = genparticles.at(ip);
        if (std::abs(genparticle.id)!=23) continue;
        float tmpdr = bestLeptonicZcand.deltaR(genparticle.p4);
        if (mindr<0.f || tmpdr<mindr){
          mindr = tmpdr;
          igen = ip;
        }
      }
      if (igen>=0){
        bestLeptonicZ_bestGenZParticle_mass = genparticles.at(igen).m();
        bestLeptonicZ_bestGenZParticle_dR = mindr;
      }
    }

    const bool hasZLL = (fabs(bestLeptonicZcand.m() - PDGHelpers::Zmass)>=5.f*PDGHelpers::Zwidth);
    MELAParticle const* leadingNonZLepton = nullptr;
    for (auto const& part:electrons){
      if (hasZLL && bestLeptonicZcand.hasDaughter(&part)) continue;
      if (!leadingNonZLepton || leadingNonZLepton->pt()<part.pt()) leadingNonZLepton = &part;
    }
    for (auto const& part:muons){
      if (hasZLL && bestLeptonicZcand.hasDaughter(&part)) continue;
      if (!leadingNonZLepton || leadingNonZLepton->pt()<part.pt()) leadingNonZLepton = &part;
    }
    if (leadingNonZLepton) transversemass_lep_MET = sqrt(2.f*leadingNonZLepton->pt()*pfmet*(1.f - cos(leadingNonZLepton->phi()-pfmetPhi)));

    size_t const nak4jets = ak4jets_pt->size();
    std::vector<MELAParticle> ak4jets; ak4jets.reserve(nak4jets);
    std::vector<std::pair<MELAParticle, float>> ak4jet_ptres_pairs; ak4jet_ptres_pairs.reserve(nak4jets);
    for (size_t ip=0; ip<nak4jets; ip++){
      if (!HelperFunctions::test_bit(ak4jets_selectionBits->at(ip), AK4JetSelectionHelpers::kTightID) || ak4jets_pt->at(ip)<AK4JetSelectionHelpers::ptThr_skim_preselection || fabs(ak4jets_eta->at(ip))>=4.7) continue;
      TLorentzVector v; v.SetPtEtaPhiM(ak4jets_pt->at(ip), ak4jets_eta->at(ip), ak4jets_phi->at(ip), ak4jets_mass->at(ip));
      ak4jets.emplace_back(0, v);
      ak4jet_ptres_pairs.emplace_back(ak4jets.back(), ak4jets_estimatedPtResolution->at(ip));

      ak4jets_preselected_pt.push_back(v.Pt());
      ak4jets_preselected_eta.push_back(v.Eta());
      ak4jets_preselected_phi.push_back(v.Phi());
      ak4jets_preselected_mass.push_back(v.M());
      ak4jets_preselected_relPtResolution.push_back(ak4jets_estimatedPtResolution->at(ip)/v.Pt());
      ak4jets_preselected_deepCSVc.push_back(ak4jets_deepCSVc->at(ip));
      ak4jets_preselected_deepCSVb.push_back(ak4jets_deepCSVb->at(ip));
      ak4jets_preselected_deepCSVbb.push_back(ak4jets_deepCSVbb->at(ip));
    }
    nak4jets_preselected = ak4jets.size();

    std::vector<MELAParticle*> ak4jets_preselected_bestGenParticles;
    matchRecoToGen(ak4jets, genHardParticles, ak4jets_preselected_bestGenParticles);
    for (size_t ip=0; ip<nak4jets_preselected; ip++){
      auto const& part=ak4jets_preselected_bestGenParticles.at(ip);
      if (!part){
        ak4jets_preselected_progenitorV_Id.push_back(-1);
        ak4jets_preselected_matchedParton_Id.push_back(-9000);
        ak4jets_preselected_matchedParton_DeltaR.push_back(-1);
      }
      else{
        ak4jets_preselected_matchedParton_Id.push_back(part->id);
        ak4jets_preselected_matchedParton_DeltaR.push_back(part->deltaR(ak4jets.at(ip)));
        if (part->getNMothers()>0){
          ak4jets_preselected_progenitorV_Id.push_back(part->getMother(0)->id);
        }
        else{
          ak4jets_preselected_progenitorV_Id.push_back(-9000);
        }
      }
    }

    size_t const nak8jets = ak8jets_pt->size();
    std::vector<MELAParticle> ak8jets; ak8jets.reserve(nak8jets);
    for (size_t ip=0; ip<nak8jets; ip++){
      TLorentzVector v; v.SetPtEtaPhiM(ak8jets_pt->at(ip), ak8jets_eta->at(ip), ak8jets_phi->at(ip), ak8jets_mass->at(ip));
      ak8jets.emplace_back(0, v);
    }

    size_t const ntftops = tftops_pt->size();
    std::vector<MELAParticle> tftops; tftops.reserve(ntftops);
    MELAParticle* bestTFTop = nullptr;
    for (size_t ip=0; ip<ntftops; ip++){
      TLorentzVector v; v.SetPtEtaPhiM(tftops_pt->at(ip), tftops_eta->at(ip), tftops_phi->at(ip), tftops_mass->at(ip));
      tftops.emplace_back(0, v);
      if (bestTFTop_disc<0.f || tftops_disc->at(ip)>bestTFTop_disc){
        bestTFTop_disc = tftops_disc->at(ip);
        bestTFTop_mass = v.M();
        bestTFTop = &(tftops.back());
      }
    }
    if (bestTFTop){
      float mindr=-1;
      int igen=-1;
      for (size_t ip=0; ip<ngenparticles; ip++){
        auto const& genparticle = genparticles.at(ip);
        if (std::abs(genparticle.id)!=6) continue;
        float tmpdr = bestTFTop->deltaR(genparticle.p4);
        if (mindr<0.f || tmpdr<mindr){
          mindr = tmpdr;
          igen = ip;
        }
      }
      if (igen>=0){
        bestTFTop_bestGenTopParticle_id = genparticles.at(igen).id;
        bestTFTop_bestGenTopParticle_mass = genparticles.at(igen).m();
        bestTFTop_bestGenTopParticle_dR = mindr;
      }
    }


    FitResultSummary bestZfit;
    FitResultSummary bestT1fit;
    FitResultSummary bestT2fit;

    // Do the different fits
    //doFits_N_3_4(
    //  ak4jet_ptres_pairs, bestLeptonicZcand,
    //  bestZfit, bestT1fit
    //);

    hadVfit_jetIndices = bestZfit.jetindices;
    hadVfit_jetnuisances = bestZfit.jetnuisances;
    hadVfit_status = bestZfit.fitstatus;
    hadVfit_VId = bestZfit.objid;
    hadVfit_mV_prefit = bestZfit.Xmass_prefit;
    hadVfit_mV_postfit = bestZfit.Xmass;
    if (!hadVfit_jetIndices.empty()){
      TLorentzVector bestHadVfit_p4(0, 0, 0, 0);
      TLorentzVector bestHadVfit_postfit_p4(0, 0, 0, 0);
      for (unsigned int jin=0; jin<hadT1fit_jetIndices.size(); jin++){
        auto const& ak4jet = ak4jets.at(hadT1fit_jetIndices.at(jin));
        bestHadVfit_p4 += ak4jet.p4;
        bestHadVfit_postfit_p4 += ak4jet.p4 * pow(ak4jet_ptres_pairs.at(hadT1fit_jetIndices.at(jin)).second/ak4jet.pt(), hadT1fit_jetnuisances.at(jin));
      }

      // Match hadronic Z to gen. Z or W
      {
        float mindr=-1;
        int igen=-1;
        for (size_t ip=0; ip<ngenparticles; ip++){
          auto const& genparticle = genparticles.at(ip);
          if (std::abs(genparticle.id)!=23 && std::abs(genparticle.id)!=24) continue;
          float tmpdr = bestHadVfit_p4.DeltaR(genparticle.p4);
          if (mindr<0.f || tmpdr<mindr){
            mindr = tmpdr;
            igen = ip;
          }
        }
        if (igen>=0){
          hadVfit_bestGenZParticle_id = genparticles.at(igen).id;
          hadVfit_bestGenZParticle_mass = genparticles.at(igen).m();
          hadVfit_bestGenZParticle_dR = mindr;
        }
      }
      {
        float mindr=-1;
        int igen=-1;
        for (size_t ip=0; ip<ngenparticles; ip++){
          auto const& genparticle = genparticles.at(ip);
          if (std::abs(genparticle.id)!=23 && std::abs(genparticle.id)!=24) continue;
          float tmpdr = bestHadVfit_postfit_p4.DeltaR(genparticle.p4);
          if (mindr<0.f || tmpdr<mindr){
            mindr = tmpdr;
            igen = ip;
          }
        }
        if (igen>=0){
          hadVfit_postfit_bestGenZParticle_id = genparticles.at(igen).id;
          hadVfit_postfit_bestGenZParticle_mass = genparticles.at(igen).m();
          hadVfit_postfit_bestGenZParticle_dR = mindr;
        }
      }
    }

    hadT1fit_jetIndices = bestT1fit.jetindices;
    hadT1fit_jetnuisances = bestT1fit.jetnuisances;
    hadT1fit_status = bestT1fit.fitstatus;
    hadT1fit_VId = bestT1fit.objid;
    hadT1fit_mW_prefit = bestT1fit.Vmass_prefit;
    hadT1fit_mW_postfit = bestT1fit.Vmass;
    hadT1fit_mtop_prefit = bestT1fit.Xmass_prefit;
    hadT1fit_mtop_postfit = bestT1fit.Xmass;
    if (!hadT1fit_jetIndices.empty()){
      TLorentzVector bestT1fit_p4(0, 0, 0, 0);
      TLorentzVector bestT1fit_postfit_p4(0, 0, 0, 0);
      for (unsigned int jin=0; jin<hadT1fit_jetIndices.size();jin++){
        auto const& ak4jet = ak4jets.at(hadT1fit_jetIndices.at(jin));
        bestT1fit_p4 += ak4jet.p4;
        bestT1fit_postfit_p4 += ak4jet.p4 * pow(ak4jet_ptres_pairs.at(hadT1fit_jetIndices.at(jin)).second/ak4jet.pt(), hadT1fit_jetnuisances.at(jin));
      }
      for (unsigned int itop=0; itop<tftops.size(); itop++){
        MELAParticle const& tftop = tftops.at(itop);
        float tmpdr = tftop.deltaR(bestT1fit_p4);
        if (bestHadT1fitMatchTFTop_dR<0.f || tmpdr<bestHadT1fitMatchTFTop_dR){
          bestHadT1fitMatchTFTop_mass=tftop.m();
          bestHadT1fitMatchTFTop_disc=tftops_disc->at(itop);
          bestHadT1fitMatchTFTop_dR=tmpdr;
        }
      }
      {
        float mindr=-1;
        int igen=-1;
        for (size_t ip=0; ip<ngenparticles; ip++){
          auto const& genparticle = genparticles.at(ip);
          if (std::abs(genparticle.id)!=6) continue;
          float tmpdr = bestT1fit_p4.DeltaR(genparticle.p4);
          if (mindr<0.f || tmpdr<mindr){
            mindr = tmpdr;
            igen = ip;
          }
        }
        if (igen>=0){
          hadT1fit_bestGenTopParticle_id = genparticles.at(igen).id;
          hadT1fit_bestGenTopParticle_mass = genparticles.at(igen).m();
          hadT1fit_bestGenTopParticle_dR = mindr;
        }
      }
      {
        float mindr=-1;
        int igen=-1;
        for (size_t ip=0; ip<ngenparticles; ip++){
          auto const& genparticle = genparticles.at(ip);
          if (std::abs(genparticle.id)!=6) continue;
          float tmpdr = bestT1fit_postfit_p4.DeltaR(genparticle.p4);
          if (mindr<0.f || tmpdr<mindr){
            mindr = tmpdr;
            igen = ip;
          }
        }
        if (igen>=0){
          hadT1fit_postfit_bestGenTopParticle_id = genparticles.at(igen).id;
          hadT1fit_postfit_bestGenTopParticle_mass = genparticles.at(igen).m();
          hadT1fit_postfit_bestGenTopParticle_dR = mindr;
        }
      }
    }

    // Fill the output tree
    tout.Fill();

    // Fill the output histograms
    for (size_t ivar=0; ivar<nvars; ivar++){
      auto& var=plotvars.at(ivar);
      var.reset();

      if (var.name.Contains("metcut") && pfmet<150.f) continue;

      if (var.name.BeginsWith("recomet")) var.setVal(pfmet, genweights);
      else if (var.name.BeginsWith("hadT1fit_mW_postfit")){ if (hadT1fit_mW_postfit<0.f){ continue; } var.setVal(hadT1fit_mW_postfit, genweights); }
      else if (var.name.BeginsWith("hadT1fit_mW_prefit")){ if (hadT1fit_mW_prefit<0.f){ continue; } var.setVal(hadT1fit_mW_prefit, genweights); }
      else if (var.name.BeginsWith("nak4jets")){
        size_t nc=0;
        for (size_t ip=0; ip<nak4jets; ip++){
          auto const& jet=ak4jets.at(ip);
          if (!HelperFunctions::test_bit(ak4jets_selectionBits->at(ip), AK4JetSelectionHelpers::kTightID) || jet.pt()<30. || fabs(jet.eta())>=4.7) continue;
          bool isBTagged = (!var.name.Contains("btagged") || ((ak4jets_deepCSVb->at(ip)+ak4jets_deepCSVbb->at(ip))>btagWP));
          if (isBTagged) nc++;
        }
        var.setVal(nc, genweights);
      }
      else if (var.name.BeginsWith("leadingptak4jet")){
        float v=-99;
        for (size_t ip=0; ip<nak4jets; ip++){
          auto const& jet=ak4jets.at(ip);
          if (!HelperFunctions::test_bit(ak4jets_selectionBits->at(ip), AK4JetSelectionHelpers::kTightID) || fabs(jet.eta())>=4.7) continue;
          bool isBTagged = (!var.name.Contains("btagged") || ((ak4jets_deepCSVb->at(ip)+ak4jets_deepCSVbb->at(ip))>btagWP));
          if (isBTagged){
            if (var.name.BeginsWith("leadingptak4jet_pt")) v=jet.pt();
            break;
          }
        }
        if (v==-99.f) continue;
        var.setVal(v, genweights);
      }
      else if (var.name.BeginsWith("ntftops")){
        size_t nc=ntftops;
        var.setVal(nc, genweights);
      }
      else if (var.name.BeginsWith("highestdisctftop")){
        float v=-99;
        int itop=-1;
        for (size_t ip=0; ip<ntftops; ip++){
          auto const& top=tftops.at(ip);
          if (itop<0 || tftops_disc->at(ip)>v){
            v=tftops_disc->at(ip);
            itop=ip;
          }
        }
        if (itop>=0){
          auto const& top=tftops.at(itop);
          if (var.name.BeginsWith("highestdisctftop_mass")) v=top.m();
        }
        if (v==-99.f) continue;
        var.setVal(v, genweights);
      }
      else if (var.name.BeginsWith("nak8jets")){
        size_t nc=nak8jets;
        var.setVal(nc, genweights);
      }
      else if (var.name.BeginsWith("leadingptak8jet")){
        float v=-99;
        for (size_t ip=0; ip<nak8jets; ip++){
          auto const& jet=ak8jets.at(ip);
          //if (fabs(jet.eta())>=4.7) continue;
          if (var.name.BeginsWith("leadingptak8jet_pt")) v=jet.pt();
          break;
        }
        if (v==-99.f) continue;
        var.setVal(v, genweights);
      }
      else if (var.name.BeginsWith("highesttopdiscak8jet")){
        float v=-99;
        int ijet=-1;
        for (size_t ip=0; ip<nak8jets; ip++){
          auto const& jet=ak8jets.at(ip);
          //if (fabs(jet.eta())>=4.7) continue;
          if (ijet<0 || ak8jets_deepdisc_top->at(ip)>v){ v=ak8jets_deepdisc_top->at(ip); ijet=ip; }
          break;
        }
        if (ijet>=0){
          auto const& jet=ak8jets.at(ijet);
          if (var.name.BeginsWith("highesttopdiscak8jet_mass")) v=jet.m();
          //else if (var.name.BeginsWith("highesttopdiscak8jet_disc")) v=ak8jets_deepdisc_top->at(ip); // v is already set
        }
        if (v==-99.f) continue;
        var.setVal(v, genweights);
      }
      else if (var.name.BeginsWith("nelectrons")){
        size_t nc=0;
        for (size_t ip=0; ip<nelectrons; ip++){
          auto const& part=electrons.at(ip);
          if (!HelperFunctions::test_bit(electrons_selectionBits->at(ip), ElectronSelectionHelpers::kMediumIDReco) || part.pt()<20. || fabs(part.eta())>=2.5) continue;
          nc++;
        }
        var.setVal(nc, genweights);
      }
      else if (var.name.BeginsWith("nmuons")){
        size_t nc=0;
        for (size_t ip=0; ip<nmuons; ip++){
          auto const& part=muons.at(ip);
          if (!HelperFunctions::test_bit(muons_selectionBits->at(ip), MuonSelectionHelpers::kMediumIDReco) || part.pt()<20. || fabs(part.eta())>=2.4) continue;
          nc++;
        }
        var.setVal(nc, genweights);
      }

      for (int iwt=0; iwt<WeightVariables::nWeightTypes; iwt++) histColl.at(ivar).at(iwt)->Fill(var.val, var.weights.at(iwt));
    }
  }

  foutput->WriteTObject(&tout); // Write the output tree

  for (auto& col:histColl){
    float ymin = std::numeric_limits<float>::max();
    float ymax = std::numeric_limits<float>::min();
    for (auto& h:col){
      if (sumCtrWgt[WeightVariables::wCentral]>0.f) h->Scale(1.f/sumCtrWgt[WeightVariables::wCentral]);
      for (int ix=1; ix<=h->GetNbinsX(); ix++){
        float bc = h->GetBinContent(ix);
        float be = h->GetBinError(ix);
        if (h==col.at(0) && (theXSec*totalLumi)>0.f){
          be = sqrt(fabs(bc)/(theXSec*totalLumi));
          h->SetBinError(ix, be);
        }
        ymin=std::min(bc-be, ymin);
        ymax=std::max(bc+be, ymax);
      }
    }
    for (auto& h:col){
      TString hname = h->GetName();
      float scaleFacUp=1.2;
      if (hname.Contains("y") || hname.Contains("leadingetaj_eta")) scaleFacUp=1.5;
      h->GetYaxis()->SetRangeUser(ymin*(ymin>=0.f ? 0.8 : 1.2), ymax*scaleFacUp);
      foutput->WriteTObject(h);
    }

    for (unsigned int ptype=0; ptype<2; ptype++){
      TString cext = "";
      std::vector<TH1F*> hplottable; hplottable.reserve(col.size());
      float ymin_ptype = std::numeric_limits<float>::max();
      float ymax_ptype = std::numeric_limits<float>::min();
      for (int iwt=0; iwt<WeightVariables::nWeightTypes; iwt++){
        TH1F* h=col.at(iwt);
        hplottable.push_back((TH1F*) h->Clone(Form("%s_ptype%i", h->GetName(), ptype)));
        if (ptype==1) hplottable.back()->Scale(sumCtrWgt[WeightVariables::wCentral]/sumCtrWgt[iwt]);
        if (ptype!=0){
          for (int ix=1; ix<=hplottable.back()->GetNbinsX(); ix++){
            float bc = hplottable.back()->GetBinContent(ix);
            float be = hplottable.back()->GetBinError(ix);
            ymin_ptype=std::min(bc-be, ymin_ptype);
            ymax_ptype=std::max(bc+be, ymax_ptype);
          }
        }
      }
      if (ptype!=0){
        for (auto& h:hplottable){
          TString hname = h->GetName();
          float scaleFacUp=1.2;
          if (hname.Contains("y") || hname.Contains("leadingetaj_eta")) scaleFacUp=1.5;
          h->GetYaxis()->SetRangeUser(ymin_ptype*(ymin_ptype>=0.f ? 0.8 : 1.2), ymax_ptype*scaleFacUp);
        }
      }
      if (ptype==1) cext = "_XsecNorm";

      {
        // Plot the histograms
        TString canvasname = sampleLabel + "_" + col.at(0)->GetName() + "_LHEWeights" + cext;
        TCanvas canvas(canvasname, "", 8, 30, 800, 800);
        canvas.cd();
        gStyle->SetOptStat(0);
        canvas.SetFillColor(0);
        canvas.SetBorderMode(0);
        canvas.SetBorderSize(2);
        canvas.SetTickx(1);
        canvas.SetTicky(1);
        canvas.SetLeftMargin(0.17);
        canvas.SetRightMargin(0.05);
        canvas.SetTopMargin(0.07);
        canvas.SetBottomMargin(0.13);
        canvas.SetFrameFillStyle(0);
        canvas.SetFrameBorderMode(0);
        canvas.SetFrameFillStyle(0);
        canvas.SetFrameBorderMode(0);

        size_t nLegend=0;
        std::vector<bool> doPlot;
        std::vector<bool> addToLegend;
        for (size_t iwt=0; iwt<col.size(); iwt++){
          WeightVariables::WeightType type = (WeightVariables::WeightType) iwt;
          bool isPlottable=(iwt==0) || type==WeightVariables::wCentral_Default || !(
            iwt>(int) WeightVariables::wAsMZDn
            ||
            type==WeightVariables::wPSUp || type==WeightVariables::wPSDn
            ||
            type==WeightVariables::wPDFUp_Default || type==WeightVariables::wPDFDn_Default
            );
          int ds = getHistogramDashByWeightType((WeightVariables::WeightType) iwt);
          doPlot.push_back(isPlottable);
          if (ds!=7 && isPlottable){ addToLegend.push_back(true); nLegend++; }
          else addToLegend.push_back(false);
        }

        const float legend_minX = 0.50;
        const float legend_maxY = 0.90;
        const float legend_maxX = 0.90;
        float legend_minY = legend_maxY - 0.05f*float(nLegend);

        TLegend legend(legend_minX, legend_minY, legend_maxX, legend_maxY);
        legend.SetBorderSize(0);
        legend.SetTextFont(42);
        legend.SetTextSize(0.04);
        legend.SetLineColor(0);
        legend.SetLineStyle(0);
        legend.SetLineWidth(1);
        legend.SetFillColor(0);
        legend.SetFillStyle(0);

        for (unsigned int iwt=0; iwt<hplottable.size(); iwt++){
          TString strLabel=WeightVariables::getWeightLabel((WeightVariables::WeightType) iwt, true);
          if (strLabel.Contains(" up")) HelperFunctions::replaceString<TString, const TString>(strLabel, " up", "");
          if (strLabel.Contains(" dn")) HelperFunctions::replaceString<TString, const TString>(strLabel, " dn", "");
          if (addToLegend.at(iwt)) legend.AddEntry(hplottable.at(iwt), strLabel, "l");
          if (doPlot.at(iwt)){
            if (iwt==0) hplottable.at(iwt)->Draw("e1p");
            hplottable.at(iwt)->Draw("histsame");
          }
        }
        legend.Draw("same");

        //canvas.SetLogy();
        canvas.RedrawAxis();
        canvas.Modified();
        canvas.Update();
        //canvas.SaveAs(stroutputcore+canvasname+".png");
        canvas.SaveAs(stroutputcore+canvasname+".pdf");
        foutput->WriteTObject(&canvas);
        canvas.Close();
      }

      if (ptype!=1){
        // Plot the histograms
        TString canvasname = sampleLabel + "_" + col.at(0)->GetName() + "_PythiaWeights" + cext;
        TCanvas canvas(canvasname, "", 8, 30, 800, 800);
        canvas.cd();
        gStyle->SetOptStat(0);
        canvas.SetFillColor(0);
        canvas.SetBorderMode(0);
        canvas.SetBorderSize(2);
        canvas.SetTickx(1);
        canvas.SetTicky(1);
        canvas.SetLeftMargin(0.17);
        canvas.SetRightMargin(0.05);
        canvas.SetTopMargin(0.07);
        canvas.SetBottomMargin(0.13);
        canvas.SetFrameFillStyle(0);
        canvas.SetFrameBorderMode(0);
        canvas.SetFrameFillStyle(0);
        canvas.SetFrameBorderMode(0);

        size_t nLegend=0;
        std::vector<bool> doPlot;
        std::vector<bool> addToLegend;
        for (size_t iwt=0; iwt<col.size(); iwt++){
          WeightVariables::WeightType type = (WeightVariables::WeightType) iwt;
          bool isPlottable=(iwt==0) || !(
            iwt<=(int) WeightVariables::wAsMZDn
            ||
            type==WeightVariables::wCentral_Default
            ||
            type==WeightVariables::wPSUp || type==WeightVariables::wPSDn
            ||
            type==WeightVariables::wPDFUp_Default || type==WeightVariables::wPDFDn_Default
            );
          int ds = getHistogramDashByWeightType((WeightVariables::WeightType) iwt);
          doPlot.push_back(isPlottable);
          if (ds!=7 && isPlottable){ addToLegend.push_back(true); nLegend++; }
          else addToLegend.push_back(false);
        }

        const float legend_minX = 0.50;
        const float legend_maxY = 0.90;
        const float legend_maxX = 0.90;
        float legend_minY = legend_maxY - 0.05f*float(nLegend);

        TLegend legend(legend_minX, legend_minY, legend_maxX, legend_maxY);
        legend.SetBorderSize(0);
        legend.SetTextFont(42);
        legend.SetTextSize(0.04);
        legend.SetLineColor(0);
        legend.SetLineStyle(0);
        legend.SetLineWidth(1);
        legend.SetFillColor(0);
        legend.SetFillStyle(0);

        for (unsigned int iwt=0; iwt<hplottable.size(); iwt++){
          TString strLabel=WeightVariables::getWeightLabel((WeightVariables::WeightType) iwt, true);
          if (strLabel.Contains(" up")) HelperFunctions::replaceString<TString, const TString>(strLabel, " up", "");
          if (strLabel.Contains(" dn")) HelperFunctions::replaceString<TString, const TString>(strLabel, " dn", "");
          if (addToLegend.at(iwt)) legend.AddEntry(hplottable.at(iwt), strLabel, "l");
          if (doPlot.at(iwt)){
            if (iwt==0) hplottable.at(iwt)->Draw("e1p");
            hplottable.at(iwt)->Draw("histsame");
          }
        }
        legend.Draw("same");

        //canvas.SetLogy();
        canvas.RedrawAxis();
        canvas.Modified();
        canvas.Update();
        //canvas.SaveAs(stroutputcore+canvasname+".png");
        canvas.SaveAs(stroutputcore+canvasname+".pdf");
        foutput->WriteTObject(&canvas);
        canvas.Close();
      }

      for (auto*& h:hplottable) delete h;
    }
  }



  for (auto& col:histColl){ for (auto& h:col) delete h; }
  foutput->Close();
  finput->Close();
}
