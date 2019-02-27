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
#include "TChain.h"
#include "TStyle.h"
#include "TLegend.h"


#define DATA_CLASS(t, name) \
  t* name=nullptr; \
  tree->SetBranchStatus(#name, 1); \
  tree->SetBranchAddress(#name, &name);
#define DATA_SIMPLE(t, name, ...) \
  t name=0; \
  tree->SetBranchStatus(#name, 1); \
  tree->SetBranchAddress(#name, &name);
#define DATA_CLASS_CONDITIONAL(t, name, cond) \
  t* name=nullptr; \
  if (cond){ \
    tree->SetBranchStatus(#name, 1); \
    tree->SetBranchAddress(#name, &name); \
  }
#define DATA_SIMPLE_CONDITIONAL(t, name, cond) \
  t name=0; \
  if (cond){ \
    tree->SetBranchStatus(#name, 1); \
    tree->SetBranchAddress(#name, &name); \
  }
#define BOOK_BRANCH(t, name, outtree) \
  t name; (outtree).Branch(#name, &name);
#define BOOK_BRANCH_CONDITIONAL(t, name, outtree, cond) \
  t name; if (cond){ (outtree).Branch(#name, &name); }


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

void get2DParallelAndPerpendicularComponents(TVector3 axis, TVector3 ref, float& parallel, float& perp){
  TVector3 unitAxis = TVector3(axis.X(), axis.Y(), 0).Unit();
  TVector3 refPerp = TVector3(ref.X(), ref.Y(), 0);
  parallel = unitAxis.Dot(refPerp);
  perp = unitAxis.Cross(refPerp).Z();
}

void createGammaTrees(bool doData=true){
  gStyle->SetOptStat(0);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  SampleHelpers::theDataYear=2017;
  SampleHelpers::theDataPeriod="2017";
  SampleHelpers::theDataVersion=SampleHelpers::kCMSSW_9_4_X;
  float const btagWP = BtagHelpers::getBtagWP(BtagHelpers::kDeepCSV_Medium);
  float const totalLumi=150.f*1000.f;

  std::vector<TString> sampleList={
    "Run2017B-31Mar2018-v1_SinglePhoton",
    "Run2017C-31Mar2018-v1_SinglePhoton",
    "Run2017D-31Mar2018-v1_SinglePhoton",
    "Run2017E-31Mar2018-v1_SinglePhoton",
    "Run2017F-31Mar2018-v1_SinglePhoton",
    "Run2017F-09May2018-v1_SinglePhoton",

    "GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",
    "TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8",
    //"TTGamma_Dilept_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromT_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromTbar_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8",
    "WZG_TuneCP5_13TeV-amcatnlo-pythia8"
  };

  gSystem->Exec(Form("mkdir -p output/GammaTrees/%i", SampleHelpers::theDataYear));
  TString const stroutputcore = Form("output/GammaTrees/%i/[OUTFILECORE].root", SampleHelpers::theDataYear);
  TString const strinputcore = "/hadoop/cms/store/user/usarica/STOP_1L/Samples_2017/[DATE]/[OUTFILECORE]";
  for (auto const& strSample:sampleList){
    bool const isData = (strSample.BeginsWith("Run"));
    if ((!doData && isData) || (doData && !isData)) continue;
    TString strinput = strinputcore;
    HelperFunctions::replaceString<TString, const TString>(strinput, TString("[OUTFILECORE]"), strSample);
    switch (SampleHelpers::theDataYear){
    case 2017:
    {
      if (isData) HelperFunctions::replaceString<TString, const TString>(strinput, TString("[DATE]"), TString("190220"));
      else HelperFunctions::replaceString<TString, const TString>(strinput, TString("[DATE]"), TString("190224"));
      break;
    }
    default:
      assert(0);
    }

    TString stroutput = stroutputcore;
    HelperFunctions::replaceString<TString, const TString>(stroutput, TString("[OUTFILECORE]"), strSample);

    TFile* foutput = TFile::Open(stroutput, "recreate");
    TTree tout("AnalysisTree", "");

    TChain* tree = new TChain("SelectedTree");
    tree->Add(Form("%s/*", strinput.Data()));
    tree->SetBranchStatus("*", 0);

    // XSEC
    std::vector<float> genweights(1, 1.f);
    if (!isData){
      genweights = std::vector<float>(WeightVariables::nWeightTypes, 0.f);
      for (int iwt=WeightVariables::wCentral; iwt<WeightVariables::nWeightTypes; iwt++){
        TString bname = WeightVariables::getWeightName((WeightVariables::WeightType) iwt);
        tree->SetBranchStatus(bname, 1); tree->SetBranchAddress(bname, &(genweights.at(iwt)));
      }
    }

    // Extra weights
    DATA_SIMPLE_CONDITIONAL(float, xsec, !isData);
    DATA_SIMPLE_CONDITIONAL(float, weight_photons, !isData);
    DATA_SIMPLE_CONDITIONAL(float, weight_PU, !isData);
    DATA_SIMPLE_CONDITIONAL(float, weight_PU_SFUp, !isData);
    DATA_SIMPLE_CONDITIONAL(float, weight_PU_SFDn, !isData);
    // Event filters
    DATA_SIMPLE(bool, passHLTPath_SingleEle);
    DATA_SIMPLE(bool, passHLTPath_SingleMu);
    DATA_SIMPLE(bool, passHLTPath_SinglePhoton);
    DATA_SIMPLE(bool, passHLTPath_MET);
    DATA_SIMPLE(bool, passEventFilters);
    DATA_SIMPLE(bool, passGoodPrimaryVertex);
    // MET
    DATA_SIMPLE(float, pfmet);
    DATA_SIMPLE(float, pfmetPhi);
    DATA_SIMPLE_CONDITIONAL(float, pfmet_JECup, !isData);
    DATA_SIMPLE_CONDITIONAL(float, pfmetPhi_JECup, !isData);
    DATA_SIMPLE_CONDITIONAL(float, pfmet_JECdn, !isData);
    DATA_SIMPLE_CONDITIONAL(float, pfmetPhi_JECdn, !isData);
    // Muons
    DATA_CLASS(std::vector<int>, muons_id);
    DATA_CLASS(std::vector<long long>, muons_selectionBits);
    // Electrons
    DATA_CLASS(std::vector<int>, electrons_id);
    DATA_CLASS(std::vector<long long>, electrons_selectionBits);
    // Photons
    DATA_CLASS(std::vector<float>, photons_pt);
    DATA_CLASS(std::vector<float>, photons_eta);
    DATA_CLASS(std::vector<float>, photons_phi);
    DATA_CLASS(std::vector<float>, photons_mass);
    DATA_CLASS(std::vector<long long>, photons_selectionBits);
    // AK4Jets
    DATA_CLASS(std::vector<float>, ak4jets_pt);
    DATA_CLASS(std::vector<float>, ak4jets_eta);
    DATA_CLASS(std::vector<float>, ak4jets_phi);
    DATA_CLASS(std::vector<float>, ak4jets_mass);
    DATA_CLASS(std::vector<long long>, ak4jets_selectionBits);
    DATA_CLASS_CONDITIONAL(std::vector<float>, ak4jets_JECup, !isData);
    DATA_CLASS_CONDITIONAL(std::vector<float>, ak4jets_JECdn, !isData);
    DATA_CLASS_CONDITIONAL(std::vector<float>, ak4jets_JERup, !isData);
    DATA_CLASS_CONDITIONAL(std::vector<float>, ak4jets_JERdn, !isData);
    // AK8Jets
    DATA_CLASS(std::vector<float>, ak8jets_pt);
    DATA_CLASS(std::vector<float>, ak8jets_eta);
    DATA_CLASS(std::vector<float>, ak8jets_phi);
    DATA_CLASS(std::vector<float>, ak8jets_mass);
    DATA_CLASS(std::vector<long long>, ak8jets_selectionBits);

    // Gen. particles
    //DATA_CLASS_CONDITIONAL(std::vector<bool>, genparticles_isPromptFinalState, !isData);
    //DATA_CLASS_CONDITIONAL(std::vector<bool>, genparticles_isPromptDecayed, !isData);
    //DATA_CLASS_CONDITIONAL(std::vector<bool>, genparticles_isDirectPromptTauDecayProductFinalState, !isData);
    DATA_CLASS_CONDITIONAL(std::vector<bool>, genparticles_isHardProcess, !isData);
    //DATA_CLASS_CONDITIONAL(std::vector<bool>, genparticles_fromHardProcessFinalState, !isData);
    //DATA_CLASS_CONDITIONAL(std::vector<bool>, genparticles_fromHardProcessDecayed, !isData);
    //DATA_CLASS_CONDITIONAL(std::vector<bool>, genparticles_isDirectHardProcessTauDecayProductFinalState, !isData);
    //DATA_CLASS_CONDITIONAL(std::vector<bool>, genparticles_fromHardProcessBeforeFSR, !isData);
    //DATA_CLASS_CONDITIONAL(std::vector<bool>, genparticles_isLastCopy, !isData);
    //DATA_CLASS_CONDITIONAL(std::vector<bool>, genparticles_isLastCopyBeforeFSR, !isData);
    DATA_CLASS_CONDITIONAL(std::vector<int>, genparticles_id, !isData);
    DATA_CLASS_CONDITIONAL(std::vector<int>, genparticles_status, !isData);
    DATA_CLASS_CONDITIONAL(std::vector<float>, genparticles_px, !isData);
    DATA_CLASS_CONDITIONAL(std::vector<float>, genparticles_py, !isData);
    DATA_CLASS_CONDITIONAL(std::vector<float>, genparticles_pz, !isData);
    DATA_CLASS_CONDITIONAL(std::vector<float>, genparticles_mass, !isData);

    // Book outout tree branches
    if (!isData){
      for (int iwt=WeightVariables::wCentral; iwt<WeightVariables::nWeightTypes; iwt++){
        TString bname = WeightVariables::getWeightName((WeightVariables::WeightType) iwt);
        tout.Branch(bname, &(genweights.at(iwt)));
      }
    }
    else tout.Branch(WeightVariables::getWeightName(WeightVariables::wCentral), &(genweights.at(WeightVariables::wCentral)));
    tout.Branch("pfmet", &pfmet);
    tout.Branch("pfmetPhi", &pfmetPhi);
    if (!isData){
      tout.Branch("pfmet_JECup", &pfmet_JECup);
      tout.Branch("pfmetPhi_JECup", &pfmetPhi_JECup);
      tout.Branch("pfmet_JECdn", &pfmet_JECdn);
      tout.Branch("pfmetPhi_JECdn", &pfmetPhi_JECdn);

      tout.Branch("xsec", &xsec);
      tout.Branch("weight_photons", &weight_photons);
      tout.Branch("weight_PU", &weight_PU);
      tout.Branch("weight_PU_SFUp", &weight_PU_SFUp);
      tout.Branch("weight_PU_SFDn", &weight_PU_SFDn);
    }

    BOOK_BRANCH(unsigned int, nak4jets_preselected, tout);
    BOOK_BRANCH(float, genJetHT, tout);
    BOOK_BRANCH(float, jetHT, tout);
    BOOK_BRANCH(float, ak4jet_leadingpt_pt, tout);
    BOOK_BRANCH(float, photon_pt, tout);
    BOOK_BRANCH(float, photon_eta, tout);
    BOOK_BRANCH(float, photon_phi, tout);
    BOOK_BRANCH(float, photon_mass, tout);
    BOOK_BRANCH(float, uParallel, tout); // Component of the sum of vector transverse momenta of jets in the direction of the photon
    BOOK_BRANCH(float, uPerp, tout); // Component of the sum of vector transverse momenta of jets perp. to the direction of the photon
    BOOK_BRANCH(float, MET_Parallel, tout); // Component of MET in the direction of the photon
    BOOK_BRANCH(float, MET_Perp, tout); // Component of MET perp. to the direction of the photon
    BOOK_BRANCH_CONDITIONAL(float, uParallel_JECup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, uPerp_JECup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, uParallel_JECdn, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, uPerp_JECdn, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, uParallel_JERup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, uPerp_JERup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, uParallel_JERdn, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, uPerp_JERdn, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, MET_Parallel_JECup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, MET_Perp_JECup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, MET_Parallel_JECdn, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, MET_Perp_JECdn, tout, !isData);

    // Loop over the tree
    std::vector<float> sumCtrWgt(genweights.size(), 0.f);
    std::vector<float> sumCtrWgt_PU(genweights.size(), 0.f);
    std::vector<float> sumCtrWgt_PU_SFUp(genweights.size(), 0.f);
    std::vector<float> sumCtrWgt_PU_SFDn(genweights.size(), 0.f);
    std::vector<unsigned int> sumCtrOne(genweights.size(), 0.f);
    int const nEntries = tree->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      tree->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      if (ev==0 && !isData) tree->SetBranchStatus("xsec", 0);

      for (size_t iwt=0; iwt<genweights.size(); iwt++){
        sumCtrWgt[iwt] += genweights[iwt];
        sumCtrWgt_PU[iwt] += genweights[iwt]*weight_PU;
        sumCtrWgt_PU_SFUp[iwt] += genweights[iwt]*weight_PU_SFUp;
        sumCtrWgt_PU_SFDn[iwt] += genweights[iwt]*weight_PU_SFDn;
        sumCtrOne[iwt]++;
      }

      genJetHT = -1;

      nak4jets_preselected = 0;
      jetHT = -1;
      ak4jet_leadingpt_pt = -1;
      photon_pt = photon_eta = photon_phi = photon_mass = 0;

      // Gen particles
      size_t ngenparticles = 0;
      std::vector<MELAParticle> genparticles;
      //std::vector<MELAParticle*> genZlhe, genWlhe, genToplhe, genHardParticles;
      if (!isData){
        ngenparticles = genparticles_px->size(); genparticles.reserve(ngenparticles);
        for (size_t ip=0; ip<ngenparticles; ip++){
          const int& id = genparticles_id->at(ip);
          double E = pow(genparticles_px->at(ip), 2)+pow(genparticles_py->at(ip), 2)+pow(genparticles_pz->at(ip), 2);
          E += (genparticles_mass->at(ip)>=0.f ? 1.f : -1.f) * pow(genparticles_mass->at(ip), 2);
          E = sqrt(E);
          TLorentzVector v(genparticles_px->at(ip), genparticles_py->at(ip), genparticles_pz->at(ip), E);
          genparticles.emplace_back(id, v);
          genparticles.back().setGenStatus(genparticles_status->at(ip));
          if (genparticles_isHardProcess->at(ip)){
            if (PDGHelpers::isAQuark(id) || PDGHelpers::isAGluon(id)){
              if (genJetHT<0.) genJetHT=0;
              genJetHT += genparticles.back().pt();
            }
            //genHardParticles.push_back(&(genparticles.back()));
            //if (std::abs(id)==6) genToplhe.push_back(&(genparticles.back()));
            //else if (std::abs(id)==23) genZlhe.push_back(&(genparticles.back()));
            //else if (std::abs(id)==24) genWlhe.push_back(&(genparticles.back()));
          }
        }
      }

      bool doSkipEvent = !(
        passEventFilters
        && passGoodPrimaryVertex
        && (
        passHLTPath_SingleEle
        || passHLTPath_SingleMu
        || passHLTPath_SinglePhoton
        || passHLTPath_MET
          )
        );
      if (doSkipEvent) continue;

      // Reco particles
      size_t const nmuons = muons_selectionBits->size();
      for (size_t ip=0; ip<nmuons; ip++){
        if (HelperFunctions::test_bit(muons_selectionBits->at(ip), MuonSelectionHelpers::kVetoIDReco)){
          doSkipEvent = true;
          break;
        }
      }
      if (doSkipEvent) continue;

      size_t const nelectrons = electrons_selectionBits->size();
      for (size_t ip=0; ip<nelectrons; ip++){
        if (HelperFunctions::test_bit(electrons_selectionBits->at(ip), ElectronSelectionHelpers::kVetoIDReco)){
          doSkipEvent = true;
          break;
        }
      }
      if (doSkipEvent) continue;

      size_t const nphotons = photons_selectionBits->size();
      std::vector<MELAParticle> photons; photons.reserve(nphotons);
      bool hasSelectedPhotons = false;
      for (size_t ip=0; ip<nphotons; ip++){
        if (!HelperFunctions::test_bit(photons_selectionBits->at(ip), PhotonSelectionHelpers::kPreselection)) continue;
        TLorentzVector v; v.SetPtEtaPhiM(photons_pt->at(ip), photons_eta->at(ip), photons_phi->at(ip), photons_mass->at(ip));
        photons.emplace_back(0, v);
        if (!hasSelectedPhotons){
          hasSelectedPhotons = true;
          photon_pt = v.Pt();
          photon_eta = v.Eta();
          photon_phi = v.Phi();
          photon_mass = v.M();
        }
      }
      size_t nphotons_preselected = photons.size();
      if (nphotons_preselected!=1) doSkipEvent = true;
      if (doSkipEvent) continue;

      // MET vectors
      TLorentzVector vec_pfmet, vec_pfmet_JECup, vec_pfmet_JECdn;
      vec_pfmet.SetPtEtaPhiM(pfmet, 0, pfmetPhi, 0);
      if (!isData){
        vec_pfmet_JECup.SetPtEtaPhiM(pfmetPhi_JECup, 0, pfmetPhi_JECup, 0);
        vec_pfmet_JECdn.SetPtEtaPhiM(pfmetPhi_JECdn, 0, pfmetPhi_JECdn, 0);
      }

      // ak4 jet vectors
      TLorentzVector ak4jets_sumP(0, 0, 0, 0);
      TLorentzVector ak4jets_sumP_JECup(0, 0, 0, 0);
      TLorentzVector ak4jets_sumP_JECdn(0, 0, 0, 0);
      TLorentzVector ak4jets_sumP_JERup(0, 0, 0, 0);
      TLorentzVector ak4jets_sumP_JERdn(0, 0, 0, 0);
      size_t const nak4jets = ak4jets_selectionBits->size();
      std::vector<MELAParticle> ak4jets; ak4jets.reserve(nak4jets);
      bool hasSelectedJets = false;
      for (size_t ip=0; ip<nak4jets; ip++){
        if (!HelperFunctions::test_bit(ak4jets_selectionBits->at(ip), AK4JetSelectionHelpers::kPreselection)) continue;
        TLorentzVector v; v.SetPtEtaPhiM(ak4jets_pt->at(ip), ak4jets_eta->at(ip), ak4jets_phi->at(ip), ak4jets_mass->at(ip));
        ak4jets.emplace_back(0, v);
        if (!hasSelectedJets){
          hasSelectedJets = true;
          ak4jet_leadingpt_pt = v.Pt();
          jetHT = ak4jet_leadingpt_pt;
        }
        else{
          jetHT += v.Pt();
        }
        ak4jets_sumP += v;
        if (!isData){
          ak4jets_sumP_JECup += v*ak4jets_JECup->at(ip);
          ak4jets_sumP_JECdn += v*ak4jets_JECdn->at(ip);
          ak4jets_sumP_JERup += v*ak4jets_JERup->at(ip);
          ak4jets_sumP_JERdn += v*ak4jets_JERdn->at(ip);
        }
      }
      nak4jets_preselected = ak4jets.size();

      TVector3 photonAxis = photons.front().p4.Vect();
      get2DParallelAndPerpendicularComponents(photonAxis, ak4jets_sumP.Vect(), uParallel, uPerp);
      get2DParallelAndPerpendicularComponents(photonAxis, vec_pfmet.Vect(), MET_Parallel, MET_Perp);
      if (!isData){
        get2DParallelAndPerpendicularComponents(photonAxis, ak4jets_sumP_JECup.Vect(), uParallel_JECup, uPerp_JECup);
        get2DParallelAndPerpendicularComponents(photonAxis, ak4jets_sumP_JECdn.Vect(), uParallel_JECdn, uPerp_JECdn);
        get2DParallelAndPerpendicularComponents(photonAxis, ak4jets_sumP_JERup.Vect(), uParallel_JERup, uPerp_JERup);
        get2DParallelAndPerpendicularComponents(photonAxis, ak4jets_sumP_JERdn.Vect(), uParallel_JERdn, uPerp_JERdn);
        get2DParallelAndPerpendicularComponents(photonAxis, vec_pfmet_JECup.Vect(), MET_Parallel_JECup, MET_Perp_JECup);
        get2DParallelAndPerpendicularComponents(photonAxis, vec_pfmet_JECdn.Vect(), MET_Parallel_JECdn, MET_Perp_JECdn);
      }

      tout.Fill();
    }

    TTree tout_md("metadata", "");
    tout_md.Branch("nAllEvents", &(sumCtrOne.at(0)));
    MELAout << "Number of output events: " << tout.GetEntries() << " / " << sumCtrOne[0] << endl;
    if (!isData){
      for (size_t iwt=0; iwt<genweights.size(); iwt++){
        MELAout << "Sum of weights: " << sumCtrWgt[iwt] << " / " << sumCtrOne[iwt] << endl;
        MELAout << "Sum of weights with PU: " << sumCtrWgt_PU[iwt] << " | " << sumCtrWgt_PU_SFUp[iwt] << " | " << sumCtrWgt_PU_SFDn[iwt] << endl;
        for (int iwt=WeightVariables::wCentral; iwt<WeightVariables::nWeightTypes; iwt++){
          TString bname = WeightVariables::getWeightName((WeightVariables::WeightType) iwt);
          tout_md.Branch(bname, &(sumCtrWgt.at(iwt)));
          tout_md.Branch(bname+"_PU", &(sumCtrWgt_PU.at(iwt)));
          tout_md.Branch(bname+"_PU_SFUp", &(sumCtrWgt_PU_SFUp.at(iwt)));
          tout_md.Branch(bname+"_PU_SFDn", &(sumCtrWgt_PU_SFDn.at(iwt)));
        }
      }
    }
    tout_md.Fill();
    foutput->WriteTObject(&tout_md);

    delete tree;
    foutput->WriteTObject(&tout);
    foutput->Close();
  }

}
