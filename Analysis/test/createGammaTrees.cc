#include <cassert>
#include "common_includes.h"
#include "PDGHelpers.h"
#include "MELAParticle.h"
#include "RooMsgService.h"
#include "RooRelBW2ProngPdf.h"
#include "RooRelBW3ProngPdf.h"
#include "RooArgList.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooDataSet.h"
#include "RooVoigtian.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TNumericUtil.hh"
#include <HiggsAnalysis/CombinedLimit/interface/AsymPow.h>
#include <limits>
#include "TMatrixDSym.h"
#include "TChain.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLegend.h"


using namespace RooFit;


#define DATA_CLASS(t, name) \
  t* name=nullptr; \
  tree->SetBranchStatus(#name, 1); \
  tree->SetBranchAddress(#name, &name);
#define DATA_SIMPLE(t, name) \
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
  ExtendedBinning binning;
  float val;
  std::vector<float> weights;

  Variable() : binning(), val(0), weights(WeightVariables::nWeightTypes, 0.f){}
  Variable(Variable const& other) : name(other.name), title(other.title), binning(other.binning), val(other.val), weights(WeightVariables::nWeightTypes, 0.f){}
  Variable(TString name_, TString title_) : name(name_), title(title_), binning(title), val(0), weights(WeightVariables::nWeightTypes, 0.f){}
  Variable(TString name_, TString title_, unsigned int nbins_, float min_, float max_) : name(name_), title(title_), binning(nbins_, min_, max_, title), val(0), weights(WeightVariables::nWeightTypes, 0.f){}
  Variable(TString name_, TString title_, ExtendedBinning const& binning_) : name(name_), title(title_), binning(binning_), val(0), weights(WeightVariables::nWeightTypes, 0.f){}

  void setVal(float v, std::vector<float> const& wgts){
    weights=wgts;
    val=v;
    if (binning.isValid()){
      unsigned int const nbins = binning.getNbins();
      float const min = binning.getMin();
      float const max = binning.getMax();
      if (val>=max) val=max-binning.getBinWidth(nbins-1)/2.;
      if (val<=min) val=min+binning.getBinWidth(0)/2.;
    }
  }
  void setVal(float v, float const& wgt){
    std::vector<float> wgts(1, wgt);
    this->setVal(v, wgts);
  }

  void reset(){
    if (binning.isValid()){
      float const min = binning.getMin();
      val=min - fabs(min);
    }
    else val=0;
    weights=std::vector<float>(weights.size(), 0.f);
  }

  // Proxy functions
  double* getBinning(){ return binning.getBinning(); }
  const double* getBinning() const{ return binning.getBinning(); }
  unsigned int getNbins() const{ return binning.getNbins(); }

};

void get2DParallelAndPerpendicularComponents(TVector3 axis, TVector3 ref, float& parallel, float& perp){
  TVector3 unitAxis = TVector3(axis.X(), axis.Y(), 0).Unit();
  TVector3 refPerp = TVector3(ref.X(), ref.Y(), 0);
  parallel = unitAxis.Dot(refPerp);
  perp = unitAxis.Cross(refPerp).Z();
}

void createGammaTrees(TString strSampleSet){
  gStyle->SetOptStat(0);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  SampleHelpers::theDataYear=2017;
  SampleHelpers::theDataPeriod="2017";
  SampleHelpers::theDataVersion=SampleHelpers::kCMSSW_9_4_X;
  //float const btagWP = BtagHelpers::getBtagWP(BtagHelpers::kDeepCSV_Medium);
  //float const totalLumi=150.f*1000.f;

  // For data only
  unordered_map<TString, DuplicateEventHandler> dupEventHandles;

  std::vector<TString> sampleList={
    "Run2017B-31Mar2018-v1_SinglePhoton",
    "Run2017C-31Mar2018-v1_SinglePhoton",
    "Run2017D-31Mar2018-v1_SinglePhoton",
    "Run2017E-31Mar2018-v1_SinglePhoton",
    "Run2017F-31Mar2018-v1_SinglePhoton",
    "Run2017F-09May2018-v1_SinglePhoton",
    "Run2017B-31Mar2018-v1_SingleElectron",
    "Run2017C-31Mar2018-v1_SingleElectron",
    "Run2017D-31Mar2018-v1_SingleElectron",
    "Run2017E-31Mar2018-v1_SingleElectron",
    "Run2017F-31Mar2018-v1_SingleElectron",
    "Run2017F-09May2018-v1_SingleElectron",
    "Run2017B-31Mar2018-v1_SingleMuon",
    "Run2017C-31Mar2018-v1_SingleMuon",
    "Run2017D-31Mar2018-v1_SingleMuon",
    "Run2017E-31Mar2018-v1_SingleMuon",
    "Run2017F-31Mar2018-v1_SingleMuon",
    "Run2017F-09May2018-v1_SingleMuon",
    "Run2017B-31Mar2018-v1_MET",
    "Run2017C-31Mar2018-v1_MET",
    "Run2017D-31Mar2018-v1_MET",
    "Run2017E-31Mar2018-v1_MET",
    "Run2017F-31Mar2018-v1_MET",
    "Run2017F-09May2018-v1_MET",

    "GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",

    "QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8",

    "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    "TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8",
    //"TTGamma_Dilept_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromT_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromTbar_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8",

    "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WZG_TuneCP5_13TeV-amcatnlo-pythia8",

    "ZJetsToNuNu_HT-100To200_13TeV-madgraph",
    "ZJetsToNuNu_HT-200To400_13TeV-madgraph",
    "ZJetsToNuNu_HT-400To600_13TeV-madgraph",
    "ZJetsToNuNu_HT-600To800_13TeV-madgraph",
    "ZJetsToNuNu_HT-800To1200_13TeV-madgraph",
    "ZJetsToNuNu_HT-1200To2500_13TeV-madgraph",
    "ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph"
  };
  {
    std::vector<TString> tmplist;
    for (TString const& strSample:sampleList){
      bool const isData = (strSample.BeginsWith("Run"));
      if (
        (strSampleSet == "Data" && isData)
        ||
        (strSampleSet == "MC" && !isData)
        ||
        (strSample.BeginsWith(strSampleSet))
        ) tmplist.push_back(strSample);
    }
    std::swap(sampleList, tmplist);
  }
  for (TString const& strSample:sampleList){
    if (strSample.Contains("Run")){
      TString key, junk;
      HelperFunctions::splitOption(strSample, key, junk, '_');
      dupEventHandles[key]=DuplicateEventHandler();
      MELAout << "Created a new duplicate event handle for key " << key << endl;
    }
  }

  gSystem->Exec(Form("mkdir -p output/GammaTrees/%i", SampleHelpers::theDataYear));
  TString const stroutputcore = Form("output/GammaTrees/%i", SampleHelpers::theDataYear);
  TString const strinputcore = "/nfs-7/userdata/usarica/STOP-1L/Samples_2017/[DATE]";
  for (auto const& strSample:sampleList){
    bool const isData = (strSample.BeginsWith("Run"));
    TString strinput = strinputcore;
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

    DuplicateEventHandler* dupEventHandle=nullptr;
    if (isData){
      TString key, junk;
      HelperFunctions::splitOption(strSample, key, junk, '_');
      dupEventHandle = &(dupEventHandles[key]);
    }

    TString stroutput = stroutputcore + "/" + strSample + ".root";

    TFile* foutput = TFile::Open(stroutput, "recreate");
    TTree tout("AnalysisTree", "");

    TChain* tree = new TChain("SelectedTree");
    std::vector<TString> inputList = SampleHelpers::lsdir(strinput);
    for (auto const& sdir:inputList){
      if (sdir.Contains(strSample)){
        TString strtmp = Form("%s/%s/*", strinput.Data(), sdir.Data());
        tree->Add(strtmp);
        MELAout << "Adding " << strtmp << endl;
      }
    }

    MELAout << "Disabling all branches before enabling them back again..." << endl;
    tree->SetBranchStatus("*", 0);

    MELAout << "Setting branch addresses.." << endl;
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
    DATA_SIMPLE_CONDITIONAL(RunNumber_t, RunNumber, isData);
    DATA_SIMPLE_CONDITIONAL(Lumisection_t, LumiSection, isData);
    DATA_SIMPLE_CONDITIONAL(EventNumber_t, EventNumber, isData);
    DATA_SIMPLE_CONDITIONAL(float, xsec, !isData);
    DATA_SIMPLE_CONDITIONAL(float, genMET, !isData);
    DATA_SIMPLE_CONDITIONAL(float, genMETPhi, !isData);
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

    MELAout << "Booking output tree branches..." << endl;
    // Book outout tree branches
    if (!isData){
      for (int iwt=WeightVariables::wCentral; iwt<WeightVariables::nWeightTypes; iwt++){
        TString bname = WeightVariables::getWeightName((WeightVariables::WeightType) iwt);
        tout.Branch(bname, &(genweights.at(iwt)));
      }
    }
    //else tout.Branch(WeightVariables::getWeightName(WeightVariables::wCentral), &(genweights.at(WeightVariables::wCentral)));
    tout.Branch("pfmet", &pfmet);
    tout.Branch("pfmetPhi", &pfmetPhi);
    tout.Branch("passHLTPath_SinglePhoton", &passHLTPath_SinglePhoton);
    tout.Branch("passHLTPath_SingleEle", &passHLTPath_SingleEle);
    tout.Branch("passHLTPath_MET", &passHLTPath_MET);
    tout.Branch("passHLTPath_SingleMu", &passHLTPath_SingleMu);
    if (!isData){
      tout.Branch("pfmet_JECup", &pfmet_JECup);
      tout.Branch("pfmetPhi_JECup", &pfmetPhi_JECup);
      tout.Branch("pfmet_JECdn", &pfmet_JECdn);
      tout.Branch("pfmetPhi_JECdn", &pfmetPhi_JECdn);

      tout.Branch("xsec", &xsec);
      tout.Branch("genMET", &genMET);
      tout.Branch("genMETPhi", &genMETPhi);
      tout.Branch("weight_photons", &weight_photons);
      tout.Branch("weight_PU", &weight_PU);
      tout.Branch("weight_PU_SFUp", &weight_PU_SFUp);
      tout.Branch("weight_PU_SFDn", &weight_PU_SFDn);
    }

    BOOK_BRANCH(unsigned int, nak4jets_preselected, tout);
    BOOK_BRANCH(float, jetHT, tout);
    BOOK_BRANCH(float, ak4jets_leadingpt_pt, tout);
    BOOK_BRANCH(float, photon_pt, tout);
    BOOK_BRANCH(float, photon_eta, tout);
    BOOK_BRANCH(float, photon_phi, tout);
    BOOK_BRANCH(float, photon_mass, tout);
    BOOK_BRANCH(float, uParallel, tout); // Component of the sum of vector transverse momenta of jets in the direction of the photon
    BOOK_BRANCH(float, uPerp, tout); // Component of the sum of vector transverse momenta of jets perp. to the direction of the photon
    BOOK_BRANCH(float, MET_Parallel, tout); // Component of MET in the direction of the photon
    BOOK_BRANCH(float, MET_Perp, tout); // Component of MET perp. to the direction of the photon
    BOOK_BRANCH_CONDITIONAL(float, genJetHT, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(unsigned int, nak4jets_preselected_JECup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, jetHT_JECup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, ak4jets_leadingpt_pt_JECup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(unsigned int, nak4jets_preselected_JECdn, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, jetHT_JECdn, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, ak4jets_leadingpt_pt_JECdn, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(unsigned int, nak4jets_preselected_JERup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, jetHT_JERup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, ak4jets_leadingpt_pt_JERup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(unsigned int, nak4jets_preselected_JERdn, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, jetHT_JERdn, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, ak4jets_leadingpt_pt_JERdn, tout, !isData);
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
    BOOK_BRANCH_CONDITIONAL(float, MET_Parallel_JERup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, MET_Perp_JERup, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, MET_Parallel_JERdn, tout, !isData);
    BOOK_BRANCH_CONDITIONAL(float, MET_Perp_JERdn, tout, !isData);

    // Loop over the tree
    std::vector<float> sumCtrWgt(genweights.size(), 0.f);
    std::vector<float> sumCtrWgt_PU(genweights.size(), 0.f);
    std::vector<float> sumCtrWgt_PU_SFUp(genweights.size(), 0.f);
    std::vector<float> sumCtrWgt_PU_SFDn(genweights.size(), 0.f);
    std::vector<unsigned int> sumCtrOne(genweights.size(), 0.f);
    int const nEntries = tree->GetEntries();
    MELAout << "Starting to loop over " << nEntries << " events" << endl;
    for (int ev=0; ev<nEntries; ev++){
      tree->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      if (ev==0 && !isData) tree->SetBranchStatus("xsec", 0);
      if (isData && !dupEventHandle->isUnique(RunNumber, LumiSection, EventNumber)) continue;

      for (size_t iwt=0; iwt<genweights.size(); iwt++){
        sumCtrWgt[iwt] += genweights[iwt];
        sumCtrWgt_PU[iwt] += genweights[iwt]*weight_PU;
        sumCtrWgt_PU_SFUp[iwt] += genweights[iwt]*weight_PU_SFUp;
        sumCtrWgt_PU_SFDn[iwt] += genweights[iwt]*weight_PU_SFDn;
        sumCtrOne[iwt]++;
      }


      genJetHT = -1;
      photon_pt = photon_eta = photon_phi = photon_mass = 0;
      nak4jets_preselected = nak4jets_preselected_JECup = nak4jets_preselected_JECdn = nak4jets_preselected_JERup = nak4jets_preselected_JERdn = 0;
      jetHT = jetHT_JECup = jetHT_JECdn = jetHT_JERup = jetHT_JERdn = -1;
      ak4jets_leadingpt_pt = ak4jets_leadingpt_pt_JECup = ak4jets_leadingpt_pt_JECdn = ak4jets_leadingpt_pt_JERup = ak4jets_leadingpt_pt_JERdn = -1;


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
        vec_pfmet_JECup.SetPtEtaPhiM(pfmet_JECup, 0, pfmetPhi_JECup, 0);
        vec_pfmet_JECdn.SetPtEtaPhiM(pfmet_JECdn, 0, pfmetPhi_JECdn, 0);
      }

      // ak4 jet vectors
      TLorentzVector ak4jets_sumP(0, 0, 0, 0);
      TLorentzVector ak4jets_sumP_JECup(0, 0, 0, 0);
      TLorentzVector ak4jets_sumP_JECdn(0, 0, 0, 0);
      TLorentzVector ak4jets_sumP_JERup(0, 0, 0, 0);
      TLorentzVector ak4jets_sumP_JERdn(0, 0, 0, 0);
      size_t const nak4jets = ak4jets_selectionBits->size();
      // Nominal, JECup, JECdn, JERup, JERdn
      for (unsigned int js=0; js<5; js++){
        if (isData && js!=0) break;

        AK4JetSelectionHelpers::SelectionBits selFlag;
        TLorentzVector* ak4jets_sumP_ptr = nullptr;
        unsigned int* nak4jets_preselected_ptr = nullptr;
        float* ak4jets_leadingpt_pt_ptr = nullptr;
        float* jetHT_ptr = nullptr;
        switch (js){
        case 1:
          selFlag = AK4JetSelectionHelpers::kPreselection_JECUp;
          ak4jets_sumP_ptr = &ak4jets_sumP_JECup;
          nak4jets_preselected_ptr = &nak4jets_preselected_JECup;
          ak4jets_leadingpt_pt_ptr = &ak4jets_leadingpt_pt_JECup;
          jetHT_ptr = &jetHT_JECup;
          break;
        case 2:
          selFlag = AK4JetSelectionHelpers::kPreselection_JECDn;
          ak4jets_sumP_ptr = &ak4jets_sumP_JECdn;
          nak4jets_preselected_ptr = &nak4jets_preselected_JECdn;
          ak4jets_leadingpt_pt_ptr = &ak4jets_leadingpt_pt_JECdn;
          jetHT_ptr = &jetHT_JECdn;
          break;
        case 3:
          selFlag = AK4JetSelectionHelpers::kPreselection_JERUp;
          ak4jets_sumP_ptr = &ak4jets_sumP_JERup;
          nak4jets_preselected_ptr = &nak4jets_preselected_JERup;
          ak4jets_leadingpt_pt_ptr = &ak4jets_leadingpt_pt_JERup;
          jetHT_ptr = &jetHT_JERup;
          break;
        case 4:
          selFlag = AK4JetSelectionHelpers::kPreselection_JERDn;
          ak4jets_sumP_ptr = &ak4jets_sumP_JERdn;
          nak4jets_preselected_ptr = &nak4jets_preselected_JERdn;
          ak4jets_leadingpt_pt_ptr = &ak4jets_leadingpt_pt_JERdn;
          jetHT_ptr = &jetHT_JERdn;
          break;
        default:
          selFlag = AK4JetSelectionHelpers::kPreselection;
          ak4jets_sumP_ptr = &ak4jets_sumP;
          nak4jets_preselected_ptr = &nak4jets_preselected;
          ak4jets_leadingpt_pt_ptr = &ak4jets_leadingpt_pt;
          jetHT_ptr = &jetHT;
          break;
        }
        bool hasSelectedJets = false;
        for (size_t ip=0; ip<nak4jets; ip++){
          if (!HelperFunctions::test_bit(ak4jets_selectionBits->at(ip), selFlag)) continue;
          float mult = 1;
          switch (js){
          case 1:
            mult = ak4jets_JECup->at(ip);
            break;
          case 2:
            mult = ak4jets_JECdn->at(ip);
            break;
          case 3:
            mult = ak4jets_JERup->at(ip);
            break;
          case 4:
            mult = ak4jets_JERdn->at(ip);
            break;
          }
          TLorentzVector v; v.SetPtEtaPhiM(ak4jets_pt->at(ip)*mult, ak4jets_eta->at(ip), ak4jets_phi->at(ip), ak4jets_mass->at(ip)*mult);
          if (!hasSelectedJets){
            hasSelectedJets = true;
            *jetHT_ptr = *ak4jets_leadingpt_pt_ptr = v.Pt();
          }
          else{
            *jetHT_ptr += v.Pt();
          }
          *ak4jets_sumP_ptr += v;
          *nak4jets_preselected_ptr += 1;
        }
      }

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
        MET_Parallel_JERup = MET_Parallel_JERdn = MET_Parallel;
        MET_Perp_JERup = MET_Perp_JERdn = MET_Perp;
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


void getCorrections_DataMC(){
  std::vector<TString> sampleList_Data={
    "Run2017B-31Mar2018-v1",
    "Run2017C-31Mar2018-v1",
    "Run2017D-31Mar2018-v1",
    "Run2017E-31Mar2018-v1",
    "Run2017F-31Mar2018-v1",
    "Run2017F-09May2018-v1"
  };
  std::vector<TString> sampleList_MC={
    "GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",

    "QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8",

    //"TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    "TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8",
    //"TTGamma_Dilept_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromT_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromTbar_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8",

    "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WZG_TuneCP5_13TeV-amcatnlo-pythia8",

    "ZJetsToNuNu_HT-100To200_13TeV-madgraph",
    "ZJetsToNuNu_HT-200To400_13TeV-madgraph",
    "ZJetsToNuNu_HT-400To600_13TeV-madgraph",
    "ZJetsToNuNu_HT-600To800_13TeV-madgraph",
    "ZJetsToNuNu_HT-800To1200_13TeV-madgraph",
    "ZJetsToNuNu_HT-1200To2500_13TeV-madgraph",
    "ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph"
  };

  SampleHelpers::theDataYear=2017;
  SampleHelpers::theDataPeriod="2017";
  SampleHelpers::theDataVersion=SampleHelpers::kCMSSW_9_4_X;
  float const btagWP = BtagHelpers::getBtagWP(BtagHelpers::kDeepCSV_Medium);

  TString const strinputcore = Form("output/GammaTrees/%i/[OUTFILECORE].root", SampleHelpers::theDataYear);
  TString const stroutputcore = Form("output/GammaTrees/%i/plots/scale_[OUTFILECORE]_MC.root", SampleHelpers::theDataYear);
  gSystem->Exec(Form("mkdir -p output/GammaTrees/%i/plots", SampleHelpers::theDataYear));

  // Establish the variables
  std::vector<Variable> varlist={
    Variable("nak4jets_preselected", "N_{jets}", 10, 0, 10),
    Variable("jetHT", "Jet H_{T} (GeV)", ExtendedBinning({ -1, 0, 30, 40, 50, 60, 75, 90, 150, 300, 600, 1000, 3000 })),
    Variable("photon_pt", "p_{T}^{#gamma} (GeV)", ExtendedBinning({ 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 300, 400, 600, 1000 }))
  };
  auto& xvar = varlist.at(0);
  auto& yvar = varlist.at(1);
  auto& zvar = varlist.at(2);

  // Record MC histograms
  unsigned int nSysts=7;
  vector<vector<TH3F>> hMClist(sampleList_MC.size(), vector<TH3F>());
  for (size_t is=0; is<sampleList_MC.size(); is++){
    auto& v=hMClist.at(is);
    v.reserve(nSysts);
    for (unsigned int isyst=0; isyst<nSysts; isyst++){
      TString hname = xvar.name + "_" + yvar.name + "_" + zvar.name + "_" + sampleList_MC.at(is);
      if (isyst==0) hname += "_Nominal";
      else if (isyst==1) hname += "_JECup";
      else if (isyst==2) hname += "_JECdn";
      else if (isyst==3) hname += "_JERup";
      else if (isyst==4) hname += "_JERdn";
      else if (isyst==5) hname += "_PUup";
      else if (isyst==6) hname += "_PUdn";
      v.emplace_back(
        hname, "",
        xvar.getNbins(), xvar.getBinning(),
        yvar.getNbins(), yvar.getBinning(),
        zvar.getNbins(), zvar.getBinning()
      );
      v.back().Sumw2();
      v.back().GetXaxis()->SetTitle(xvar.title);
      v.back().GetYaxis()->SetTitle(yvar.title);
      v.back().GetZaxis()->SetTitle(zvar.title);
    }
  }
  for (size_t is=0; is<sampleList_MC.size(); is++){
    TString strinput = strinputcore;
    HelperFunctions::replaceString<TString, const TString>(strinput, "[OUTFILECORE]", sampleList_MC.at(is));
    MELAout << "Opening file " << strinput << endl;
    TFile* finput = TFile::Open(strinput, "read");
    TTree* tree = (TTree*) finput->Get("AnalysisTree"); tree->SetBranchStatus("*", 0);

    std::vector<float> genweights(WeightVariables::nWeightTypes, 0);
    for (int iwt=WeightVariables::wCentral; iwt<WeightVariables::nWeightTypes; iwt++){
      TString bname = WeightVariables::getWeightName((WeightVariables::WeightType) iwt);
      tree->SetBranchStatus(bname, 1);
      tree->SetBranchAddress(bname, &(genweights.at(iwt)));
    }
    DATA_SIMPLE(float, xsec);
    DATA_SIMPLE(float, weight_photons);
    DATA_SIMPLE(float, weight_PU);
    DATA_SIMPLE(float, weight_PU_SFUp);
    DATA_SIMPLE(float, weight_PU_SFDn);

    DATA_SIMPLE(unsigned int, nak4jets_preselected);
    DATA_SIMPLE(unsigned int, nak4jets_preselected_JECup);
    DATA_SIMPLE(unsigned int, nak4jets_preselected_JECdn);
    DATA_SIMPLE(unsigned int, nak4jets_preselected_JERup);
    DATA_SIMPLE(unsigned int, nak4jets_preselected_JERdn);
    DATA_SIMPLE(float, photon_pt);
    DATA_SIMPLE(float, jetHT);
    DATA_SIMPLE(float, jetHT_JECup);
    DATA_SIMPLE(float, jetHT_JECdn);
    DATA_SIMPLE(float, jetHT_JERup);
    DATA_SIMPLE(float, jetHT_JERdn);
    const int nEntries = tree->GetEntries();
    bool firstEvent=false;
    for (int ev=0; ev<nEntries; ev++){
      tree->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      for (unsigned int isyst=0; isyst<nSysts; isyst++){
        float totwgt = genweights.at(0)*weight_photons*xsec;
        if (isyst<5) totwgt *= weight_PU;
        else if (isyst==5) totwgt *= weight_PU_SFUp;
        else if (isyst==6) totwgt *= weight_PU_SFDn;
        else assert(0);

        unsigned int* nak4jets_preselected_ptr = &nak4jets_preselected;
        float* jetHT_ptr = &jetHT;
        switch (isyst){
        case 1:
          nak4jets_preselected_ptr = &nak4jets_preselected_JECup;
          jetHT_ptr = &jetHT_JECup;
          break;
        case 2:
          nak4jets_preselected_ptr = &nak4jets_preselected_JECdn;
          jetHT_ptr = &jetHT_JECdn;
          break;
        case 3:
          nak4jets_preselected_ptr = &nak4jets_preselected_JERup;
          jetHT_ptr = &jetHT_JERup;
          break;
        case 4:
          nak4jets_preselected_ptr = &nak4jets_preselected_JERdn;
          jetHT_ptr = &jetHT_JERdn;
          break;
        }

        if (firstEvent){
          MELAout << "First event, systematic " << isyst << ":\n";
          MELAout << "nak4jets_preselected: " << *nak4jets_preselected_ptr << " (address: " << nak4jets_preselected_ptr << ")" << endl;
          MELAout << "jetHT: " << *jetHT_ptr << " (address: " << jetHT_ptr << ")" << endl;
          MELAout << "photon: " << photon_pt << " (address: " << &photon_pt << ")" << endl;
        }

        xvar.setVal(*nak4jets_preselected_ptr, 1);
        yvar.setVal(*jetHT_ptr, 1);
        zvar.setVal(photon_pt, 1);

        if (firstEvent){
          MELAout << "Filling " << hMClist.at(is).at(isyst).GetName() << endl;
          MELAout << "\t- Fill variables: " << xvar.val << ", " << yvar.val << ", " << zvar.val << endl;
          MELAout << "\t- Fill weight: " << totwgt << endl;
        }

        hMClist.at(is).at(isyst).Fill(xvar.val, yvar.val, zvar.val, totwgt);

        xvar.reset();
        yvar.reset();
        zvar.reset();
      }
      if (firstEvent) firstEvent=false;
    }
    // Scale MC histograms
    {
      TTree* tree_md = (TTree*) finput->Get("metadata");
      float /*sumWgts=0, */sumWgts_PU=0, sumWgts_PU_SFUp=0, sumWgts_PU_SFDn=0;
      TString bname = WeightVariables::getWeightName(WeightVariables::wCentral);
      //tree_md->SetBranchAddress(bname, &sumWgts);
      tree_md->SetBranchAddress(bname+"_PU", &sumWgts_PU);
      tree_md->SetBranchAddress(bname+"_PU_SFUp", &sumWgts_PU_SFUp);
      tree_md->SetBranchAddress(bname+"_PU_SFDn", &sumWgts_PU_SFDn);
      tree_md->GetEntry(0);
      if (sumWgts_PU<=0.f) MELAout << "ERROR: sumWgts_PU = " << sumWgts_PU << endl;
      else{
        for (unsigned int isyst=0; isyst<nSysts; isyst++){
          auto& h = hMClist.at(is).at(isyst);
          float* sumWgts_ptr = &sumWgts_PU;
          if (isyst==5) sumWgts_ptr = &sumWgts_PU_SFUp;
          else if (isyst==6) sumWgts_ptr = &sumWgts_PU_SFDn;
          h.Scale(1./(*sumWgts_ptr));
          MELAout << "Scaled " << h.GetName() << " by 1 / " << (*sumWgts_ptr) << endl;
        }
      }
    }

    finput->Close();
  }

  for (auto const& sample_data:sampleList_Data){
    TString stroutput = stroutputcore;
    HelperFunctions::replaceString<TString, const TString>(stroutput, "[OUTFILECORE]", sample_data);
    MELAout << "Creating output file " << stroutput << endl;
    TFile* foutput = TFile::Open(stroutput, "recreate");

    // Record data histograms
    TH3F hdata = TH3F(
      xvar.name+"_"+yvar.name+"_"+zvar.name+"_Data", "",
      xvar.getNbins(), xvar.getBinning(),
      yvar.getNbins(), yvar.getBinning(),
      zvar.getNbins(), zvar.getBinning()
    );
    hdata.GetXaxis()->SetTitle(xvar.title);
    hdata.GetYaxis()->SetTitle(yvar.title);
    hdata.GetZaxis()->SetTitle(zvar.title);
    {
      TString strinput = strinputcore;
      TChain* tree = new TChain("AnalysisTree");
      HelperFunctions::replaceString<TString, const TString>(strinput, "/[OUTFILECORE].root", "");
      std::vector<TString> inputList = SampleHelpers::lsdir(strinput);
      //MELAout << "Input directories: " << inputList << endl;
      for (auto const& sfile:inputList){
        if (sfile.Contains(sample_data)){
          TString strtmp = Form("%s/%s", strinput.Data(), sfile.Data());
          tree->Add(strtmp);
          MELAout << "Adding " << strtmp << endl;
        }
      }
      DATA_SIMPLE(unsigned int, nak4jets_preselected);
      DATA_SIMPLE(float, photon_pt);
      DATA_SIMPLE(float, jetHT);
      int nEntries = tree->GetEntries();
      for (int ev=0; ev<nEntries; ev++){
        tree->GetEntry(ev);
        HelperFunctions::progressbar(ev, nEntries);

        xvar.setVal(nak4jets_preselected, 1);
        yvar.setVal(jetHT, 1);
        zvar.setVal(photon_pt, 1);

        hdata.Fill(xvar.val, yvar.val, zvar.val, 1);

        xvar.reset();
        yvar.reset();
        zvar.reset();
      }
      delete tree;
    }
    {
      foutput->WriteTObject(&hdata);
      {
        int iz=hdata.GetZaxis()->FindBin(130);
        float dataYield = hdata.Integral(1, hdata.GetNbinsX(), 1, hdata.GetNbinsY(), iz, hdata.GetNbinsZ());
        MELAout << "Data yield above pTgamma 130: " << dataYield << endl;
      }

      for (unsigned int isyst=0; isyst<nSysts; isyst++){
        TString strappend;
        if (isyst==0) strappend="_Nominal";
        else if (isyst==1) strappend="_JECup";
        else if (isyst==2) strappend="_JECdn";
        else if (isyst==3) strappend="_JERup";
        else if (isyst==4) strappend="_JERdn";
        else if (isyst==5) strappend="_PUup";
        else if (isyst==6) strappend="_PUdn";
        TH3F* htotMC_syst = (TH3F*) hdata.Clone(xvar.name+"_"+yvar.name+"_"+zvar.name+"_totMC"+strappend);
        htotMC_syst->Reset("ICESM");
        htotMC_syst->Sumw2();
        for (auto const& v:hMClist) htotMC_syst->Add(&(v.at(isyst)));
        float dataYield = hdata.Integral(1, hdata.GetNbinsX(), 1, hdata.GetNbinsY(), 1, hdata.GetNbinsZ());
        float MCyield = htotMC_syst->Integral(1, htotMC_syst->GetNbinsX(), 1, htotMC_syst->GetNbinsY(), 1, htotMC_syst->GetNbinsZ());
        float data_MC_scale = dataYield/MCyield;
        MELAout << "Data/MC = " << dataYield << " / " << MCyield << " = " << data_MC_scale << endl;
        //htotMC_syst->Scale(data_MC_scale);
        foutput->WriteTObject(htotMC_syst);
        {
          int iz=htotMC_syst->GetZaxis()->FindBin(130);
          float MCyieldn = htotMC_syst->Integral(1, htotMC_syst->GetNbinsX(), 1, htotMC_syst->GetNbinsY(), iz, htotMC_syst->GetNbinsZ());
          MELAout << "MC yield above pTgamma 130: " << MCyieldn << endl;
        }
        delete htotMC_syst;
      }
    }

    foutput->Close();
  }
}

void getResidualCorrections_DataMC(){
  std::vector<TString> sampleList_Data={
    "Run2017B-31Mar2018-v1",
    "Run2017C-31Mar2018-v1",
    "Run2017D-31Mar2018-v1",
    "Run2017E-31Mar2018-v1",
    "Run2017F-31Mar2018-v1",
    "Run2017F-09May2018-v1"
  };
  std::vector<TString> sampleList_MC={
    "GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",

    "QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8",

    //"TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    "TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8",
    //"TTGamma_Dilept_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromT_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromTbar_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8",

    "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WZG_TuneCP5_13TeV-amcatnlo-pythia8",

    "ZJetsToNuNu_HT-100To200_13TeV-madgraph",
    "ZJetsToNuNu_HT-200To400_13TeV-madgraph",
    "ZJetsToNuNu_HT-400To600_13TeV-madgraph",
    "ZJetsToNuNu_HT-600To800_13TeV-madgraph",
    "ZJetsToNuNu_HT-800To1200_13TeV-madgraph",
    "ZJetsToNuNu_HT-1200To2500_13TeV-madgraph",
    "ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph"
  };

  SampleHelpers::theDataYear=2017;
  SampleHelpers::theDataPeriod="2017";
  SampleHelpers::theDataVersion=SampleHelpers::kCMSSW_9_4_X;
  float const btagWP = BtagHelpers::getBtagWP(BtagHelpers::kDeepCSV_Medium);

  TString const strinputcore = Form("output/GammaTrees/%i/[OUTFILECORE].root", SampleHelpers::theDataYear);
  TString const strinputscalecore = Form("output/GammaTrees/%i/plots/scale_[OUTFILECORE]_MC.root", SampleHelpers::theDataYear);
  TString const stroutputcore = Form("output/GammaTrees/%i/plots/scale_residual_[OUTFILECORE]_MC.root", SampleHelpers::theDataYear);
  gSystem->Exec(Form("mkdir -p output/GammaTrees/%i/plots", SampleHelpers::theDataYear));

  // Establish the variables
  std::vector<Variable> varlist={
    Variable("nak4jets_preselected", "N_{jets}", ExtendedBinning({ 0, 1, 2, 3, 5, 10 })),
    Variable("abs_uPerp", "|u_{perp}| (GeV)", ExtendedBinning({ 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200 })),
    Variable("uParallel", "u_{//} (GeV)", ExtendedBinning({ -500, -150, -100, -50, -30, -20, -10, 0, 20, 70, 200 })),
  };
  auto& xvar = varlist.at(0);
  auto& yvar = varlist.at(1);
  auto& zvar = varlist.at(2);

  // Record MC histograms
  const size_t nSysts=7;
  for (auto const& sample_data:sampleList_Data){
    TString stroutput = stroutputcore;
    HelperFunctions::replaceString<TString, const TString>(stroutput, "[OUTFILECORE]", sample_data);
    MELAout << "Creating output file " << stroutput << endl;
    TFile* foutput = TFile::Open(stroutput, "recreate");

    // Record data histograms
    TH3F hdata = TH3F(
      xvar.name+"_"+yvar.name+"_"+zvar.name+"_Data", "",
      xvar.getNbins(), xvar.getBinning(),
      yvar.getNbins(), yvar.getBinning(),
      zvar.getNbins(), zvar.getBinning()
      );
    hdata.GetXaxis()->SetTitle(xvar.title);
    hdata.GetYaxis()->SetTitle(yvar.title);
    hdata.GetZaxis()->SetTitle(zvar.title);
    {
      TString strinput = strinputcore;
      TChain* tree = new TChain("AnalysisTree");
      HelperFunctions::replaceString<TString, const TString>(strinput, "/[OUTFILECORE].root", "");
      std::vector<TString> inputList = SampleHelpers::lsdir(strinput);
      //MELAout << "Input directories: " << inputList << endl;
      for (auto const& sfile:inputList){
        if (sfile.Contains(sample_data)){
          TString strtmp = Form("%s/%s", strinput.Data(), sfile.Data());
          tree->Add(strtmp);
          MELAout << "Adding " << strtmp << endl;
        }
      }
      DATA_SIMPLE(unsigned int, nak4jets_preselected);
      DATA_SIMPLE(float, uParallel);
      DATA_SIMPLE(float, uPerp);
      DATA_SIMPLE(float, photon_pt);
      int nEntries = tree->GetEntries();
      for (int ev=0; ev<nEntries; ev++){
        tree->GetEntry(ev);
        HelperFunctions::progressbar(ev, nEntries);

        if (photon_pt<130.f) continue;

        xvar.setVal(nak4jets_preselected, 1);
        yvar.setVal(fabs(uPerp), 1);
        zvar.setVal(uParallel, 1);

        hdata.Fill(xvar.val, yvar.val, zvar.val, 1);

        xvar.reset();
        yvar.reset();
        zvar.reset();
      }
      delete tree;

      foutput->WriteTObject(&hdata);
    } // End data

    // Process the MC
    vector<vector<TH3F>> hMClist(sampleList_MC.size(), vector<TH3F>());
    {
      TString strinputscale = strinputscalecore;
      HelperFunctions::replaceString<TString, const TString>(strinputscale, "[OUTFILECORE]", sample_data);
      TFile* finput_scale = TFile::Open(strinputscale, "read");
      std::vector<TH3F const*> hscale;
      hscale.push_back((TH3F const*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_Nominal"));
      hscale.push_back((TH3F const*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_JECup"));
      hscale.push_back((TH3F const*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_JECdn"));
      hscale.push_back((TH3F const*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_JERup"));
      hscale.push_back((TH3F const*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_JERdn"));
      hscale.push_back((TH3F const*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_PUup"));
      hscale.push_back((TH3F const*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_PUdn"));
      hscale.push_back((TH3F const*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_Data"));
      foutput->cd();

      for (size_t is=0; is<sampleList_MC.size(); is++){
        auto& v=hMClist.at(is);
        v.reserve(nSysts);
        for (unsigned int isyst=0; isyst<nSysts; isyst++){
          TString hname = xvar.name + "_" + yvar.name + "_" + sampleList_MC.at(is);
          if (isyst==0) hname += "_Nominal";
          else if (isyst==1) hname += "_JECup";
          else if (isyst==2) hname += "_JECdn";
          else if (isyst==3) hname += "_JERup";
          else if (isyst==4) hname += "_JERdn";
          else if (isyst==5) hname += "_PUup";
          else if (isyst==6) hname += "_PUdn";
          v.emplace_back(
            hname, "",
            xvar.getNbins(), xvar.getBinning(),
            yvar.getNbins(), yvar.getBinning(),
            zvar.getNbins(), zvar.getBinning()
          );
          v.back().Sumw2();
          v.back().GetXaxis()->SetTitle(xvar.title);
          v.back().GetYaxis()->SetTitle(yvar.title);
          v.back().GetZaxis()->SetTitle(zvar.title);
        }
      }
      for (size_t is=0; is<sampleList_MC.size(); is++){
        TString strinput = strinputcore;
        HelperFunctions::replaceString<TString, const TString>(strinput, "[OUTFILECORE]", sampleList_MC.at(is));
        MELAout << "Opening file " << strinput << endl;
        TFile* finput = TFile::Open(strinput, "read");
        TTree* tree = (TTree*) finput->Get("AnalysisTree"); tree->SetBranchStatus("*", 0);

        std::vector<float> genweights(WeightVariables::nWeightTypes, 0);
        for (int iwt=WeightVariables::wCentral; iwt<WeightVariables::nWeightTypes; iwt++){
          TString bname = WeightVariables::getWeightName((WeightVariables::WeightType) iwt);
          tree->SetBranchStatus(bname, 1);
          tree->SetBranchAddress(bname, &(genweights.at(iwt)));
        }
        DATA_SIMPLE(float, xsec);
        DATA_SIMPLE(float, weight_photons);
        DATA_SIMPLE(float, weight_PU);
        DATA_SIMPLE(float, weight_PU_SFUp);
        DATA_SIMPLE(float, weight_PU_SFDn);

        DATA_SIMPLE(float, photon_pt);

        DATA_SIMPLE(unsigned int, nak4jets_preselected);
        DATA_SIMPLE(unsigned int, nak4jets_preselected_JECup);
        DATA_SIMPLE(unsigned int, nak4jets_preselected_JECdn);
        DATA_SIMPLE(unsigned int, nak4jets_preselected_JERup);
        DATA_SIMPLE(unsigned int, nak4jets_preselected_JERdn);

        DATA_SIMPLE(float, uParallel); DATA_SIMPLE(float, uPerp);
        DATA_SIMPLE(float, uParallel_JECup); DATA_SIMPLE(float, uPerp_JECup);
        DATA_SIMPLE(float, uParallel_JECdn); DATA_SIMPLE(float, uPerp_JECdn);
        DATA_SIMPLE(float, uParallel_JERup); DATA_SIMPLE(float, uPerp_JERup);
        DATA_SIMPLE(float, uParallel_JERdn); DATA_SIMPLE(float, uPerp_JERdn);

        DATA_SIMPLE(float, jetHT);
        DATA_SIMPLE(float, jetHT_JECup);
        DATA_SIMPLE(float, jetHT_JECdn);
        DATA_SIMPLE(float, jetHT_JERup);
        DATA_SIMPLE(float, jetHT_JERdn);

        int nEntries = tree->GetEntries();
        bool firstEvent=false;
        std::vector<float> totalWeights_prescaling(nSysts, 0);
        std::vector<float> totalWeights_postscaling(nSysts, 0);
        for (int ev=0; ev<nEntries; ev++){
          tree->GetEntry(ev);
          HelperFunctions::progressbar(ev, nEntries);

          if (photon_pt<130.f) continue;

          for (unsigned int isyst=0; isyst<nSysts; isyst++){
            float totwgt = genweights.at(0)*weight_photons*xsec;
            if (isyst<5) totwgt *= weight_PU;
            else if (isyst==5) totwgt *= weight_PU_SFUp;
            else if (isyst==6) totwgt *= weight_PU_SFDn;
            else assert(0);
            totalWeights_prescaling.at(isyst) += totwgt;

            unsigned int* nak4jets_preselected_ptr = &nak4jets_preselected;
            float* uParallel_ptr = &uParallel;
            float* uPerp_ptr = &uPerp;
            float* jetHT_ptr = &jetHT;
            switch (isyst){
            case 1:
              nak4jets_preselected_ptr = &nak4jets_preselected_JECup;
              uParallel_ptr = &uParallel_JECup;
              uPerp_ptr = &uPerp_JECup;
              jetHT_ptr = &jetHT_JECup;
              break;
            case 2:
              nak4jets_preselected_ptr = &nak4jets_preselected_JECdn;
              uParallel_ptr = &uParallel_JECdn;
              uPerp_ptr = &uPerp_JECdn;
              jetHT_ptr = &jetHT_JECdn;
              break;
            case 3:
              nak4jets_preselected_ptr = &nak4jets_preselected_JERup;
              uParallel_ptr = &uParallel_JERup;
              uPerp_ptr = &uPerp_JERup;
              jetHT_ptr = &jetHT_JERup;
              break;
            case 4:
              nak4jets_preselected_ptr = &nak4jets_preselected_JERdn;
              uParallel_ptr = &uParallel_JERdn;
              uPerp_ptr = &uPerp_JERdn;
              jetHT_ptr = &jetHT_JERdn;
              break;
            }

            if (firstEvent){
              MELAout << "First event, systematic " << isyst << ":\n";
              MELAout << "nak4jets_preselected: " << *nak4jets_preselected_ptr << " (address: " << nak4jets_preselected_ptr << ")" << endl;
              MELAout << "jetHT: " << *jetHT_ptr << " (address: " << jetHT_ptr << ")" << endl;
              MELAout << "uParallel: " << *uParallel_ptr << " (address: " << uParallel_ptr << ")" << endl;
              MELAout << "uPerp: " << *uPerp_ptr << " (address: " << uPerp_ptr << ")" << endl;
            }

            {
              TH3F const*& hscale_MC = hscale.at(isyst);
              int bx_MC = hscale_MC->GetXaxis()->FindBin(*nak4jets_preselected_ptr);
              int by_MC = hscale_MC->GetYaxis()->FindBin(*jetHT_ptr);
              int bz_MC = hscale_MC->GetZaxis()->FindBin(photon_pt);
              if (bx_MC<=0) bx_MC=1;
              else if (bx_MC>hscale_MC->GetNbinsX()) bx_MC=hscale_MC->GetNbinsX();
              if (by_MC<=0) by_MC=1;
              else if (by_MC>hscale_MC->GetNbinsY()) by_MC=hscale_MC->GetNbinsY();
              if (bz_MC<=0) bz_MC=1;
              else if (bz_MC>hscale_MC->GetNbinsZ()) bz_MC=hscale_MC->GetNbinsZ();
              float mcval = hscale_MC->GetBinContent(bx_MC, by_MC, bz_MC);
              if (mcval==0.f){
                MELAerr << "MC val = 0 in bins ( " << bx_MC << ", " << by_MC << ", " << bz_MC << " ) | [";
                MELAerr << *nak4jets_preselected_ptr << ", " << *jetHT_ptr << ", " << photon_pt << " ]" << endl;
              }

              if (firstEvent){
                MELAout << "MC histogram: " << hscale_MC->GetName() << endl;
                MELAout << "MC histogram bins = " << bx_MC << "," << by_MC << "," << bz_MC << endl;
                MELAout << "MC histogram value = " << mcval << endl;
              }

              TH3F const*& hscale_Data = hscale.back();
              int bx_Data = hscale_Data->GetXaxis()->FindBin(*nak4jets_preselected_ptr);
              int by_Data = hscale_Data->GetYaxis()->FindBin(*jetHT_ptr);
              int bz_Data = hscale_Data->GetZaxis()->FindBin(photon_pt);
              if (bx_Data<=0) bx_Data=1;
              else if (bx_Data>hscale_Data->GetNbinsX()) bx_Data=hscale_Data->GetNbinsX();
              if (by_Data<=0) by_Data=1;
              else if (by_Data>hscale_Data->GetNbinsY()) by_Data=hscale_Data->GetNbinsY();
              if (bz_Data<=0) bz_Data=1;
              else if (bz_Data>hscale_Data->GetNbinsZ()) bz_Data=hscale_Data->GetNbinsZ();
              float dataval = hscale_Data->GetBinContent(bx_Data, by_Data, bz_Data);

              if (firstEvent){
                MELAout << "Data histogram: " << hscale_Data->GetName() << endl;
                MELAout << "Data histogram bins = " << bx_Data << "," << by_Data << "," << bz_Data << endl;
                MELAout << "Data histogram value = " << dataval << endl;
              }

              totwgt *= (mcval<=0.f ? 0.f : dataval/mcval);
            }

            xvar.setVal(*nak4jets_preselected_ptr, 1);
            yvar.setVal(fabs(*uPerp_ptr), 1);
            zvar.setVal(*uParallel_ptr, 1);

            if (firstEvent){
              MELAout << "Filling " << hMClist.at(is).at(isyst).GetName() << endl;
              MELAout << "\t- Fill variables: " << xvar.val << ", " << yvar.val << ", " << zvar.val << endl;
              MELAout << "\t- Fill weight: " << totwgt << endl;
            }

            hMClist.at(is).at(isyst).Fill(xvar.val, yvar.val, zvar.val, totwgt);
            totalWeights_postscaling.at(isyst) += totwgt;

            xvar.reset();
            yvar.reset();
            zvar.reset();
          }
          if (firstEvent) firstEvent=false;
        }
        // Scale MC histograms
        {
          TTree* tree_md = (TTree*) finput->Get("metadata");
          float /*sumWgts=0, */sumWgts_PU=0, sumWgts_PU_SFUp=0, sumWgts_PU_SFDn=0;
          TString bname = WeightVariables::getWeightName(WeightVariables::wCentral);
          //tree_md->SetBranchAddress(bname, &sumWgts);
          tree_md->SetBranchAddress(bname+"_PU", &sumWgts_PU);
          tree_md->SetBranchAddress(bname+"_PU_SFUp", &sumWgts_PU_SFUp);
          tree_md->SetBranchAddress(bname+"_PU_SFDn", &sumWgts_PU_SFDn);
          tree_md->GetEntry(0);
          if (sumWgts_PU<=0.f) MELAout << "ERROR: sumWgts_PU = " << sumWgts_PU << endl;
          else{
            for (unsigned int isyst=0; isyst<nSysts; isyst++){
              auto& h = hMClist.at(is).at(isyst);
              float* sumWgts_ptr = &sumWgts_PU;
              if (isyst==5) sumWgts_ptr = &sumWgts_PU_SFUp;
              else if (isyst==6) sumWgts_ptr = &sumWgts_PU_SFDn;
              h.Scale(1./(*sumWgts_ptr));
              MELAout << "Scaled " << h.GetName() << " by 1 / " << (*sumWgts_ptr) << endl;
              totalWeights_prescaling.at(isyst) /= (*sumWgts_ptr);
              totalWeights_postscaling.at(isyst) /= (*sumWgts_ptr);
            }
          }
        }

        finput->Close();

        MELAout << "\t- Sum of pre-scaling weights: " << totalWeights_prescaling << endl;
        MELAout << "\t- Sum of post-scaling weights: " << totalWeights_postscaling << endl;
      }

      finput_scale->Close();
    } // End MC processing
    { // Get total MC
      float dataYield = hdata.Integral(1, hdata.GetNbinsX(), 1, hdata.GetNbinsY(), 1, hdata.GetNbinsZ());
      for (unsigned int isyst=0; isyst<nSysts; isyst++){
        TString strappend;
        if (isyst==0) strappend="_Nominal";
        else if (isyst==1) strappend="_JECup";
        else if (isyst==2) strappend="_JECdn";
        else if (isyst==3) strappend="_JERup";
        else if (isyst==4) strappend="_JERdn";
        else if (isyst==5) strappend="_PUup";
        else if (isyst==6) strappend="_PUdn";
        TH3F* htotMC_syst = (TH3F*) hdata.Clone(xvar.name+"_"+yvar.name+"_"+zvar.name+"_totMC"+strappend);
        htotMC_syst->Reset("ICESM");
        //htotMC_syst->Sumw2();
        for (auto const& v:hMClist) htotMC_syst->Add(&(v.at(isyst)));
        float MCyield = htotMC_syst->Integral(1, htotMC_syst->GetNbinsX(), 1, htotMC_syst->GetNbinsY(), 1, htotMC_syst->GetNbinsZ());
        float data_MC_scale = dataYield/MCyield;
        MELAout << "Data/MC = " << dataYield << " / " << MCyield << " = " << data_MC_scale << endl;
        //htotMC_syst->Scale(data_MC_scale);
        foutput->WriteTObject(htotMC_syst);

        // Get conditional histograms in nak4jets_preselected
        HelperFunctions::conditionalizeHistogram<TH3F>(htotMC_syst, 0, nullptr, false); htotMC_syst->SetName(Form("%s_Conditional", htotMC_syst->GetName()));
        foutput->WriteTObject(htotMC_syst);

        delete htotMC_syst;
      }
      HelperFunctions::conditionalizeHistogram<TH3F>(&hdata, 0, nullptr, false); hdata.SetName(Form("%s_Conditional", hdata.GetName()));
      foutput->WriteTObject(&hdata);
    }

    foutput->Close();
  }
}

void getMCHistograms_1D(
  TString const& strinputcore, std::vector<TString> const& sampleList_MC, std::vector<Variable>& varlist, std::vector<std::vector<std::vector<TH1F>>>& hMClist,
  std::vector<TH3F*>* hscale=nullptr,
  std::vector<TH3F*>* hresscale=nullptr,
  TTree* tin=nullptr,
  METCorrectionHandler* metCorrector=nullptr
){
  const size_t nSysts = 7;
  if (hscale) assert(hscale.size()==nSysts+1);
  if (hresscale) assert(hresscale.size()==nSysts+1);
  hMClist = std::vector<std::vector<std::vector<TH1F>>>(sampleList_MC.size(), vector<vector<TH1F>>(nSysts, vector<TH1F>()));

  for (size_t is=0; is<sampleList_MC.size(); is++){
    auto& vv=hMClist.at(is);
    for (unsigned int isyst=0; isyst<nSysts; isyst++){
      auto& v=vv.at(isyst);
      v.reserve(varlist.size());
      for (auto& var:varlist){
        TString hname = var.name + "_" + sampleList_MC.at(is);
        if (isyst==0) hname += "_Nominal";
        else if (isyst==1) hname += "_JECup";
        else if (isyst==2) hname += "_JECdn";
        else if (isyst==3) hname += "_JERup";
        else if (isyst==4) hname += "_JERdn";
        else if (isyst==5) hname += "_PUup";
        else if (isyst==6) hname += "_PUdn";
        v.emplace_back(hname, "", var.getNbins(), var.getBinning());
        v.back().Sumw2();
        v.back().GetXaxis()->SetTitle(var.title);
        v.back().GetYaxis()->SetTitle("Events");
      }
    }
  }
  for (size_t is=0; is<sampleList_MC.size(); is++){
    TString strinput = strinputcore;
    HelperFunctions::replaceString<TString, const TString>(strinput, "[OUTFILECORE]", sampleList_MC.at(is));
    MELAout << "Opening file " << strinput << endl;
    TFile* finput = TFile::Open(strinput, "read");
    TTree* tree = (TTree*) finput->Get("AnalysisTree"); tree->SetBranchStatus("*", 0);

    // Scale MC histograms
    vector<float> overallWeight(nSysts, 1);
    {
      TTree* tree_md = (TTree*) finput->Get("metadata");
      float /*sumWgts=0, */sumWgts_PU=0, sumWgts_PU_SFUp=0, sumWgts_PU_SFDn=0;
      TString bname = WeightVariables::getWeightName(WeightVariables::wCentral);
      //tree_md->SetBranchAddress(bname, &sumWgts);
      tree_md->SetBranchAddress(bname+"_PU", &sumWgts_PU);
      tree_md->SetBranchAddress(bname+"_PU_SFUp", &sumWgts_PU_SFUp);
      tree_md->SetBranchAddress(bname+"_PU_SFDn", &sumWgts_PU_SFDn);
      tree_md->GetEntry(0);
      if (sumWgts_PU<=0.f) MELAout << "ERROR: sumWgts_PU = " << sumWgts_PU << endl;
      else{
        for (unsigned int isyst=0; isyst<nSysts; isyst++){
          float* sumWgts_ptr = &sumWgts_PU;
          if (isyst==5) sumWgts_ptr = &sumWgts_PU_SFUp;
          else if (isyst==6) sumWgts_ptr = &sumWgts_PU_SFDn;
          overallWeight.at(isyst) = 1./(*sumWgts_ptr);
        }
      }
    }

    std::vector<float> genweights(WeightVariables::nWeightTypes, 0);
    for (int iwt=WeightVariables::wCentral; iwt<WeightVariables::nWeightTypes; iwt++){
      TString bname = WeightVariables::getWeightName((WeightVariables::WeightType) iwt);
      tree->SetBranchStatus(bname, 1);
      tree->SetBranchAddress(bname, &(genweights.at(iwt)));
    }
    DATA_SIMPLE(float, xsec);
    DATA_SIMPLE(float, weight_photons);
    DATA_SIMPLE(float, weight_PU);
    DATA_SIMPLE(float, weight_PU_SFUp);
    DATA_SIMPLE(float, weight_PU_SFDn);

    DATA_SIMPLE(float, genMET);
    DATA_SIMPLE(float, genMETPhi);

    DATA_SIMPLE(float, photon_pt);

    DATA_SIMPLE(unsigned int, nak4jets_preselected);
    DATA_SIMPLE(unsigned int, nak4jets_preselected_JECup);
    DATA_SIMPLE(unsigned int, nak4jets_preselected_JECdn);
    DATA_SIMPLE(unsigned int, nak4jets_preselected_JERup);
    DATA_SIMPLE(unsigned int, nak4jets_preselected_JERdn);

    DATA_SIMPLE(float, uParallel); DATA_SIMPLE(float, uPerp);
    DATA_SIMPLE(float, uParallel_JECup); DATA_SIMPLE(float, uPerp_JECup);
    DATA_SIMPLE(float, uParallel_JECdn); DATA_SIMPLE(float, uPerp_JECdn);
    DATA_SIMPLE(float, uParallel_JERup); DATA_SIMPLE(float, uPerp_JERup);
    DATA_SIMPLE(float, uParallel_JERdn); DATA_SIMPLE(float, uPerp_JERdn);

    DATA_SIMPLE(float, jetHT);
    DATA_SIMPLE(float, jetHT_JECup);
    DATA_SIMPLE(float, jetHT_JECdn);
    DATA_SIMPLE(float, jetHT_JERup);
    DATA_SIMPLE(float, jetHT_JERdn);

    DATA_SIMPLE(float, pfmet); DATA_SIMPLE(float, pfmet_JECup); DATA_SIMPLE(float, pfmet_JECdn);
    DATA_SIMPLE(float, pfmetPhi); DATA_SIMPLE(float, pfmetPhi_JECup); DATA_SIMPLE(float, pfmetPhi_JECdn);
    DATA_SIMPLE(float, MET_Parallel); DATA_SIMPLE(float, MET_Perp);
    DATA_SIMPLE(float, MET_Parallel_JECup); DATA_SIMPLE(float, MET_Perp_JECup);
    DATA_SIMPLE(float, MET_Parallel_JECdn); DATA_SIMPLE(float, MET_Perp_JECdn);

    float totwgt_Nominal, totwgt_JECup, totwgt_JECdn, totwgt_JERup, totwgt_JERdn, totwgt_PUup, totwgt_PUdn;
    if (tin && is==0){
      tin->Branch("nak4jets_preselected", &nak4jets_preselected);
      tin->Branch("nak4jets_preselected_JECup", &nak4jets_preselected_JECup); tin->Branch("nak4jets_preselected_JECdn", &nak4jets_preselected_JECdn);
      tin->Branch("nak4jets_preselected_JERup", &nak4jets_preselected_JERup); tin->Branch("nak4jets_preselected_JERdn", &nak4jets_preselected_JERdn);
      tin->Branch("MET_Parallel", &MET_Parallel);
      tin->Branch("MET_Parallel_JECup", &MET_Parallel_JECup); tin->Branch("MET_Parallel_JECdn", &MET_Parallel_JECdn);
      tin->Branch("MET_Perp", &MET_Perp);
      tin->Branch("MET_Perp_JECup", &MET_Perp_JECup); tin->Branch("MET_Perp_JECdn", &MET_Perp_JECdn);
      tin->Branch("totwgt", &totwgt_Nominal);
      tin->Branch("totwgt_JECup", &totwgt_JECup); tin->Branch("totwgt_JECdn", &totwgt_JECdn);
      tin->Branch("totwgt_JERup", &totwgt_JERup); tin->Branch("totwgt_JERdn", &totwgt_JERdn);
      tin->Branch("totwgt_PUup", &totwgt_PUup); tin->Branch("totwgt_PUdn", &totwgt_PUdn);
    }
    else if (tin){
      tin->ResetBranchAddresses();
      tin->SetBranchAddress("nak4jets_preselected", &nak4jets_preselected);
      tin->SetBranchAddress("nak4jets_preselected_JECup", &nak4jets_preselected_JECup); tin->SetBranchAddress("nak4jets_preselected_JECdn", &nak4jets_preselected_JECdn);
      tin->SetBranchAddress("nak4jets_preselected_JERup", &nak4jets_preselected_JERup); tin->SetBranchAddress("nak4jets_preselected_JERdn", &nak4jets_preselected_JERdn);
      tin->SetBranchAddress("MET_Parallel", &MET_Parallel);
      tin->SetBranchAddress("MET_Parallel_JECup", &MET_Parallel_JECup); tin->SetBranchAddress("MET_Parallel_JECdn", &MET_Parallel_JECdn);
      tin->SetBranchAddress("MET_Perp", &MET_Perp);
      tin->SetBranchAddress("MET_Perp_JECup", &MET_Perp_JECup); tin->SetBranchAddress("MET_Perp_JECdn", &MET_Perp_JECdn);
      tin->SetBranchAddress("totwgt", &totwgt_Nominal);
      tin->SetBranchAddress("totwgt_JECup", &totwgt_JECup); tin->SetBranchAddress("totwgt_JECdn", &totwgt_JECdn);
      tin->SetBranchAddress("totwgt_JERup", &totwgt_JERup); tin->SetBranchAddress("totwgt_JERdn", &totwgt_JERdn);
      tin->SetBranchAddress("totwgt_PUup", &totwgt_PUup); tin->SetBranchAddress("totwgt_PUdn", &totwgt_PUdn);
    }

    int nEntries = tree->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      tree->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      float photon_phi = pfmetPhi - atan2(MET_Perp, MET_Parallel);

      // If a corrector is available, correct the MET before looping over systematics
      METObject metobj;
      metobj.extras.met
        = metobj.extras.met_METup = metobj.extras.met_METdn
        = metobj.extras.met_JERup = metobj.extras.met_JERdn
        = metobj.extras.met_PUup = metobj.extras.met_PUdn
        = pfmet;
      metobj.extras.met_JECup = pfmet_JECup;
      metobj.extras.met_JECdn = pfmet_JECdn;
      metobj.extras.phi
        = metobj.extras.phi_METup = metobj.extras.phi_METdn
        = metobj.extras.phi_JERup = metobj.extras.phi_JERdn
        = metobj.extras.phi_PUup = metobj.extras.phi_PUdn
        = pfmetPhi;
      metobj.extras.phi_JECup = pfmetPhi_JECup;
      metobj.extras.phi_JECdn = pfmetPhi_JECdn;
      if (metCorrector) metCorrector->correctMET(genMET, genMETPhi, &metobj, false);

      for (unsigned int isyst=0; isyst<nSysts; isyst++){
        float totwgt = genweights.at(0)*weight_photons*xsec*overallWeight.at(isyst);
        if (isyst<5) totwgt *= weight_PU;
        else if (isyst==5) totwgt *= weight_PU_SFUp;
        else if (isyst==6) totwgt *= weight_PU_SFDn;
        else assert(0);

        unsigned int* nak4jets_preselected_ptr = &nak4jets_preselected;
        float* uParallel_ptr = &uParallel;
        float* uPerp_ptr = &uPerp;
        float* jetHT_ptr = &jetHT;
        float* pfmet_ptr = &pfmet;
        float* pfmet_corrected_ptr = &metobj.extras.met;
        float* pfmetPhi_corrected_ptr = &metobj.extras.phi;
        float* MET_Parallel_ptr = &MET_Parallel;
        float* MET_Perp_ptr = &MET_Perp;
        switch (isyst){
        case 1:
          nak4jets_preselected_ptr = &nak4jets_preselected_JECup;
          uParallel_ptr = &uParallel_JECup;
          uPerp_ptr = &uPerp_JECup;
          jetHT_ptr = &jetHT_JECup;
          pfmet_ptr = &pfmet_JECup;
          pfmet_corrected_ptr = &metobj.extras.met_JECup;
          pfmetPhi_corrected_ptr = &metobj.extras.phi_JECup;
          MET_Parallel_ptr = &MET_Parallel_JECup;
          MET_Perp_ptr = &MET_Perp_JECup;
          break;
        case 2:
          nak4jets_preselected_ptr = &nak4jets_preselected_JECdn;
          uParallel_ptr = &uParallel_JECdn;
          uPerp_ptr = &uPerp_JECdn;
          jetHT_ptr = &jetHT_JECdn;
          pfmet_ptr = &pfmet_JECdn;
          pfmet_corrected_ptr = &metobj.extras.met_JECdn;
          pfmetPhi_corrected_ptr = &metobj.extras.phi_JECdn;
          MET_Parallel_ptr = &MET_Parallel_JECdn;
          MET_Perp_ptr = &MET_Perp_JECdn;
          break;
        case 3:
          nak4jets_preselected_ptr = &nak4jets_preselected_JERup;
          uParallel_ptr = &uParallel_JERup;
          uPerp_ptr = &uPerp_JERup;
          jetHT_ptr = &jetHT_JERup;
          pfmet_corrected_ptr = &metobj.extras.met_JERup;
          pfmetPhi_corrected_ptr = &metobj.extras.phi_JERup;
          break;
        case 4:
          nak4jets_preselected_ptr = &nak4jets_preselected_JERdn;
          uParallel_ptr = &uParallel_JERdn;
          uPerp_ptr = &uPerp_JERdn;
          jetHT_ptr = &jetHT_JERdn;
          pfmet_corrected_ptr = &metobj.extras.met_JERdn;
          pfmetPhi_corrected_ptr = &metobj.extras.phi_JERdn;
          break;
        case 5:
          pfmet_corrected_ptr = &metobj.extras.met_PUup;
          pfmetPhi_corrected_ptr = &metobj.extras.phi_PUup;
          break;
        case 6:
          pfmet_corrected_ptr = &metobj.extras.met_PUdn;
          pfmetPhi_corrected_ptr = &metobj.extras.phi_PUdn;
          break;
        }

        float MET_Perp_corrected=0;
        float MET_Parallel_corrected=0;
        get2DParallelAndPerpendicularComponents(
          TVector3(photon_pt*cos(photon_phi), photon_pt*sin(photon_phi), 0),
          TVector3((*pfmet_corrected_ptr)*cos(*pfmetPhi_corrected_ptr), (*pfmet_corrected_ptr)*sin(*pfmetPhi_corrected_ptr), 0),
          MET_Parallel_corrected, MET_Perp_corrected
        );

        if (hscale){
          int bx_MC = hscale->at(isyst)->GetXaxis()->FindBin(*nak4jets_preselected_ptr);
          int by_MC = hscale->at(isyst)->GetYaxis()->FindBin(*jetHT_ptr);
          int bz_MC = hscale->at(isyst)->GetZaxis()->FindBin(photon_pt);
          if (bx_MC<=0) bx_MC=1;
          else if (bx_MC>hscale->at(isyst)->GetNbinsX()) bx_MC=hscale->at(isyst)->GetNbinsX();
          if (by_MC<=0) by_MC=1;
          else if (by_MC>hscale->at(isyst)->GetNbinsY()) by_MC=hscale->at(isyst)->GetNbinsY();
          if (bz_MC<=0) bz_MC=1;
          else if (bz_MC>hscale->at(isyst)->GetNbinsZ()) bz_MC=hscale->at(isyst)->GetNbinsZ();
          float mcval = hscale->at(isyst)->GetBinContent(bx_MC, by_MC, bz_MC);

          int bx_Data = hscale->back()->GetXaxis()->FindBin(*nak4jets_preselected_ptr);
          int by_Data = hscale->back()->GetYaxis()->FindBin(*jetHT_ptr);
          int bz_Data = hscale->back()->GetZaxis()->FindBin(photon_pt);
          if (bx_Data<=0) bx_Data=1;
          else if (bx_Data>hscale->back()->GetNbinsX()) bx_Data=hscale->back()->GetNbinsX();
          if (by_Data<=0) by_Data=1;
          else if (by_Data>hscale->back()->GetNbinsY()) by_Data=hscale->back()->GetNbinsY();
          if (bz_Data<=0) bz_Data=1;
          else if (bz_Data>hscale->back()->GetNbinsZ()) bz_Data=hscale->back()->GetNbinsZ();
          float dataval = hscale->back()->GetBinContent(bx_Data, by_Data, bz_Data);

          totwgt *= (mcval<=0.f ? 0.f : dataval/mcval);
        }
        if (hresscale){
          int bx_MC = hresscale->at(isyst)->GetXaxis()->FindBin(*nak4jets_preselected_ptr);
          int by_MC = hresscale->at(isyst)->GetYaxis()->FindBin(fabs(*uPerp_ptr));
          int bz_MC = hresscale->at(isyst)->GetZaxis()->FindBin(*uParallel_ptr);
          if (bx_MC<=0) bx_MC=1;
          else if (bx_MC>hresscale->at(isyst)->GetNbinsX()) bx_MC=hresscale->at(isyst)->GetNbinsX();
          if (by_MC<=0) by_MC=1;
          else if (by_MC>hresscale->at(isyst)->GetNbinsY()) by_MC=hresscale->at(isyst)->GetNbinsY();
          if (bz_MC<=0) bz_MC=1;
          else if (bz_MC>hresscale->at(isyst)->GetNbinsZ()) bz_MC=hresscale->at(isyst)->GetNbinsZ();
          float mcval = hresscale->at(isyst)->GetBinContent(bx_MC, by_MC, bz_MC);

          int bx_Data = hresscale->back()->GetXaxis()->FindBin(*nak4jets_preselected_ptr);
          int by_Data = hresscale->back()->GetYaxis()->FindBin(fabs(*uPerp_ptr));
          int bz_Data = hresscale->back()->GetZaxis()->FindBin(*uParallel_ptr);
          if (bx_Data<=0) bx_Data=1;
          else if (bx_Data>hresscale->back()->GetNbinsX()) bx_Data=hresscale->back()->GetNbinsX();
          if (by_Data<=0) by_Data=1;
          else if (by_Data>hresscale->back()->GetNbinsY()) by_Data=hresscale->back()->GetNbinsY();
          if (bz_Data<=0) bz_Data=1;
          else if (bz_Data>hresscale->back()->GetNbinsZ()) bz_Data=hresscale->back()->GetNbinsZ();
          float dataval = hresscale->back()->GetBinContent(bx_Data, by_Data, bz_Data);

          totwgt *= (mcval<=0.f ? 0.f : dataval/mcval);
        }

        switch (isyst){
        case 1:
          totwgt_JECup = totwgt;
          break;
        case 2:
          totwgt_JECdn = totwgt;
          break;
        case 3:
          totwgt_JERup = totwgt;
          break;
        case 4:
          totwgt_JERdn = totwgt;
          break;
        case 5:
          totwgt_PUup = totwgt;
          break;
        case 6:
          totwgt_PUdn = totwgt;
          break;
        default:
          totwgt_Nominal = totwgt;
          break;
        }

        for (size_t iv=0; iv<varlist.size(); iv++){
          auto& var = varlist.at(iv);
          if (var.name == "nak4jets_preselected") var.setVal(*nak4jets_preselected_ptr, 1);
          else{
            if (*nak4jets_preselected_ptr==0) continue;
            if (var.name == "photon_pt") var.setVal(photon_pt, 1);
            else{
              if (photon_pt<130.f) continue;
              if (var.name == "pfmet") var.setVal(*pfmet_ptr, 1);
              else if (var.name == "pfmet_corrected") var.setVal(*pfmet_corrected_ptr, 1);
              else if (var.name == "MET_Parallel") var.setVal(*MET_Parallel_ptr, 1);
              else if (var.name == "MET_Perp") var.setVal(*MET_Perp_ptr, 1);
              else if (var.name == "MET_Parallel_corrected") var.setVal(MET_Parallel_corrected, 1);
              else if (var.name == "MET_Perp_corrected") var.setVal(MET_Perp_corrected, 1);
              else if (var.name == "uParallel") var.setVal(*uParallel_ptr, 1);
              else if (var.name == "uPerp") var.setVal(*uPerp_ptr, 1);
              else if (var.name == "uT") var.setVal(sqrt(pow(*uPerp_ptr, 2)+pow(*uParallel_ptr, 2)), 1);
              else if (var.name == "jetHT") var.setVal(*jetHT_ptr, 1);
              else if (var.name == "totParallel") var.setVal(photon_pt + (*uParallel_ptr) + (*MET_Parallel_ptr), 1);
              else if (var.name == "totPerp") var.setVal((*uPerp_ptr) + (*MET_Perp_ptr), 1);
            }
          }
          hMClist.at(is).at(isyst).at(iv).Fill(var.val, totwgt);
          var.reset();
        }
      }

      if (tin && photon_pt>=130.f) tin->Fill();
    }

    finput->Close();
  }
}


void plotGammaTrees_DataMC(){
  std::vector<TString> sampleList_Data={
    "Run2017B-31Mar2018-v1",
    "Run2017C-31Mar2018-v1",
    "Run2017D-31Mar2018-v1",
    "Run2017E-31Mar2018-v1",
    "Run2017F-31Mar2018-v1",
    "Run2017F-09May2018-v1"
  };
  std::vector<TString> sampleList_MC={
    "GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",

    "QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8",

    "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    "TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8",
    //"TTGamma_Dilept_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromT_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromTbar_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8",

    "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WZG_TuneCP5_13TeV-amcatnlo-pythia8",

    "ZJetsToNuNu_HT-100To200_13TeV-madgraph",
    "ZJetsToNuNu_HT-200To400_13TeV-madgraph",
    "ZJetsToNuNu_HT-400To600_13TeV-madgraph",
    "ZJetsToNuNu_HT-600To800_13TeV-madgraph",
    "ZJetsToNuNu_HT-800To1200_13TeV-madgraph",
    "ZJetsToNuNu_HT-1200To2500_13TeV-madgraph",
    "ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph"
  };

  SampleHelpers::theDataYear=2017;
  SampleHelpers::theDataPeriod="2017";
  SampleHelpers::theDataVersion=SampleHelpers::kCMSSW_9_4_X;
  float const btagWP = BtagHelpers::getBtagWP(BtagHelpers::kDeepCSV_Medium);

  TString const strinputcore = Form("output/GammaTrees/%i/[OUTFILECORE].root", SampleHelpers::theDataYear);
  TString const stroutputcore = Form("output/GammaTrees/%i/plots/compare_[OUTFILECORE]_MC.root", SampleHelpers::theDataYear);
  gSystem->Exec(Form("mkdir -p output/GammaTrees/%i/plots", SampleHelpers::theDataYear));

  // Establish the variables
  std::vector<Variable> varlist={
    Variable("nak4jets_preselected", "N_{jets}", 10, 0, 10),
    Variable("pfmet", "E_{T}^{miss} (GeV)", 100, 0, 200),
    Variable("photon_pt", "p_{T}^{#gamma} (GeV)", 70, 50, 400),
    Variable("uParallel", "u_{//} (GeV)", 100, -500, 200),
    Variable("uPerp", "u_{#perp} (GeV)", 100, -200, 200),
    Variable("totParallel", "p_{T}^{#gamma}+u_{//}+E_{T,//}^{miss} (GeV)", 100, -200, 200),
    Variable("totPerp", "u_{#perp}+E_{T,#perp}^{miss} (GeV)", 100, -200, 200)
  };

  const unsigned int nSysts=5;
  // Record MC histograms
  vector<vector<vector<TH1F>>> hMClist;
  getMCHistograms_1D(strinputcore, sampleList_MC, varlist, hMClist);

  for (auto const& sample_data:sampleList_Data){
    TString stroutput = stroutputcore;
    HelperFunctions::replaceString<TString, const TString>(stroutput, "[OUTFILECORE]", sample_data);
    MELAout << "Creating output file " << stroutput << endl;
    TFile* foutput = TFile::Open(stroutput, "recreate");

    // Record data histograms
    vector<TH1F> hdatalist; hdatalist.reserve(varlist.size());
    for (auto const& var:varlist){
      hdatalist.emplace_back(
        var.name+"_Data", "",
        var.getNbins(), var.getBinning()
      );
      hdatalist.back().GetXaxis()->SetTitle(var.title);
      hdatalist.back().GetYaxis()->SetTitle("Events");
    }
    {
      TString strinput = strinputcore;
      TChain* tree = new TChain("AnalysisTree");
      HelperFunctions::replaceString<TString, const TString>(strinput, "/[OUTFILECORE].root", "");
      std::vector<TString> inputList = SampleHelpers::lsdir(strinput);
      //MELAout << "Input directories: " << inputList << endl;
      for (auto const& sfile:inputList){
        if (sfile.Contains(sample_data)){
          TString strtmp = Form("%s/%s", strinput.Data(), sfile.Data());
          tree->Add(strtmp);
          MELAout << "Adding " << strtmp << endl;
        }
      }
      DATA_SIMPLE(unsigned int, nak4jets_preselected);
      DATA_SIMPLE(float, pfmet);
      //DATA_SIMPLE(float, pfmetPhi);
      //DATA_SIMPLE(float, genJetHT);
      //DATA_SIMPLE(float, jetHT);
      //DATA_SIMPLE(float, ak4jet_leadingpt_pt);
      DATA_SIMPLE(float, photon_pt);
      //DATA_SIMPLE(float, photon_eta);
      //DATA_SIMPLE(float, photon_phi);
      //DATA_SIMPLE(float, photon_mass);
      DATA_SIMPLE(float, uParallel);
      DATA_SIMPLE(float, uPerp);
      DATA_SIMPLE(float, MET_Parallel);
      DATA_SIMPLE(float, MET_Perp);
      int nEntries = tree->GetEntries();
      for (int ev=0; ev<nEntries; ev++){
        tree->GetEntry(ev);
        HelperFunctions::progressbar(ev, nEntries);

        for (size_t iv=0; iv<varlist.size(); iv++){
          auto& var = varlist.at(iv);
          if (var.name == "nak4jets_preselected") var.setVal(nak4jets_preselected, 1);
          else{
            if (nak4jets_preselected==0) continue;
            if (var.name == "pfmet") var.setVal(pfmet, 1);
            else if (var.name == "photon_pt") var.setVal(photon_pt, 1);
            else if (var.name == "uParallel") var.setVal(uParallel, 1);
            else if (var.name == "uPerp") var.setVal(uPerp, 1);
            else if (var.name == "totParallel") var.setVal(photon_pt+uParallel+MET_Parallel, 1);
            else if (var.name == "totPerp") var.setVal(uPerp+MET_Perp, 1);
          }

          hdatalist.at(iv).Fill(var.val, var.weights.at(0));
          var.reset();
        }
      }

      delete tree;
    }
    for (size_t iv=0; iv<varlist.size();iv++){
      auto& var = varlist.at(iv);

      auto const& hdata=hdatalist.at(iv);
      foutput->WriteTObject(&hdata);

      for (unsigned int isyst=0; isyst<nSysts; isyst++){
        TString strappend;
        if (isyst==0) strappend="_Nominal";
        else if (isyst==1) strappend="_JECup";
        else if (isyst==2) strappend="_JECdn";
        TH1F* htotMC_syst = (TH1F*) hdata.Clone(var.name+"_totMC"+strappend);
        htotMC_syst->Reset("ICESM");
        htotMC_syst->Sumw2();
        for (auto const& vs:hMClist) htotMC_syst->Add(&(vs.at(isyst).at(iv)));
        float dataYield = hdata.Integral(1, hdata.GetNbinsX());
        float MCyield = htotMC_syst->Integral(1, htotMC_syst->GetNbinsX());
        float data_MC_scale = dataYield/MCyield;
        MELAout << "Data/MC = " << dataYield << " / " << MCyield << " = " << data_MC_scale << endl;
        htotMC_syst->Scale(data_MC_scale);
        foutput->WriteTObject(htotMC_syst);
        delete htotMC_syst;
      }
    }

    foutput->Close();
  }
}

void plotGammaTrees_Scaled_DataMC(){
  std::vector<TString> sampleList_Data={
    "Run2017B-31Mar2018-v1",
    "Run2017C-31Mar2018-v1",
    "Run2017D-31Mar2018-v1",
    "Run2017E-31Mar2018-v1",
    "Run2017F-31Mar2018-v1",
    "Run2017F-09May2018-v1"
  };
  std::vector<TString> sampleList_MC={
    "GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",

    "QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8",

    "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    "TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8",
    //"TTGamma_Dilept_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromT_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromTbar_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8",

    "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WZG_TuneCP5_13TeV-amcatnlo-pythia8",

    "ZJetsToNuNu_HT-100To200_13TeV-madgraph",
    "ZJetsToNuNu_HT-200To400_13TeV-madgraph",
    "ZJetsToNuNu_HT-400To600_13TeV-madgraph",
    "ZJetsToNuNu_HT-600To800_13TeV-madgraph",
    "ZJetsToNuNu_HT-800To1200_13TeV-madgraph",
    "ZJetsToNuNu_HT-1200To2500_13TeV-madgraph",
    "ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph"
  };

  SampleHelpers::theDataYear=2017;
  SampleHelpers::theDataPeriod="2017";
  SampleHelpers::theDataVersion=SampleHelpers::kCMSSW_9_4_X;
  float const btagWP = BtagHelpers::getBtagWP(BtagHelpers::kDeepCSV_Medium);

  TString const strinputcore = Form("output/GammaTrees/%i/[OUTFILECORE].root", SampleHelpers::theDataYear);
  TString const strinputscalecore = Form("output/GammaTrees/%i/plots/scale_[OUTFILECORE]_MC.root", SampleHelpers::theDataYear);
  TString const strinputresscalecore = Form("output/GammaTrees/%i/plots/scale_residual_[OUTFILECORE]_MC.root", SampleHelpers::theDataYear);
  TString const stroutputcore = Form("output/GammaTrees/%i/plots/compare_scaled_[OUTFILECORE]_MC.root", SampleHelpers::theDataYear);
  gSystem->Exec(Form("mkdir -p output/GammaTrees/%i/plots", SampleHelpers::theDataYear));

  // Establish the variables
  std::vector<Variable> varlist={
    Variable("nak4jets_preselected", "N_{jets}", 10, 0, 10),
    Variable("pfmet", "E_{T}^{miss} (GeV)", 100, 0, 200),
    Variable("MET_Parallel", "E_{T,//}^{miss} (GeV)", 150, -300, 300),
    Variable("MET_Perp", "E_{T,perp}^{miss} (GeV)", 100, -200, 200),
    Variable("photon_pt", "p_{T}^{#gamma} (GeV)", 70, 50, 400),
    Variable("uParallel", "u_{//} (GeV)", 100, -500, 200),
    Variable("uPerp", "u_{perp} (GeV)", 100, -200, 200),
    Variable("uT", "u_{T} (GeV)", 100, 0, 600),
    Variable("jetHT", "Jet H_{T} (GeV)", 100, 0, 1000),
    Variable("totParallel", "p_{T}^{#gamma}+u_{//}+E_{T,//}^{miss} (GeV)", 100, -200, 200),
    Variable("totPerp", "u_{perp}+E_{T,perp}^{miss} (GeV)", 100, -200, 200)
  };

  for (auto const& sample_data:sampleList_Data){
    TString stroutput = stroutputcore;
    HelperFunctions::replaceString<TString, const TString>(stroutput, "[OUTFILECORE]", sample_data);
    MELAout << "Creating output file " << stroutput << endl;
    TFile* foutput = TFile::Open(stroutput, "recreate");

    TTree tout_MC("FinalTree_MC", "");
    TTree tout_Data("FinalTree_Data", "");

    // Record data histograms
    vector<TH1F> hdatalist; hdatalist.reserve(varlist.size());
    for (auto const& var:varlist){
      hdatalist.emplace_back(
        var.name+"_Data", "",
        var.getNbins(), var.getBinning()
      );
      hdatalist.back().GetXaxis()->SetTitle(var.title);
      hdatalist.back().GetYaxis()->SetTitle("Events");
    }
    {
      TString strinput = strinputcore;
      TChain* tree = new TChain("AnalysisTree");
      HelperFunctions::replaceString<TString, const TString>(strinput, "/[OUTFILECORE].root", "");
      std::vector<TString> inputList = SampleHelpers::lsdir(strinput);
      //MELAout << "Input directories: " << inputList << endl;
      for (auto const& sfile:inputList){
        if (sfile.Contains(sample_data)){
          TString strtmp = Form("%s/%s", strinput.Data(), sfile.Data());
          tree->Add(strtmp);
          MELAout << "Adding " << strtmp << endl;
        }
      }
      DATA_SIMPLE(unsigned int, nak4jets_preselected);
      DATA_SIMPLE(float, pfmet);
      DATA_SIMPLE(float, photon_pt);
      DATA_SIMPLE(float, uParallel);
      DATA_SIMPLE(float, uPerp);
      DATA_SIMPLE(float, jetHT);
      DATA_SIMPLE(float, MET_Parallel);
      DATA_SIMPLE(float, MET_Perp);

      tout_Data.Branch("nak4jets_preselected", &nak4jets_preselected);
      tout_Data.Branch("MET_Parallel", &MET_Parallel);
      tout_Data.Branch("MET_Perp", &MET_Perp);

      int nEntries = tree->GetEntries();
      size_t nDataValid=0;
      for (int ev=0; ev<nEntries; ev++){
        tree->GetEntry(ev);
        HelperFunctions::progressbar(ev, nEntries);

        for (size_t iv=0; iv<varlist.size(); iv++){
          auto& var = varlist.at(iv);
          if (var.name == "nak4jets_preselected") var.setVal(nak4jets_preselected, 1);
          else{
            if (nak4jets_preselected==0) continue;
            if (var.name == "photon_pt") var.setVal(photon_pt, 1);
            else{
              if (photon_pt<130.f) continue;
              if (var.name == "pfmet") var.setVal(pfmet, 1);
              else if (var.name == "MET_Parallel") var.setVal(MET_Parallel, 1);
              else if (var.name == "MET_Perp") var.setVal(MET_Perp, 1);
              else if (var.name == "uParallel") var.setVal(uParallel, 1);
              else if (var.name == "uPerp") var.setVal(uPerp, 1);
              else if (var.name == "uT") var.setVal(sqrt(pow(uPerp, 2)+pow(uParallel, 2)), 1);
              else if (var.name == "jetHT") var.setVal(jetHT, 1);
              else if (var.name == "totParallel") var.setVal(photon_pt+uParallel+MET_Parallel, 1);
              else if (var.name == "totPerp") var.setVal(uPerp+MET_Perp, 1);
            }
          }

          hdatalist.at(iv).Fill(var.val, var.weights.at(0));
          var.reset();
        }

        if (photon_pt>=130.f){
          tout_Data.Fill();
          if (!(nak4jets_preselected<=0 || MET_Perp>300. || MET_Perp<-300.)) nDataValid++;
        }
      }
      MELAout << "Valid data points: " << nDataValid << " / " << nEntries << endl;
      delete tree;
    }

    // Record MC histograms
    // Scaled with 2D histograms
    vector<vector<vector<TH1F>>> hMClist;

    TString strinputscale = strinputscalecore;
    HelperFunctions::replaceString<TString, const TString>(strinputscale, "[OUTFILECORE]", sample_data);
    TFile* finput_scale = TFile::Open(strinputscale, "read");
    std::vector<TH3F*> hscale;
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_Nominal"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_JECup"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_JECdn"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_JERup"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_JERdn"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_PUup"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_PUdn"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_Data"));

    TString strinputresscale = strinputresscalecore;
    HelperFunctions::replaceString<TString, const TString>(strinputresscale, "[OUTFILECORE]", sample_data);
    TFile* finput_resscale = TFile::Open(strinputresscale, "read");
    std::vector<TH3F*> hresscale;
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_Nominal_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_JECup_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_JECdn_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_JERup_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_JERdn_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_PUup_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_PUdn_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_Data_Conditional"));

    foutput->cd();
    getMCHistograms_1D(strinputcore, sampleList_MC, varlist, hMClist, &hscale, &hresscale, &tout_MC);
    finput_resscale->Close();
    finput_scale->Close();

    for (size_t iv=0; iv<varlist.size(); iv++){
      auto& var = varlist.at(iv);

      auto const& hdata=hdatalist.at(iv);
      foutput->WriteTObject(&hdata);

      const size_t nSysts = 7;
      for (size_t isyst=0; isyst<nSysts; isyst++){
        TString strappend;
        if (isyst==0) strappend="_Nominal";
        else if (isyst==1) strappend="_JECup";
        else if (isyst==2) strappend="_JECdn";
        else if (isyst==3) strappend="_JERup";
        else if (isyst==4) strappend="_JERdn";
        else if (isyst==5) strappend="_PUup";
        else if (isyst==6) strappend="_PUdn";
        TH1F* htotMC_syst = (TH1F*) hdata.Clone(var.name+"_totMC"+strappend);
        htotMC_syst->Reset("ICESM");
        htotMC_syst->Sumw2();
        for (auto const& vs:hMClist) htotMC_syst->Add(&(vs.at(isyst).at(iv)));
        float dataYield = hdata.Integral(1, hdata.GetNbinsX());
        float MCyield = htotMC_syst->Integral(1, htotMC_syst->GetNbinsX());
        float data_MC_scale = dataYield/MCyield;
        MELAout << "Data/MC = " << dataYield << " / " << MCyield << " = " << data_MC_scale << endl;
        //htotMC_syst->Scale(data_MC_scale);
        foutput->WriteTObject(htotMC_syst);
        delete htotMC_syst;
      }
    }

    foutput->WriteTObject(&tout_Data);
    foutput->WriteTObject(&tout_MC);

    foutput->Close();
  }
}

void getFitCovarianceMatrix(RooFitResult const* fitResult, RooArgList const& ordered_args, TMatrixDSym& res){
  if (!fitResult) return;
  if (!(fitResult->status()==0 || fitResult->status()==4)) return;

  const int nFinalDimCovMat = ordered_args.getSize();

  const RooArgList pars = fitResult->floatParsFinal();
  const int nFinalPars = pars.getSize();
  TMatrixDSym mat_tmp = fitResult->covarianceMatrix();
  const int nDimCovMat = mat_tmp.GetNcols();
  if (nFinalDimCovMat<nDimCovMat) MELAout << "getFitCovarianceMatrix: Not all fit parameters are included in the ordered_args argument!" << endl;
  if (nFinalPars!=nDimCovMat){ MELAout << "getFitCovarianceMatrix: nFinalPars!=nDimCovMat! No matrix can be returned" << endl; return; }

  res.ResizeTo(nFinalDimCovMat, nFinalDimCovMat); for (int ix=0; ix<nFinalDimCovMat; ix++){ for (int iy=0; iy<nFinalDimCovMat; iy++) res[ix][iy] = 0; }
  std::vector<int> order(nFinalDimCovMat, -1);
  for (int ip=0; ip<nFinalDimCovMat; ip++){
    RooAbsArg const* target_par = ordered_args.at(ip);
    for (int jp=0; jp<nFinalPars; jp++){
      RooAbsArg const* test_par = pars.at(jp);
      if (TString(test_par->GetName())==target_par->GetName()){ order.at(ip)=jp; break; }
    }
  }

  for (int ix=0; ix<nFinalDimCovMat; ix++){
    for (int iy=0; iy<nFinalDimCovMat; iy++){
      int const& ip = order.at(ix);
      int const& jp = order.at(iy);
      if (ip<0 || jp<0) continue;
      res[ix][iy] = mat_tmp[ip][jp];
    }
  }
}

void fitFinalGammaTrees(){
  std::vector<TString> sampleList_Data={
    "Run2017B-31Mar2018-v1",
    "Run2017C-31Mar2018-v1",
    "Run2017D-31Mar2018-v1",
    "Run2017E-31Mar2018-v1",
    "Run2017F-31Mar2018-v1",
    "Run2017F-09May2018-v1"
  };

  SampleHelpers::theDataYear=2017;
  SampleHelpers::theDataPeriod="2017";
  SampleHelpers::theDataVersion=SampleHelpers::kCMSSW_9_4_X;

  TString const strinputcore = Form("output/GammaTrees/%i/plots/compare_scaled_[OUTFILECORE]_MC.root", SampleHelpers::theDataYear);
  TString const stroutputtxtcore = Form("output/GammaTrees/%i/plots/fitparameters_[OUTFILECORE]_MC.txt", SampleHelpers::theDataYear);

  RooRealVar xvar("MET_Perp", "", 0, -300, 300); xvar.setBins(150);
  RooRealVar wgtvar("weight", "", 1, -10, 10); wgtvar.removeMin(); wgtvar.removeMax();

  RooConstVar g1_mean("g1_mean", "", 0);
  RooRealVar g1_sigma("g1_sigma", "", 15, 10, 20);
  RooGaussian g1_pdf("g1_pdf", "", xvar, g1_mean, g1_sigma);
  RooRealVar g1_frac("g1_frac", "", 0, 0, 1);

  RooConstVar g2_mean("g2_mean", "", 0);
  RooRealVar g2_sigma("g2_sigma", "", 30, 20, 40);
  RooGaussian g2_pdf("g2_pdf", "", xvar, g2_mean, g2_sigma);
  RooRealVar g2_frac("g2_frac", "", 1, 0, 1);

  RooConstVar g3_mean("g3_mean", "", 0);
  RooRealVar g3_sigma("g3_sigma", "", 50, 40, 100);
  RooGaussian g3_pdf("g3_pdf", "", xvar, g3_mean, g3_sigma);

  RooAddPdf pdf("pdf", "", RooArgList(g1_pdf, g2_pdf, g3_pdf), RooArgList(g1_frac, g2_frac), true);

  for (auto const& sample_data:sampleList_Data){
    TString strinput = strinputcore;
    HelperFunctions::replaceString<TString, const TString>(strinput, "[OUTFILECORE]", sample_data);
    MELAout << "Opening input file " << strinput << endl;
    TFile* finput = TFile::Open(strinput, "read");
    std::vector<TTree*> trees ={
      (TTree*) finput->Get("FinalTree_Data"),
      (TTree*) finput->Get("FinalTree_MC")
    };

    TString stroutputtxt = stroutputtxtcore;
    HelperFunctions::replaceString<TString, const TString>(stroutputtxt, "[OUTFILECORE]", sample_data);
    // Open only when writing final parameters, but open here first to rewrite it.
    MELAout.open(stroutputtxt.Data(), std::ios_base::out);
    MELAout.close();

    for (size_t it=0; it<trees.size(); it++){
      TTree*& tree = trees.at(it);

      float totwgt=1;
      unsigned int nak4jets_preselected;
      float MET_Perp;
      float totwgt_JECup=1;
      unsigned int nak4jets_preselected_JECup;
      float MET_Perp_JECup;
      float totwgt_JECdn=1;
      unsigned int nak4jets_preselected_JECdn;
      float MET_Perp_JECdn;
      float totwgt_JERup=1;
      unsigned int nak4jets_preselected_JERup;
      float totwgt_JERdn=1;
      unsigned int nak4jets_preselected_JERdn;
      float totwgt_PUup=1;
      float totwgt_PUdn=1;
      tree->SetBranchAddress("nak4jets_preselected", &nak4jets_preselected); // Need to select Njets>0 from this tree
      tree->SetBranchAddress("MET_Perp", &MET_Perp); // = xvar
      if (it==1){
        tree->SetBranchAddress("totwgt", &totwgt); // = wgtvar

        tree->SetBranchAddress("totwgt_JECup", &totwgt_JECup); // = wgtvar
        tree->SetBranchAddress("nak4jets_preselected_JECup", &nak4jets_preselected_JECup); // Need to select Njets>0 from this tree
        tree->SetBranchAddress("MET_Perp_JECup", &MET_Perp_JECup); // = xvar
        tree->SetBranchAddress("totwgt_JECdn", &totwgt_JECdn); // = wgtvar
        tree->SetBranchAddress("nak4jets_preselected_JECdn", &nak4jets_preselected_JECdn); // Need to select Njets>0 from this tree
        tree->SetBranchAddress("MET_Perp_JECdn", &MET_Perp_JECdn); // = xvar
        tree->SetBranchAddress("totwgt_JERup", &totwgt_JERup); // = wgtvar
        tree->SetBranchAddress("nak4jets_preselected_JERup", &nak4jets_preselected_JERup); // Need to select Njets>0 from this tree
        tree->SetBranchAddress("totwgt_JERdn", &totwgt_JERdn); // = wgtvar
        tree->SetBranchAddress("nak4jets_preselected_JERdn", &nak4jets_preselected_JERdn); // Need to select Njets>0 from this tree
        tree->SetBranchAddress("totwgt_PUup", &totwgt_PUup); // = wgtvar
        tree->SetBranchAddress("totwgt_PUdn", &totwgt_PUdn); // = wgtvar
      }
      const int nSysts=(it==0 ? 1 : 7);
      for (int isyst=0; isyst<nSysts; isyst++){
        RooArgSet treevars(xvar, wgtvar);
        RooDataSet fit_data("fit_data", "", treevars, WeightVar(wgtvar));

        const int nEntries = tree->GetEntries();
        size_t nValidEntries = 0;
        size_t nValidEntries_wgtPos = 0;
        float sumWgts=0;
        float sumWgts_Pos=0;
        for (int ev=0; ev<nEntries; ev++){
          tree->GetEntry(ev);

          unsigned int* nak4jets_preselected_ptr = &nak4jets_preselected;
          float* MET_Perp_ptr = &MET_Perp;
          float* totwgt_ptr = &totwgt;
          switch (isyst){
          case 1:
            nak4jets_preselected_ptr = &nak4jets_preselected_JECup;
            MET_Perp_ptr = &MET_Perp_JECup;
            totwgt_ptr = &totwgt_JECup;
            break;
          case 2:
            nak4jets_preselected_ptr = &nak4jets_preselected_JECdn;
            MET_Perp_ptr = &MET_Perp_JECdn;
            totwgt_ptr = &totwgt_JECdn;
            break;
          case 3:
            nak4jets_preselected_ptr = &nak4jets_preselected_JERup;
            totwgt_ptr = &totwgt_JERup;
            break;
          case 4:
            nak4jets_preselected_ptr = &nak4jets_preselected_JERdn;
            totwgt_ptr = &totwgt_JERdn;
            break;
          case 5:
            totwgt_ptr = &totwgt_PUup;
            break;
          case 6:
            totwgt_ptr = &totwgt_PUdn;
            break;
          }

          if (*nak4jets_preselected_ptr>0){
            sumWgts += *totwgt_ptr;
            nValidEntries++;
            if (*totwgt_ptr>0.f){
              sumWgts_Pos += *totwgt_ptr;
              nValidEntries_wgtPos++;
            }
          }
          if (*nak4jets_preselected_ptr<=0 || *MET_Perp_ptr>xvar.getMax() || *MET_Perp_ptr<xvar.getMin()) continue;

          xvar.setVal(*MET_Perp_ptr);
          wgtvar.setVal(*totwgt_ptr);
          fit_data.add(treevars, *totwgt_ptr);
        }
        MELAout << "Number of all valid entries in the data set: " << nValidEntries << " / " << nEntries << " (sum of weights: " << sumWgts << ")" << endl;
        MELAout << "\t- Number of valid entries in the data set with positive weight: " << nValidEntries_wgtPos << " / " << nEntries << " (sum of weights: " << sumWgts_Pos << ", fraction: " << sumWgts_Pos/sumWgts << ")" << endl;
        fit_data.Print("v");

        { // Do the fit
          g1_sigma.setVal(15);
          g2_sigma.setVal(30);
          g3_sigma.setVal(50);
          if (it==0){
            g1_frac.setConstant(false);
            g1_frac.setVal(0);

            g2_frac.setConstant(false);
            g2_frac.setVal(1);
          }
          else{
            g1_frac.setConstant(true);
            g2_frac.setConstant(true);
          }

          RooLinkedList cmdList;
          RooCmdArg saveArg = RooFit::Save(true); cmdList.Add((TObject*) &saveArg);
          //RooCmdArg splitRangeArg = RooFit::SplitRange(true); cmdList.Add((TObject*) &splitRangeArg);
          RooCmdArg sumw2Arg = RooFit::SumW2Error(true);
          if (it==1) cmdList.Add((TObject*) &sumw2Arg);
          RooCmdArg hesseArg = RooFit::Hesse(true);// cmdList.Add((TObject*) &hesseArg);
          RooCmdArg initialhesseArg = RooFit::InitialHesse(true);// cmdList.Add((TObject*) &initialhesseArg);
          RooCmdArg minosArg = RooFit::Minos(true);// cmdList.Add((TObject*) &minosArg);
          //RooCmdArg minimizerArg = RooFit::Minimizer("Minuit", "migrad"); cmdList.Add((TObject*)&minimizerArg);
          RooCmdArg minimizerStrategyArg = RooFit::Strategy((it==0 ? 2 : 1)); cmdList.Add((TObject*) &minimizerStrategyArg);
          RooCmdArg cpuArg = RooFit::NumCPU(4, 0); cmdList.Add((TObject*) &cpuArg);
          // Misc. options
          RooCmdArg timerArg = RooFit::Timer(true); cmdList.Add((TObject*) &timerArg);
          //RooCmdArg printlevelArg = RooFit::PrintLevel(3); cmdList.Add((TObject*) &printlevelArg);
          RooCmdArg printlevelArg = RooFit::PrintLevel(-1); cmdList.Add((TObject*) &printlevelArg);
          RooCmdArg printerrorsArg = RooFit::PrintEvalErrors(-1); cmdList.Add((TObject*) &printerrorsArg);

          RooFitResult* fitResult=nullptr;
          int fitStatus=-1;
          unsigned int itry=0;
          constexpr unsigned int ntries=10;
          bool doImprove=false;
          constexpr bool applyImprovement=false;
          while (fitStatus!=0 || doImprove){
            if (applyImprovement && doImprove){
              MELAout << "Improving the fit result with a re-trial." << endl;
              cmdList.Add((TObject*) &hesseArg);
              cmdList.Add((TObject*) &initialhesseArg);
              cmdList.Add((TObject*) &minosArg);
            }
            delete fitResult;
            fitResult = pdf.fitTo(fit_data, cmdList);
            fitStatus = fitResult->status();
            cout << "****************************" << endl;
            MELAout << "Fitted parameters:\n";
            MELAout << "\t- Status: " << fitStatus << endl;
            MELAout << "\t- Mean 1: " << g1_mean.getVal()/* << " +- " << g1_mean.getError()*/ << endl;
            MELAout << "\t- Sigma 1: " << g1_sigma.getVal() << " +- " << g1_sigma.getError() << endl;
            MELAout << "\t- Frac 1: " << g1_frac.getVal() << " +- " << g1_frac.getError() << endl;
            MELAout << "\t- Mean 2: " << g2_mean.getVal()/* << " +- " << g2_mean.getError()*/ << endl;
            MELAout << "\t- Sigma 2: " << g2_sigma.getVal() << " +- " << g2_sigma.getError() << endl;
            MELAout << "\t- Frac 2: " << g2_frac.getVal() << " +- " << g2_frac.getError() << endl;
            MELAout << "\t- Mean 3: " << g3_mean.getVal()/* << " +- " << g3_mean.getError()*/ << endl;
            MELAout << "\t- Sigma 3: " << g3_sigma.getVal() << " +- " << g3_sigma.getError() << endl;
            cout << "****************************" << endl;
            if (applyImprovement && !doImprove && fitStatus==0) doImprove=true;
            else{
              if (!doImprove){
                itry++;
                if (itry==ntries) break;
              }
              else doImprove=false;
            }
          }
          if (fitStatus==0 || fitStatus==4){
            TString systname;
            TString systlabel;
            switch (isyst){
            case 1:
              systname="JECup";
              systlabel="JEC up";
              break;
            case 2:
              systname="JECdn";
              systlabel="JEC down";
              break;
            case 3:
              systname="JERup";
              systlabel="JER up";
              break;
            case 4:
              systname="JERdn";
              systlabel="JER down";
              break;
            case 5:
              systname="PUup";
              systlabel="PU up";
              break;
            case 6:
              systname="PUdn";
              systlabel="PU down";
              break;
            default:
              systname="Nominal";
              systlabel="nominal";
              break;
            }

            TMatrixDSym covMat;
            getFitCovarianceMatrix(fitResult, RooArgList(g1_sigma, g1_frac, g2_sigma, g2_frac, g3_sigma), covMat);
            MELAout << "****************************" << endl;
            MELAout << "Final fit properties for systematic " << systname << " in " << (it==0 ? "data" : "MC") << ":" << endl;
            fitResult->Print("v");
            MELAout.open(stroutputtxt.Data(), std::ios_base::app);
            MELAout << "****************************" << endl;
            MELAout << "Final fitted parameters for systematic " << systname << " in " << (it==0 ? "data" : "MC") << ":" << endl;
            if (it==0) MELAout << "\t- Nevents: " << nValidEntries << endl;
            MELAout << "\t- Mean 1: " << g1_mean.getVal()/* << " +- " << g1_mean.getError()*/ << endl;
            MELAout << "\t- Sigma 1: " << g1_sigma.getVal() << " +- " << g1_sigma.getError() << endl;
            MELAout << "\t- Frac 1: " << g1_frac.getVal() << " +- " << g1_frac.getError() << endl;
            MELAout << "\t- Mean 2: " << g2_mean.getVal()/* << " +- " << g2_mean.getError()*/ << endl;
            MELAout << "\t- Sigma 2: " << g2_sigma.getVal() << " +- " << g2_sigma.getError() << endl;
            MELAout << "\t- Frac 2: " << g2_frac.getVal() << " +- " << g2_frac.getError() << endl;
            MELAout << "\t- Mean 3: " << g3_mean.getVal()/* << " +- " << g3_mean.getError()*/ << endl;
            MELAout << "\t- Sigma 3: " << g3_sigma.getVal() << " +- " << g3_sigma.getError() << endl;
            //MELAout << "Covariance matrix:" << endl;
            //for (int ix=0; ix<covMat.GetNrows(); ix++){
            //  for (int iy=0; iy<covMat.GetNcols(); iy++){
            //    MELAout << covMat[ix][iy] << " ";
            //  }
            //  MELAout << endl;
            //}
            MELAout << "****************************" << endl;
            MELAout.close();

            RooPlot fit_plot(xvar, xvar.getMin(), xvar.getMax(), 150);

            fit_data.plotOn(&fit_plot, LineColor(kBlack), MarkerColor(kBlack), MarkerStyle(30), LineWidth(2), Name("Data"), XErrorSize(0)/*, Rescale(rescale_factor)*/);
            pdf.plotOn(&fit_plot, LineColor(kRed), LineWidth(2), Name("FitPdf")/*, Normalization(rescale_factor, RooAbsPdf::Relative)*/);

            fit_plot.SetTitle("");
            fit_plot.SetXTitle("E^{miss}_{T,perp} (GeV)");
            fit_plot.SetYTitle("Events");
            fit_plot.SetNdivisions(505, "X");
            fit_plot.SetLabelFont(42, "X");
            fit_plot.SetLabelOffset(0.007, "X");
            fit_plot.SetLabelSize(0.04, "X");
            fit_plot.SetTitleSize(0.06, "X");
            fit_plot.SetTitleOffset(0.9, "X");
            fit_plot.SetTitleFont(42, "X");
            fit_plot.SetNdivisions(505, "Y");
            fit_plot.SetLabelFont(42, "Y");
            fit_plot.SetLabelOffset(0.007, "Y");
            fit_plot.SetLabelSize(0.04, "Y");
            fit_plot.SetTitleSize(0.06, "Y");
            fit_plot.SetTitleOffset(1.2, "Y");
            fit_plot.SetTitleFont(42, "Y");

            TString canvasname = Form("fit_%s_%s", sample_data.Data(), (it==0 ? "Data" : "MC"));
            if (it==1) canvasname += "_" + systname;
            TCanvas can(canvasname, "", 8, 30, 800, 800);
            gStyle->SetOptStat(0);
            can.SetFillColor(0);
            can.SetBorderMode(0);
            can.SetBorderSize(2);
            can.SetTickx(1);
            can.SetTicky(1);
            can.SetLeftMargin(0.17);
            can.SetRightMargin(0.05);
            can.SetTopMargin(0.07);
            can.SetBottomMargin(0.13);
            can.SetFrameFillStyle(0);
            can.SetFrameBorderMode(0);
            can.SetFrameFillStyle(0);
            can.SetFrameBorderMode(0);

            TLegend legend(0.20, 0.90-0.15, 0.50, 0.90);
            legend.SetBorderSize(0);
            legend.SetTextFont(42);
            legend.SetTextSize(0.03);
            legend.SetLineColor(1);
            legend.SetLineStyle(1);
            legend.SetLineWidth(1);
            legend.SetFillColor(0);
            legend.SetFillStyle(0);

            TPaveText pavetext(0.15, 0.93, 0.85, 1, "brNDC");
            pavetext.SetBorderSize(0);
            pavetext.SetFillStyle(0);
            pavetext.SetTextAlign(12);
            pavetext.SetTextFont(42);
            pavetext.SetTextSize(0.045);
            TText* text = pavetext.AddText(0.025, 0.45, "#font[61]{CMS}");
            text->SetTextSize(0.044);
            if (it==0){
              text = pavetext.AddText(0.165, 0.42, "#font[52]{Preliminary}");
              text->SetTextSize(0.0315);
            }
            else{
              text = pavetext.AddText(0.165, 0.42, "#font[52]{Simulation}");
              text->SetTextSize(0.0315);
            }
            TString cErgTev = Form("#font[42]{%i TeV (%s)}", theSqrts, theDataPeriod.Data());
            text = pavetext.AddText(0.87, 0.45, cErgTev);
            text->SetTextSize(0.0315);

            TString strDataTitle=(it==0 ? "Data" : "Simulation");
            TString strPdfTitle="Fit";
            fit_plot.Draw();
            //reference_hists.at(it)->Draw("histsame");
            //xcheck_hists.at(it).Draw("histsame");
            TString datalabel = strDataTitle; if (it==1) datalabel = datalabel + " (" + systlabel + ")";
            legend.AddEntry("Data", datalabel, "lp");
            legend.AddEntry("FitPdf", strPdfTitle, "l");
            //legend.AddEntry(reference_hists.at(it), "Reference hist.", "l");
            legend.Draw("same");
            pavetext.Draw();
            can.RedrawAxis();
            can.Modified();
            can.Update();
            can.SaveAs(Form("%s.pdf", can.GetName()));
            can.SaveAs(Form("%s.root", can.GetName()));
            can.Close();
          }
          delete fitResult;
        }


      }
    }
    finput->Close();
  }
}

void checkCorrectedMETDistributions(){
  std::vector<TString> sampleList_Data={
    "Run2017B-31Mar2018-v1",
    "Run2017C-31Mar2018-v1",
    "Run2017D-31Mar2018-v1",
    "Run2017E-31Mar2018-v1",
    "Run2017F-31Mar2018-v1",
    "Run2017F-09May2018-v1"
  };
  std::vector<TString> sampleList_MC={
    "GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8",
    "GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",

    "QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8",
    "QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8",

    "TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8",
    "TGJets_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8",
    //"TTGamma_Dilept_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromT_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    //"TTGamma_SingleLeptFromTbar_TuneCP5_PSweights_13TeV_madgraph_pythia8",
    "TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8",

    "WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WZG_TuneCP5_13TeV-amcatnlo-pythia8",

    "ZJetsToNuNu_HT-100To200_13TeV-madgraph",
    "ZJetsToNuNu_HT-200To400_13TeV-madgraph",
    "ZJetsToNuNu_HT-400To600_13TeV-madgraph",
    "ZJetsToNuNu_HT-600To800_13TeV-madgraph",
    "ZJetsToNuNu_HT-800To1200_13TeV-madgraph",
    "ZJetsToNuNu_HT-1200To2500_13TeV-madgraph",
    "ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph"
  };

  SampleHelpers::theDataYear=2017;
  SampleHelpers::theDataVersion=SampleHelpers::kCMSSW_9_4_X;

  METCorrectionHandler metCorrector;

  std::vector<Variable> varlist={
    Variable("pfmet", "E_{T}^{miss} (GeV)", 100, 0, 200),
    Variable("MET_Parallel", "E_{T,//}^{miss} (GeV)", 150, -300, 300),
    Variable("MET_Perp", "E_{T,perp}^{miss} (GeV)", 100, -200, 200),

    Variable("pfmet_corrected", "E_{T}^{miss} (GeV)", 100, 0, 200),
    Variable("MET_Parallel_corrected", "E_{T,//}^{miss} (GeV)", 150, -300, 300),
    Variable("MET_Perp_corrected", "E_{T,perp}^{miss} (GeV)", 100, -200, 200),
  };

  TString const strinputcore = Form("output/GammaTrees/%i/[OUTFILECORE].root", SampleHelpers::theDataYear);
  TString const strinputdatahistcore = Form("output/GammaTrees/%i/plots/compare_scaled_[OUTFILECORE]_MC.root", SampleHelpers::theDataYear);
  TString const strinputscalecore = Form("output/GammaTrees/%i/plots/scale_[OUTFILECORE]_MC.root", SampleHelpers::theDataYear);
  TString const strinputresscalecore = Form("output/GammaTrees/%i/plots/scale_residual_[OUTFILECORE]_MC.root", SampleHelpers::theDataYear);
  TString const stroutputcore = Form("output/GammaTrees/%i/plots/checkCorrectedMET_[OUTFILECORE]_MC.root", SampleHelpers::theDataYear);

  for (auto const& sample_data:sampleList_Data){
    for (TString const& period:SampleHelpers::getValidDataPeriods()){
      if (sample_data.Contains(period)) SampleHelpers::theDataPeriod = period;
    }
    if (sample_data.Contains("09May2018")) SampleHelpers::theDataPeriod += "-09May2018";
    MELAout << "Data period is " << SampleHelpers::theDataPeriod << endl;
    metCorrector.setVerbosity(TVar::DEBUG);
    metCorrector.setup();
    metCorrector.setVerbosity(TVar::ERROR);

    TString stroutput = stroutputcore;
    HelperFunctions::replaceString<TString, const TString>(stroutput, "[OUTFILECORE]", sample_data);
    MELAout << "Creating output file " << stroutput << endl;
    TFile* foutput = TFile::Open(stroutput, "recreate");

    TString strinputdatahist = strinputdatahistcore;
    HelperFunctions::replaceString<TString, const TString>(strinputdatahist, "[OUTFILECORE]", sample_data);
    TFile* finput_datahist = TFile::Open(strinputdatahist, "read");
    vector<TH1F*> hdatalist; hdatalist.reserve(varlist.size()/2); // Corrected or uncorrected in data are just the same
    for (size_t iv=0; iv<varlist.size()/2; iv++) hdatalist.push_back((TH1F*) finput_datahist->Get(Form("%s_Data", varlist.at(iv).name.Data())));

    TString strinputscale = strinputscalecore;
    HelperFunctions::replaceString<TString, const TString>(strinputscale, "[OUTFILECORE]", sample_data);
    TFile* finput_scale = TFile::Open(strinputscale, "read");
    std::vector<TH3F*> hscale;
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_Nominal"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_JECup"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_JECdn"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_JERup"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_JERdn"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_PUup"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_totMC_PUdn"));
    hscale.push_back((TH3F*) finput_scale->Get("nak4jets_preselected_jetHT_photon_pt_Data"));

    TString strinputresscale = strinputresscalecore;
    HelperFunctions::replaceString<TString, const TString>(strinputresscale, "[OUTFILECORE]", sample_data);
    TFile* finput_resscale = TFile::Open(strinputresscale, "read");
    std::vector<TH3F*> hresscale;
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_Nominal_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_JECup_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_JECdn_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_JERup_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_JERdn_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_PUup_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_totMC_PUdn_Conditional"));
    hresscale.push_back((TH3F*) finput_resscale->Get("nak4jets_preselected_abs_uPerp_uParallel_Data_Conditional"));

    foutput->cd();
    vector<vector<vector<TH1F>>> hMClist;
    getMCHistograms_1D(strinputcore, sampleList_MC, varlist, hMClist, &hscale, &hresscale, nullptr, &metCorrector);

    for (auto*& hdata:hdatalist) foutput->WriteTObject(hdata);
    if (!hMClist.empty()){
      for (size_t iv=0; iv<varlist.size(); iv++){
        auto& var = varlist.at(iv);

        TH1F*& hdata = hdatalist.at(iv % (varlist.size()/2));

        const size_t nSysts = hMClist.front().size();
        for (size_t isyst=0; isyst<nSysts; isyst++){
          TString strappend;
          if (isyst==0) strappend="_Nominal";
          else if (isyst==1) strappend="_JECup";
          else if (isyst==2) strappend="_JECdn";
          else if (isyst==3) strappend="_JERup";
          else if (isyst==4) strappend="_JERdn";
          else if (isyst==5) strappend="_PUup";
          else if (isyst==6) strappend="_PUdn";
          TH1F* htotMC_syst = (TH1F*) hdata->Clone(var.name+"_totMC"+strappend);
          htotMC_syst->Reset("ICESM");
          htotMC_syst->Sumw2();
          for (auto const& vs:hMClist) htotMC_syst->Add(&(vs.at(isyst).at(iv)));
          float dataYield = hdata->Integral(1, hdata->GetNbinsX());
          float MCyield = htotMC_syst->Integral(1, htotMC_syst->GetNbinsX());
          float data_MC_scale = dataYield/MCyield;
          MELAout << "Data/MC = " << dataYield << " / " << MCyield << " = " << data_MC_scale << endl;
          //htotMC_syst->Scale(data_MC_scale);
          foutput->WriteTObject(htotMC_syst);
          delete htotMC_syst;
        }
      }
    }

    finput_resscale->Close();
    finput_scale->Close();
    finput_datahist->Close();
    foutput->Close();
  }

}
