#include "common_includes.h"
#include "RooMsgService.h"
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
  for (auto lep=leps.cbegin(); lep!=leps.cend(); lep++){
    for (auto alep=aleps.cbegin(); alep!=aleps.cend(); alep++){
      if (PDGHelpers::getCoupledVertex(lep.id, alep.id)!=23) continue;
      if (
        (!Zpair.first || !Zpair.second)
        ||
        fabs((lep.p4 + alep.p4).M()-PDGHelpers::Zmass)<fabs((Zpair.first->p4 + Zpair.second->p4).M()-PDGHelpers::Zmass)
        ){ Zpair.first = &lep; Zpair.second.second = &alep; }
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

void testAllFits(TString fname){
  gStyle->SetOptStat(0);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  SampleHelpers::theDataYear=2018;
  SampleHelpers::theDataPeriod="2018";
  SampleHelpers::theDataVersion=SampleHelpers::kCMSSW_10_X;
  float const btagWP = BtagHelpers::getBtagWP(BtagHelpers::kDeepCSV_Medium);
  float const totalLumi=150.f*1000.f;

  RooConstVar one_rcv("one_rcv", "one_rcv", 1);
  RooConstVar zero_rcv("zero_rcv", "zero_rcv", 0);

  RooRelBW3ProngPdf::modelParameters Tpars;
  Tpars.mX = new RooConstVar("mT", "mT", 173.2);
  Tpars.gamX = new RooConstVar("gamT", "gamT", 1.4);
  Tpars.mV = new RooConstVar("mW", "mW", 80.4);
  Tpars.gamV = new RooConstVar("gamW", "gamW", 2.);

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


  TFile* foutput = TFile::Open(stroutput, "recreate");
  TTree tout("ZFitTest", "");
  tout.Branch("pfmet", &pfmet);
  tout.Branch("pfmetPhi", &pfmetPhi);
  BOOK_BRANCH(float, mT, tout);
  BOOK_BRANCH(unsigned int, nak4jets_preselected, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_pt, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_eta, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_phi, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_mass, tout);
  BOOK_BRANCH(std::vector<float>, ak4jets_preselected_relPtResolution, tout);
  BOOK_BRANCH(float, bestTFTop_mass, tout);
  BOOK_BRANCH(float, bestTFTop_disc, tout);
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

  BOOK_BRANCH(std::vector<unsigned int>, hadTfit_jetIndices, tout);
  BOOK_BRANCH(std::vector<float>, hadTfit_jetnuisance, tout);
  BOOK_BRANCH(int, hadTfit_status, tout);
  BOOK_BRANCH(int, hadTfit_VId, tout);
  BOOK_BRANCH(float, hadTfit_mW_prefit, tout);
  BOOK_BRANCH(float, hadTfit_mW_postfit, tout);
  BOOK_BRANCH(float, hadTfit_mtop_prefit, tout);
  BOOK_BRANCH(float, hadTfit_mtop_postfit, tout);
  BOOK_BRANCH(float, bestHadTFitMatchTFTop_mass, tout);
  BOOK_BRANCH(float, bestHadTFitMatchTFTop_disc, tout);
  BOOK_BRANCH(float, bestHadTFitMatchTFTop_dR, tout);

  BOOK_BRANCH(std::vector<unsigned int>, hadVfit_jetIndices, tout);
  BOOK_BRANCH(std::vector<float>, hadVfit_jetnuisance, tout);
  BOOK_BRANCH(int, hadVfit_status, tout);
  BOOK_BRANCH(int, hadVfit_VId, tout);
  BOOK_BRANCH(float, hadVfit_mV_prefit, tout);
  BOOK_BRANCH(float, hadVfit_mV_postfit, tout);


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
  plotvars.emplace_back("hadTfit_mW_prefit", "m_{qq} (pre-fit) (GeV)", 60, 70, 100);
  plotvars.emplace_back("hadTfit_mW_postfit", "m_{qq} (post-fit) (GeV)", 60, 70, 100);
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

    mT=-1;
    nak4jets_preselected=0;
    ak4jets_preselected_pt.clear();
    ak4jets_preselected_eta.clear();
    ak4jets_preselected_phi.clear();
    ak4jets_preselected_mass.clear();
    ak4jets_preselected_relPtResolution.clear();
    bestTFTop_mass=-1;
    bestTFTop_disc=-1;
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

    hadTfit_jetIndices.clear();
    hadTfit_jetnuisance.clear();
    hadTfit_status=-1;
    hadTfit_VId=-9000;
    hadTfit_mW_prefit=-1;
    hadTfit_mW_postfit=-1;
    hadTfit_mtop_prefit=-1;
    hadTfit_mtop_postfit=-1;
    bestHadTFitMatchTFTop_mass=-1;
    bestHadTFitMatchTFTop_disc=-1;
    bestHadTFitMatchTFTop_dR=-1;

    hadVfit_jetIndices.clear();
    hadVfit_jetnuisance.clear();
    hadVfit_status=-1;
    hadVfit_VId=-9000;
    hadVfit_mV_prefit=-1;
    hadVfit_mV_postfit=-1;


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


    size_t const nak4jets = ak4jets_pt->size();
    std::vector<MELAParticle> ak4jets; ak4jets.reserve(nak4jets);
    for (size_t ip=0; ip<nak4jets; ip++){
      if (!HelperFunctions::test_bit(ak4jets_selectionBits->at(ip), AK4JetSelectionHelpers::kTightID) || ak4jets_pt->at(ip)<AK4JetSelectionHelpers::ptThr_skim_preselection || fabs(ak4jets_eta->at(ip))>=4.7) continue;
      TLorentzVector v; v.SetPtEtaPhiM(ak4jets_pt->at(ip), ak4jets_eta->at(ip), ak4jets_phi->at(ip), ak4jets_mass->at(ip));
      ak4jets.emplace_back(0, v);

      ak4jets_preselected_pt.push_back(v.Pt());
      ak4jets_preselected_eta.push_back(v.Eta());
      ak4jets_preselected_phi.push_back(v.Phi());
      ak4jets_preselected_mass.push_back(v.M());
      ak4jets_preselected_relPtResolution.push_back(ak4jets_estimatedPtResolution->at(ip)/v.Pt());
    }
    nak4jets_preselected = ak4jets.size();

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


    int bestfitstatus=-999;
    float bestpdfval=-1;
    float bestconstraintpdfval=-1;
    float bestWmass = -1;
    float bestWmass_prefit = -1;
    float bestTmass = -1;
    float bestTmass_prefit = -1;
    std::vector<float> topfit_jetpt;
    std::vector<float> topfit_jeteta;
    std::vector<float> topfit_jetphi;
    std::vector<float> topfit_jetmass;
    std::vector<float> topfit_jetrelresolution;
    std::vector<float> topfit_jetnuisance;

    if (nak4jets_preselected>=3 && nak4jets_preselected<=4){
      std::vector<std::vector<int>> comblist;
      TNumericUtil::PermutationGenerator(nak4jets_preselected, 3, comblist, 0);

      for (auto const& comb:comblist){
        size_t const nprongs = comb.size();
        if (comb.at(1)>comb.at(2)) continue;
        RooRelBW3ProngPdf::modelMeasurables vars;
        std::vector<RooConstVar> ptraw, eta, phi, massraw, kappaLow, kappaHigh;
        ptraw.reserve(nprongs); eta.reserve(nprongs); phi.reserve(nprongs); massraw.reserve(nprongs); kappaLow.reserve(nprongs); kappaHigh.reserve(nprongs);
        std::vector<RooRealVar> theta_ptunc; theta_ptunc.reserve(nprongs);
        std::vector<AsymPow> asympow_mult; asympow_mult.reserve(nprongs);
        std::vector<RooFormulaVar> pt_rfv, mass_rfv;
        pt_rfv.reserve(nprongs); mass_rfv.reserve(nprongs);
        std::vector<RooGaussian> pdf_ptunc; pdf_ptunc.reserve(nprongs);
        for (unsigned int ijet=0; ijet<nprongs; ijet++){
          const int& ipart = comb.at(ijet);

          kappaLow.emplace_back(Form("kappaLow_jet%i", ijet+1), Form("kappaLow_jet%i", ijet+1), 1.f/(1.f+ak4jets_estimatedPtResolution->at(ipart)/ak4jets.at(ipart).pt()));
          kappaHigh.emplace_back(Form("kappaHigh_jet%i", ijet+1), Form("kappaHigh_jet%i", ijet+1), 1.f+ak4jets_estimatedPtResolution->at(ipart)/ak4jets.at(ipart).pt());
          theta_ptunc.emplace_back(Form("theta_ptunc_jet%i", ijet+1), Form("theta_ptunc_jet%i", ijet+1), 0., -5., 5.);
          asympow_mult.emplace_back(Form("asympow_mult_jet%i", ijet+1), Form("asympow_mult_jet%i", ijet+1), kappaLow.back(), kappaHigh.back(), theta_ptunc.back());
          pdf_ptunc.emplace_back(Form("pdf_ptunc_jet%i", ijet+1), Form("pdf_ptunc_jet%i", ijet+1), theta_ptunc.back(), zero_rcv, one_rcv);

          eta.emplace_back(Form("eta_jet%i", ijet+1), Form("eta_jet%i", ijet+1), ak4jets.at(ipart).eta());
          phi.emplace_back(Form("phi_jet%i", ijet+1), Form("phi_jet%i", ijet+1), ak4jets.at(ipart).phi());

          ptraw.emplace_back(Form("ptraw_jet%i", ijet+1), Form("ptraw_jet%i", ijet+1), ak4jets.at(ipart).pt());
          pt_rfv.emplace_back(Form("pt_rfv_jet%i", ijet+1), "(@0*@1)", RooArgList(ptraw.back(), asympow_mult.back()));

          massraw.emplace_back(Form("massraw_jet%i", ijet+1), Form("massraw_jet%i", ijet+1), ak4jets.at(ipart).m());
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

        RooRelBW3ProngPdf topWpdf("topWpdf", "topWpdf", Tpars, vars, RooRelBW3ProngPdf::kWany);
        RooProdPdf totalpdf("totalpdf", "totalpdf", RooArgList(topWpdf, constraintspdf));
        RooFormulaVar logL_rfv("logLikelihood_rfv", "-log(@0)", RooArgList(totalpdf));

        float topmass_prefit = topWpdf.getQXSq(); // Not important from which pdf you read it
        float Wmass_prefit = topWpdf.getQVSq(); // Not important from which pdf you read it
        /*
        bool const dofits = 
          (topmass_prefit>0.f && sqrt(topmass_prefit)/RooRelBW3ProngPdf::GeVunit>=100.f && sqrt(topmass_prefit)/RooRelBW3ProngPdf::GeVunit<=350.f)
          &&
          (Wmass_prefit>0.f && sqrt(Wmass_prefit)/RooRelBW3ProngPdf::GeVunit>=40.f && sqrt(Wmass_prefit)/RooRelBW3ProngPdf::GeVunit<=210.f)
          ;
        if (!dofits) continue;
        */
        RooFitResult* fitResult = nullptr;

        float constraintsStatus;
        std::vector<float> thetaValues; thetaValues.reserve(nprongs);
        float tmppdfval;
        float topmass_postfit;
        float Wmass_postfit;
        bool isFitOK=false;
        int fitstatus=-999;
        { // Do T+W minimization
          for (auto& var:theta_ptunc) var.setVal(0);
          RooMinuit minimizer(logL_rfv);
          minimizer.setPrintLevel(-1);
          minimizer.setNoWarn();
          minimizer.setStrategy(2);
          minimizer.migrad();
          fitResult = minimizer.save();
          if (fitResult){
            fitstatus=fitResult->status();
            isFitOK=(fitstatus==0);
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
          tmppdfval = totalpdf.getVal();
          topmass_postfit = topWpdf.getQXSq();
          Wmass_postfit = topWpdf.getQVSq();
          constraintsStatus = constraintspdf.getVal();
          for (size_t ijet=0; ijet<nprongs; ijet++) thetaValues.push_back(theta_ptunc.at(ijet).getVal());
          //MELAout << "\t\t- Masses for W and top: " << sqrt(Wmass_prefit)/100.f  << ", " << sqrt(topmass_prefit)/100.f << " -> " << sqrt(Wmass_postfit)/100.f << ", " << sqrt(topmass_postfit)/100.f << endl;
        }

        if (isFitOK && (bestpdfval<0.f || tmppdfval>bestpdfval)){
          bestpdfval = tmppdfval;
          bestconstraintpdfval = constraintsStatus;
          bestTmass = topmass_postfit;
          bestTmass_prefit = topmass_prefit;
          bestWmass = Wmass_postfit;
          bestWmass_prefit = Wmass_prefit;
          bestfitstatus = fitstatus;

          topfit_jetpt.clear();
          topfit_jeteta.clear();
          topfit_jetphi.clear();
          topfit_jetmass.clear();
          topfit_jetrelresolution.clear();
          topfit_jetnuisance = thetaValues;
          for (unsigned int ijet=0; ijet<nprongs; ijet++){
            topfit_jetpt.push_back(ptraw.at(ijet).getVal());
            topfit_jetmass.push_back(massraw.at(ijet).getVal());
            topfit_jeteta.push_back(eta.at(ijet).getVal());
            topfit_jetphi.push_back(phi.at(ijet).getVal());
            topfit_jetrelresolution.push_back(kappaHigh.at(ijet).getVal()-1.);
          }
        }
      }
    }

    if (bestpdfval>0.f){
      hadTfit_VId=6;
      hadTfit_mW_prefit = sqrt(bestWmass_prefit)/RooRelBW3ProngPdf::GeVunit;
      hadTfit_mW_postfit = sqrt(bestWmass)/RooRelBW3ProngPdf::GeVunit;
      hadTfit_mtop_prefit = sqrt(bestTmass_prefit)/RooRelBW3ProngPdf::GeVunit;
      hadTfit_mtop_postfit = sqrt(bestTmass)/RooRelBW3ProngPdf::GeVunit;
      hadTfit_jetpt = topfit_jetpt;
      hadTfit_jeteta = topfit_jeteta;
      hadTfit_jetphi = topfit_jetphi;
      hadTfit_jetmass = topfit_jetmass;
      hadTfit_jetrelresolution = topfit_jetrelresolution;
      hadTfit_jetnuisance = topfit_jetnuisance;
      hadTfit_status = bestfitstatus;
      //MELAout << "\t- Best masses for W and top: " << hadTfit_mW_prefit << ", " << hadTfit_mtop_prefit << " -> " << hadTfit_mW_postfit << ", " << hadTfit_mtop_postfit << endl;
      if (hadTfit_jetpt.size()>=3){
        TLorentzVector p_top_prefit;
        for (unsigned int i=0; i<3; i++){ TLorentzVector pjet; pjet.SetPtEtaPhiM(hadTfit_jetpt.at(i), hadTfit_jeteta.at(i), hadTfit_jetphi.at(i), hadTfit_jetmass.at(i)); p_top_prefit += pjet; }

        for (size_t ip=0; ip<ntftops; ip++){
          auto const& tftop = tftops.at(ip);
          float tmpdr = tftop.deltaR(p_top_prefit);
          if (bestHadTFitMatchTFTop_dR<0.f || tmpdr<bestHadTFitMatchTFTop_dR){
            bestHadTFitMatchTFTop_dR = tmpdr;
            bestHadTFitMatchTFTop_mass = tftop.m();
          }
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
      else if (var.name.BeginsWith("hadTfit_mW_postfit")){ if (hadTfit_mW_postfit<0.f){ continue; } var.setVal(hadTfit_mW_postfit, genweights); }
      else if (var.name.BeginsWith("hadTfit_mW_prefit")){ if (hadTfit_mW_prefit<0.f){ continue; } var.setVal(hadTfit_mW_prefit, genweights); }
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

  delete Tpars.mX;
  delete Tpars.gamX;
  delete Tpars.mV;
  delete Tpars.gamV;
}
