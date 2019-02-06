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

void plotJetResolution(){
  gStyle->SetOptStat(0);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  const TString strindir="output/MTD/2018/";
  TString stroutputcore=strindir+"/plots/JetResolution/";
  TString stroutput=stroutputcore + "QCD_comparison.root";
  std::vector<TString> sampleList{
    "QCD_Pt-15to20_EMEnriched_TuneCP5_13TeV_pythia8.root",
    "QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8.root",
    "QCD_Pt-170to300_EMEnriched_TuneCP5_13TeV_pythia8.root",
    "QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8.root"
  };
  size_t const nsamples = sampleList.size();
  std::vector<TString> histNames{
    "QCD_Pt-15to20_EMEnriched",
    "QCD_Pt-50to80_EMEnriched",
    "QCD_Pt-170to300_EMEnriched",
    "QCD_Pt-600to800_MuEnrichedPt5"
  };
  std::vector<TString> histTitles{
    "p^{true}_{T}: [10, 30) GeV",
    "p^{true}_{T}: [30, 100) GeV",
    "p^{true}_{T}: [100, 300) GeV",
    "p^{true}_{T} >= 300 GeV"
  };
  std::vector<int> histColor{
    (int) kBlue,
    (int) kGreen+2,
    (int) kOrange-3,
    (int) kRed
  };
  size_t const nhists = histColor.size();

  std::vector<Variable> plotvars;
  plotvars.emplace_back("mu", "<E_{jet}^{reco}/E_{jet}^{gen}-1>", 50, -5, 5);
  plotvars.emplace_back("sigma", "#sigma(E_{jet}^{reco}/E_{jet}^{gen}-1)", 100, 0, 10);
  size_t const nvars = plotvars.size();

  TFile* foutput = TFile::Open(stroutput, "recreate");
  std::vector<TProfile*> histColl(sampleList.size(), nullptr);
  std::vector<std::vector<TH1F*>> histVarColl(nvars, std::vector<TH1F*>(sampleList.size(), nullptr));
  for (size_t is=0; is<nhists; is++){
    foutput->cd();

    TProfile*& h = histColl.at(is);
    h = new TProfile(histNames.at(is), "", 50, -5, 5);
    h->SetLineColor(histColor.at(is));
    h->SetLineWidth(2);
    h->SetMarkerColor(histColor.at(is));
    h->GetXaxis()->SetTitle("#eta");
    h->GetYaxis()->SetTitle(plotvars.at(0).title);
    h->Sumw2();

    for (unsigned int ivar=0; ivar<nvars; ivar++){
      TH1F*& hh = histVarColl.at(ivar).at(is);
      hh = new TH1F(histNames.at(is)+"_"+plotvars.at(ivar).name, "", 50, -5, 5);
      hh->SetLineColor(histColor.at(is));
      hh->SetLineWidth(2);
      hh->SetMarkerColor(histColor.at(is));
      hh->GetXaxis()->SetTitle("#eta");
      hh->GetYaxis()->SetTitle(plotvars.at(ivar).title);
      hh->Sumw2();
    }
  }
  for (size_t is=0; is<nsamples; is++){
    TString const& fname = sampleList.at(is);
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

    // AK4Jets
    DATA_CLASS(std::vector<float>, ak4jets_pt)
    DATA_CLASS(std::vector<float>, ak4jets_eta)
    DATA_CLASS(std::vector<float>, ak4jets_phi)
    DATA_CLASS(std::vector<float>, ak4jets_mass)
    DATA_CLASS(std::vector<long long>, ak4jets_selectionBits)
    DATA_CLASS(std::vector<int>, ak4jets_genjetIndex)
    DATA_CLASS(std::vector<float>, ak4jets_genjetDeltaR)
    //DATA_CLASS(std::vector<float>, ak4jets_estimatedPtResolution)

    // Gen. jets
    DATA_CLASS(std::vector<float>, genjets_pt)
    DATA_CLASS(std::vector<float>, genjets_eta)
    DATA_CLASS(std::vector<float>, genjets_phi)
    DATA_CLASS(std::vector<float>, genjets_mass)

    float sumCtrWgt[WeightVariables::nWeightTypes]={ 0 };
    int const nEntries = tree->GetEntries();
    for (int ev=0; ev<nEntries; ev++){
      tree->GetEntry(ev);
      HelperFunctions::progressbar(ev, nEntries);

      if (ev==0){ theXSec=xsec; tree->SetBranchStatus("xsec", 0); }

      for (int iwt=0; iwt<WeightVariables::nWeightTypes; iwt++) sumCtrWgt[iwt] += genweights[iwt];

      size_t const ngenjets = genjets_pt->size();
      std::vector<MELAParticle> genjets; genjets.reserve(ngenjets);
      for (size_t ip=0; ip<ngenjets; ip++){
        TLorentzVector v; v.SetPtEtaPhiM(genjets_pt->at(ip), genjets_eta->at(ip), genjets_phi->at(ip), genjets_mass->at(ip));
        genjets.emplace_back(0, v);
      }

      size_t const nak4jets = ak4jets_pt->size();
      std::vector<MELAParticle> ak4jets; ak4jets.reserve(nak4jets);
      std::vector<std::pair<unsigned int, float>> ak4jets_preselected_genjetIndexDeltaR; ak4jets_preselected_genjetIndexDeltaR.reserve(nak4jets);
      for (size_t ip=0; ip<nak4jets; ip++){
        if (!HelperFunctions::test_bit(ak4jets_selectionBits->at(ip), AK4JetSelectionHelpers::kTightID)) continue;
        if (ak4jets_genjetIndex->at(ip)<0 || ak4jets_genjetDeltaR->at(ip)>0.2) continue;
        TLorentzVector v; v.SetPtEtaPhiM(ak4jets_pt->at(ip), ak4jets_eta->at(ip), ak4jets_phi->at(ip), ak4jets_mass->at(ip));
        ak4jets.emplace_back(0, v);
        ak4jets_preselected_genjetIndexDeltaR.emplace_back(ak4jets_genjetIndex->at(ip), ak4jets_genjetDeltaR->at(ip));

        MELAParticle const& ak4jet = ak4jets.back();
        MELAParticle const& genjet = genjets.at(ak4jets_preselected_genjetIndexDeltaR.back().first);
        float genpt = genjet.pt();
        float genE = genjet.t();
        int ihist=-1;
        if (genpt>=10.f && genpt<30.f) ihist=0;
        else if (genpt>=30.f && genpt<100.f) ihist=1;
        else if (genpt>=100.f && genpt<300.f) ihist=2;
        else if (genpt>300.f) ihist=3;
        if (genE>0. && ihist>=0){
          float var = ak4jet.t()/genE-1.;
          histColl.at(ihist)->Fill(ak4jet.eta(), var, genweights[0]);
        }
      }
      //size_t nak4jets_preselected = ak4jets.size();
    }
    finput->Close();
  }

  foutput->cd();
  for(size_t ivar=0;ivar<nvars;ivar++){
    float ymin = std::numeric_limits<float>::max();
    float ymax = std::numeric_limits<float>::min();

    for (size_t is=0; is<nhists; is++){
      TProfile*& h = histColl.at(is);
      TH1F*& hh = histVarColl.at(ivar).at(is);
      for (int ix=0; ix<=h->GetNbinsX()+1; ix++){
        float bc=0, be=0;
        if (plotvars.at(ivar).name=="mu"){
          bc = h->GetBinContent(ix);
          be = h->GetBinError(ix);
        }
        else if (plotvars.at(ivar).name=="sigma"){
          bc = h->GetBinError(ix)*sqrt(h->GetBinEffectiveEntries(ix));
          be = 0;
        }
        hh->SetBinContent(ix, bc);
        hh->SetBinError(ix, be);
        if (ix<1 || ix>h->GetNbinsX()) continue;
        ymin=std::min(bc-be, ymin);
        ymax=std::max(bc+be, ymax);
      }
    }
    for (auto& h:histVarColl.at(ivar)){
      float scaleFacUp=1.7;
      h->GetYaxis()->SetRangeUser(ymin*(ymin>=0.f ? 0.8 : 1.2), ymax*scaleFacUp);
    }

    {
      // Plot the histograms
      TString canvasname = plotvars.at(ivar).name+"_comparison";
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

      size_t nLegend=nhists;

      const float legend_minX = 0.50;
      const float legend_maxY = 0.90;
      const float legend_maxX = 0.90;
      float legend_minY = legend_maxY - 0.06f*float(nLegend);

      TLegend legend(legend_minX, legend_minY, legend_maxX, legend_maxY);
      legend.SetBorderSize(0);
      legend.SetTextFont(42);
      legend.SetTextSize(0.04);
      legend.SetLineColor(0);
      legend.SetLineStyle(0);
      legend.SetLineWidth(1);
      legend.SetFillColor(0);
      legend.SetFillStyle(0);

      for (size_t is=0; is<nhists; is++){
        TString strLabel=histTitles.at(is);
        legend.AddEntry(histVarColl.at(ivar).at(is), strLabel, "l");
        if (is==0) histVarColl.at(ivar).at(is)->Draw("hist");
        else histVarColl.at(ivar).at(is)->Draw("histsame");
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
  }

  for (auto& col:histVarColl){ for (auto& h:col){ foutput->WriteTObject(h); delete h; } }
  for (auto& h:histColl){ foutput->WriteTObject(h); delete h; }
  foutput->Close();
}
