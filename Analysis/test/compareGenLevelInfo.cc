#include "common_includes.h"
#include <limits>
#include "TStyle.h"
#include "TLegend.h"


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

void compareGenLevelInfo(TString fname){
  gStyle->SetOptStat(0);

  bool isWZ = (fname.Contains("WZTo"));
  bool isWZZ = (fname.Contains("WZZ"));
  bool isWWW = (fname.Contains("WWW"));
  bool isDYLLJets = (fname.Contains("DYJetsToLL"));
  bool isTT2L2Nu = (fname.Contains("TTTo2L2Nu"));
  bool isTT2L2Nu_0j_special = (fname.Contains("TT_DiLept"));
  bool isTT2L2Nu_1j_special = (fname.Contains("TTPlus1Jet_DiLept"));
  bool isZZ4L = (fname.Contains("ZZTo4L"));
  if (!fname.Contains(".root")) fname += ".root";
  const TString strindir="output/test_GenWeights/2018/";
  TString stroutputcore=strindir+"/plots/gen/";
  TString stroutput=stroutputcore;
  TString sampleLabel;
  if (isWZ) sampleLabel="WZTo3LNu_comparison";
  if (isWZZ) sampleLabel="WZZ_comparison";
  if (isWWW) sampleLabel="WWW_comparison";
  if (isDYLLJets) sampleLabel="DYJetsToLL_comparison";
  if (isTT2L2Nu) sampleLabel="TTTo2L2Nu_comparison";
  if (isTT2L2Nu_0j_special) sampleLabel="TTTo2L2Nu_0j_special_comparison";
  if (isTT2L2Nu_1j_special) sampleLabel="TTTo2L2Nu_1j_special_comparison";
  if (isZZ4L) sampleLabel="ZZTo4L_comparison";
  stroutput += sampleLabel + ".root";

  TFile* finput = TFile::Open(strindir+fname, "read");
  TTree* tree = (TTree*) finput->Get("test");
  
  std::vector<float> genweights(WeightVariables::nWeightTypes, 0.f);
  std::vector<float> oneWeight(WeightVariables::nWeightTypes, 1.f);
  for (int iwt=WeightVariables::wCentral; iwt<WeightVariables::nWeightTypes; iwt++) tree->SetBranchAddress(WeightVariables::getWeightName((WeightVariables::WeightType) iwt), &(genweights.at(iwt)));

  float genMET; tree->SetBranchAddress("genMET", &genMET);
  
  // Gen. particles
  std::vector<bool>* genparticles_isPromptFinalState=nullptr; tree->SetBranchAddress("genparticles_isPromptFinalState", &genparticles_isPromptFinalState);
  std::vector<bool>* genparticles_isPromptDecayed=nullptr; tree->SetBranchAddress("genparticles_isPromptDecayed", &genparticles_isPromptDecayed);
  std::vector<bool>* genparticles_isDirectPromptTauDecayProductFinalState=nullptr; tree->SetBranchAddress("genparticles_isDirectPromptTauDecayProductFinalState", &genparticles_isDirectPromptTauDecayProductFinalState);
  std::vector<bool>* genparticles_isHardProcess=nullptr; tree->SetBranchAddress("genparticles_isHardProcess", &genparticles_isHardProcess);
  std::vector<bool>* genparticles_fromHardProcessFinalState=nullptr; tree->SetBranchAddress("genparticles_fromHardProcessFinalState", &genparticles_fromHardProcessFinalState);
  std::vector<bool>* genparticles_fromHardProcessDecayed=nullptr; tree->SetBranchAddress("genparticles_fromHardProcessDecayed", &genparticles_fromHardProcessDecayed);
  std::vector<bool>* genparticles_isDirectHardProcessTauDecayProductFinalState=nullptr; tree->SetBranchAddress("genparticles_isDirectHardProcessTauDecayProductFinalState", &genparticles_isDirectHardProcessTauDecayProductFinalState);
  std::vector<bool>* genparticles_fromHardProcessBeforeFSR=nullptr; tree->SetBranchAddress("genparticles_fromHardProcessBeforeFSR", &genparticles_fromHardProcessBeforeFSR);
  std::vector<bool>* genparticles_isLastCopy=nullptr; tree->SetBranchAddress("genparticles_isLastCopy", &genparticles_isLastCopy);
  std::vector<bool>* genparticles_isLastCopyBeforeFSR=nullptr; tree->SetBranchAddress("genparticles_isLastCopyBeforeFSR", &genparticles_isLastCopyBeforeFSR);

  std::vector<int>* genparticles_id=nullptr; tree->SetBranchAddress("genparticles_id", &genparticles_id);
  std::vector<int>* genparticles_status=nullptr; tree->SetBranchAddress("genparticles_status", &genparticles_status);

  std::vector<float>* genparticles_px=nullptr; tree->SetBranchAddress("genparticles_px", &genparticles_px);
  std::vector<float>* genparticles_py=nullptr; tree->SetBranchAddress("genparticles_py", &genparticles_py);
  std::vector<float>* genparticles_pz=nullptr; tree->SetBranchAddress("genparticles_pz", &genparticles_pz);
  std::vector<float>* genparticles_mass=nullptr; tree->SetBranchAddress("genparticles_mass", &genparticles_mass);

  // Gen. jets
  std::vector<float>* genjets_pt=nullptr; tree->SetBranchAddress("genjets_pt", &genjets_pt);
  std::vector<float>* genjets_eta=nullptr; tree->SetBranchAddress("genjets_eta", &genjets_eta);
  std::vector<float>* genjets_phi=nullptr; tree->SetBranchAddress("genjets_phi", &genjets_phi);
  std::vector<float>* genjets_mass=nullptr; tree->SetBranchAddress("genjets_mass", &genjets_mass);

  TFile* foutput = TFile::Open(stroutput, "recreate");
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
  plotvars.emplace_back("sqrts", "#sqrt{s_{pp}} (GeV)", 150, 0, 3000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=50; var.max=1000; }
  plotvars.emplace_back("sqrts_initialqg", "#sqrt{s_{pp}} (qg initial states) (GeV)", 150, 0, 3000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=50; var.max=1000; }
  if (isWZ){ auto& var=plotvars.back(); var.nbins=100; var.max=2000; }
  plotvars.emplace_back("sqrts_initialqq", "#sqrt{s_{pp}} (qq initial states) (GeV)", 150, 0, 3000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=50; var.max=1000; }
  if (isWZ){ auto& var=plotvars.back(); var.nbins=50; var.max=1000; }
  plotvars.emplace_back("y", "y_{pp}", 40, -4, 4);
  plotvars.emplace_back("y_initialqg", "y_{pp} (qg initial states)", 40, -4, 4);
  plotvars.emplace_back("y_initialqq", "y_{pp} (qq initial states)", 40, -4, 4);
  plotvars.emplace_back("genmet", "E_{T}^{miss} (GeV)", 40, 0, 400);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=20; var.max=100; }
  plotvars.emplace_back("genmet_initialqg", "E_{T}^{miss} (qg initial states) (GeV)", 40, 0, 400);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=20; var.max=100; }
  plotvars.emplace_back("genmet_initialqq", "E_{T}^{miss} (qq initial states) (GeV)", 40, 0, 400);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=20; var.max=100; }
  plotvars.emplace_back("njets", "N_{jets}", 30, 0, 30);
  plotvars.emplace_back("njets_ptetacuts", "N_{jets} (p_{T}>30 GeV, |#eta|<4.7)", 30, 0, 30);
  plotvars.emplace_back("leadingptj", "Leading p_{T} jet, p_{T} (GeV)", 100, 0, 1000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=40; var.max=400; }
  plotvars.emplace_back("leadingptj_initialqg", "Leading p_{T} jet, p_{T} (qg initial states) (GeV)", 100, 0, 1000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=40; var.max=400; }
  plotvars.emplace_back("leadingptj_initialqq", "Leading p_{T} jet, p_{T} (qq initial states) (GeV)", 100, 0, 1000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=40; var.max=400; }
  plotvars.emplace_back("leadingetaj_eta_ptcut", "Leading |#eta| jet, |#eta| (p_{T}>30 GeV)", 50, 0, 5);
  plotvars.emplace_back("leadingetaj_pt", "Leading |#eta| jet, p_{T}", 40, 0, 200);
  plotvars.emplace_back("jetht", "Jet H_{T} (GeV)", 200, 0, 2000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=60; var.max=600; }
  plotvars.emplace_back("jetht_initialqg", "Jet H_{T} (qg initial states) (GeV)", 200, 0, 2000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=60; var.max=600; }
  plotvars.emplace_back("jetht_initialqq", "Jet H_{T} (qq initial states) (GeV)", 200, 0, 2000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=60; var.max=600; }
  plotvars.emplace_back("nemu", "N_{e,#mu}", 20, 0, 20);
  plotvars.emplace_back("nemu_initialqg", "N_{e,#mu} (qg initial states) (GeV)", 20, 0, 20);
  plotvars.emplace_back("nemu_initialqq", "N_{e,#mu} (qq initial states) (GeV)", 20, 0, 20);
  plotvars.emplace_back("nemu_acceptance", "N_{e,#mu} (p_{T}>5 GeV, |#eta|<2.4)", 12, 0, 12);
  plotvars.emplace_back("nemu_acceptance_initialqg", "N_{e,#mu} (p_{T}>5 GeV, |#eta|<2.4, qg initial states) (GeV)", 12, 0, 12);
  plotvars.emplace_back("nemu_acceptance_initialqq", "N_{e,#mu} (p_{T}>5 GeV, |#eta|<2.4, qq initial states) (GeV)", 12, 0, 12);
  size_t nvars = plotvars.size();

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

    for (int iwt=0; iwt<WeightVariables::nWeightTypes; iwt++) sumCtrWgt[iwt] += genweights[iwt];

    size_t const ngenparts = genparticles_pz->size();
    size_t Nemu=0;
    size_t Nemu_acc=0;
    std::vector<MELAParticle> genparticles; genparticles.reserve(ngenparts);
    for (size_t ip=0; ip<ngenparts; ip++){
      TLorentzVector v; v.SetXYZM(genparticles_px->at(ip), genparticles_py->at(ip), genparticles_pz->at(ip), genparticles_mass->at(ip));
      genparticles.emplace_back(genparticles_id->at(ip), v);
      if (genparticles_status->at(ip)==1 && PDGHelpers::isALepton(genparticles.back().id) && std::abs(genparticles.back().id)!=16){
        Nemu++;
        if (genparticles.back().pt()>5.f && fabs(genparticles.back().eta())<2.5) Nemu_acc++;
      }
    }
    bool isQQ=false;
    bool isQG=false;
    bool isGG=false;
    isGG = (ngenparts>=2 && PDGHelpers::isAGluon(genparticles.at(0).id) && PDGHelpers::isAGluon(genparticles.at(1).id));
    isQQ = (ngenparts>=2 && !PDGHelpers::isAGluon(genparticles.at(0).id) && !PDGHelpers::isAGluon(genparticles.at(1).id));
    isQG = (ngenparts>=2 && !(isGG || isQQ));

    size_t const ngenjets = genjets_pt->size();
    std::vector<MELAParticle> genjets; genjets.reserve(ngenjets);
    for (size_t ip=0; ip<ngenjets; ip++){
      TLorentzVector v; v.SetPtEtaPhiM(genjets_pt->at(ip), genjets_eta->at(ip), genjets_phi->at(ip), genjets_mass->at(ip));
      genjets.emplace_back(0, v);
    }

    for (size_t ivar=0; ivar<nvars; ivar++){
      auto& var=plotvars.at(ivar);
      var.reset();
      if (var.name=="sqrts" && ngenparts>=2) var.setVal((genparticles.at(0).p4+genparticles.at(1).p4).M(), genweights);
      else if (var.name=="sqrts_initialqg" && isQG && ngenparts>=2) var.setVal((genparticles.at(0).p4+genparticles.at(1).p4).M(), genweights);
      else if (var.name=="sqrts_initialqq" && isQQ && ngenparts>=2) var.setVal((genparticles.at(0).p4+genparticles.at(1).p4).M(), genweights);
      else if (var.name=="y" && ngenparts>=2) var.setVal((genparticles.at(0).p4+genparticles.at(1).p4).Rapidity(), genweights);
      else if (var.name=="y_initialqg" && isQG && ngenparts>=2) var.setVal((genparticles.at(0).p4+genparticles.at(1).p4).Rapidity(), genweights);
      else if (var.name=="y_initialqq" && isQQ && ngenparts>=2) var.setVal((genparticles.at(0).p4+genparticles.at(1).p4).Rapidity(), genweights);
      else if (var.name=="genmet") var.setVal(genMET, genweights);
      else if (var.name=="genmet_initialqg" && isQG) var.setVal(genMET, genweights);
      else if (var.name=="genmet_initialqq" && isQQ) var.setVal(genMET, genweights);
      else if (var.name=="njets") var.setVal(ngenjets, genweights);
      else if (var.name=="njets_ptetacuts"){ unsigned int count=0; for (auto const& jet:genjets){ if (jet.pt()>30 && fabs(jet.eta())<4.7) count++; } var.setVal(count, genweights); }
      else if (var.name=="leadingptj" && ngenjets>0) var.setVal(genjets.at(0).pt(), genweights);
      else if (var.name=="leadingptj_initialqg" && isQG && ngenjets>0) var.setVal(genjets.at(0).pt(), genweights);
      else if (var.name=="leadingptj_initialqq" && isQQ && ngenjets>0) var.setVal(genjets.at(0).pt(), genweights);
      else if (var.name=="leadingetaj_eta_ptcut"){ MELAParticle const* theJet=nullptr; for (auto const& jet:genjets){ if (!theJet || fabs(theJet->eta())<fabs(jet.eta())) theJet=&jet; } if (theJet && theJet->pt()>30.) var.setVal(fabs(theJet->eta()), genweights); }
      else if (var.name=="leadingetaj_pt"){ MELAParticle const* theJet=nullptr; for (auto const& jet:genjets){ if (!theJet || fabs(theJet->eta())<fabs(jet.eta())) theJet=&jet; } if (theJet) var.setVal(fabs(theJet->pt()), genweights); }
      else if (var.name=="jetht" && ngenjets>0){ float sumpt=0; for (size_t ijet=0; ijet<ngenjets; ijet++){ sumpt += genjets.at(ijet).pt(); } var.setVal(sumpt, genweights); }
      else if (var.name=="jetht_initialqg" && isQG && ngenjets>0){ float sumpt=0; for (size_t ijet=0; ijet<ngenjets; ijet++){ sumpt += genjets.at(ijet).pt(); } var.setVal(sumpt, genweights); }
      else if (var.name=="jetht_initialqq" && isQQ && ngenjets>0){ float sumpt=0; for (size_t ijet=0; ijet<ngenjets; ijet++){ sumpt += genjets.at(ijet).pt(); } var.setVal(sumpt, genweights); }
      else if (var.name=="nemu") var.setVal(Nemu, genweights);
      else if (var.name=="nemu_initialqg" && isQG) var.setVal(Nemu, genweights);
      else if (var.name=="nemu_initialqq" && isQQ) var.setVal(Nemu, genweights);
      else if (var.name=="nemu_acceptance") var.setVal(Nemu_acc, genweights);
      else if (var.name=="nemu_acceptance_initialqg" && isQG) var.setVal(Nemu_acc, genweights);
      else if (var.name=="nemu_acceptance_initialqq" && isQQ) var.setVal(Nemu_acc, genweights);
      for (int iwt=0; iwt<WeightVariables::nWeightTypes; iwt++) histColl.at(ivar).at(iwt)->Fill(var.val, var.weights.at(iwt));
    }
  }
  for (auto& col:histColl){
    float ymin = std::numeric_limits<float>::max();
    float ymax = std::numeric_limits<float>::min();
    for (auto& h:col){
      if (sumCtrWgt[WeightVariables::wCentral]>0.f) h->Scale(1.f/sumCtrWgt[WeightVariables::wCentral]);
      for (int ix=1; ix<=h->GetNbinsX(); ix++){
        float bc = h->GetBinContent(ix);
        float be = h->GetBinError(ix);
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
          if (doPlot.at(iwt)) hplottable.at(iwt)->Draw((iwt==0 ? "hist" : "histsame"));
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
          if (doPlot.at(iwt)) hplottable.at(iwt)->Draw((iwt==0 ? "hist" : "histsame"));
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
