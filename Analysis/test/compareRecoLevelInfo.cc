#include "common_includes.h"
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

void compareRecoLevelInfo(TString fname){
  gStyle->SetOptStat(0);

  SampleHelpers::theDataYear=2018;
  SampleHelpers::theDataPeriod="2018";
  SampleHelpers::theDataVersion=SampleHelpers::kCMSSW_10_X;
  float const btagWP = BtagHelpers::getBtagWP(BtagHelpers::kDeepCSV_Medium);
  float const totalLumi=150.f*1000.f;

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
  TString stroutputcore=strindir+"/plots/reco/";
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

  // Gen. mothers
  DATA_CLASS(std::vector<int>, genparticles_id)
  DATA_CLASS(std::vector<float>, genparticles_px)
  DATA_CLASS(std::vector<float>, genparticles_py)
  DATA_CLASS(std::vector<float>, genparticles_pz)
  DATA_CLASS(std::vector<float>, genparticles_mass)
  // MET
  DATA_SIMPLE(float, pfmet)
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
  DATA_CLASS(std::vector<int>, ak4jets_npfcands)
  DATA_CLASS(std::vector<float>, ak4jets_ptDistribution)
  DATA_CLASS(std::vector<float>, ak4jets_axis1)
  DATA_CLASS(std::vector<float>, ak4jets_axis2)
  DATA_CLASS(std::vector<float>, ak4jets_deepCSVb)
  DATA_CLASS(std::vector<float>, ak4jets_deepCSVbb)
  DATA_CLASS(std::vector<float>, ak4jets_deepCSVc)
  DATA_CLASS(std::vector<float>, ak4jets_deepCSVl)
  // AK8Jets
  DATA_CLASS(std::vector<float>, ak8jets_pt)
  DATA_CLASS(std::vector<float>, ak8jets_eta)
  DATA_CLASS(std::vector<float>, ak8jets_phi)
  DATA_CLASS(std::vector<float>, ak8jets_mass)
  DATA_CLASS(std::vector<long long>, ak8jets_selectionBits)
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
  plotvars.emplace_back("recomet", "E_{T}^{miss} (GeV)", 40, 0, 400);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=20; var.max=100; }
  plotvars.emplace_back("recomet_initialqq", "E_{T}^{miss} (qq initial states) (GeV)", 40, 0, 400);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=20; var.max=100; }
  plotvars.emplace_back("recomet_initialqg", "E_{T}^{miss} (qg initial states) (GeV)", 40, 0, 400);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=20; var.max=100; }
  plotvars.emplace_back("recomet_initialgg", "E_{T}^{miss} (gg initial states) (GeV)", 40, 0, 400);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=20; var.max=100; }

  plotvars.emplace_back("nak4jets", "N_{ak4}", 15, 0, 15);
  plotvars.emplace_back("nak4jets_initialqq", "N_{ak4} (qq initial states)", 15, 0, 15);
  plotvars.emplace_back("nak4jets_initialqg", "N_{ak4} (qg initial states)", 15, 0, 15);
  plotvars.emplace_back("nak4jets_initialgg", "N_{ak4} (gg initial states)", 15, 0, 15);
  plotvars.emplace_back("nak4jets_metcut", "N_{ak4} (E_{T}^{miss}>150 GeV)", 15, 0, 15);
  plotvars.emplace_back("nak4jets_metcut_initialqq", "N_{ak4} (E_{T}^{miss}>150 GeV, qq initial states)", 15, 0, 15);
  plotvars.emplace_back("nak4jets_metcut_initialqg", "N_{ak4} (E_{T}^{miss}>150 GeV, qg initial states)", 15, 0, 15);
  plotvars.emplace_back("nak4jets_metcut_initialgg", "N_{ak4} (E_{T}^{miss}>150 GeV, gg initial states)", 15, 0, 15);
  plotvars.emplace_back("nak4jets_btagged", "N_{ak4} (b-tagged)", 8, 0, 8);
  plotvars.emplace_back("nak4jets_btagged_initialqq", "N_{ak4} (b-tagged, qq initial states)", 8, 0, 8);
  plotvars.emplace_back("nak4jets_btagged_initialqg", "N_{ak4} (b-tagged, qg initial states)", 8, 0, 8);
  plotvars.emplace_back("nak4jets_btagged_initialgg", "N_{ak4} (b-tagged, gg initial states)", 8, 0, 8);
  plotvars.emplace_back("nak4jets_btagged_metcut", "N_{ak4} (b-tagged, E_{T}^{miss}>150 GeV)", 8, 0, 8);
  plotvars.emplace_back("nak4jets_btagged_metcut_initialqq", "N_{ak4} (b-tagged, E_{T}^{miss}>150 GeV, qq initial states)", 8, 0, 8);
  plotvars.emplace_back("nak4jets_btagged_metcut_initialqg", "N_{ak4} (b-tagged, E_{T}^{miss}>150 GeV, qg initial states)", 8, 0, 8);
  plotvars.emplace_back("nak4jets_btagged_metcut_initialgg", "N_{ak4} (b-tagged, E_{T}^{miss}>150 GeV, gg initial states)", 8, 0, 8);

  plotvars.emplace_back("leadingptak4jet_pt", "Leading p_{T} ak4-jet, p_{T} (GeV)", 50, 0, 1000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=20; var.max=400; }
  plotvars.emplace_back("leadingptak4jet_pt_initialqq", "Leading p_{T} ak4-jet, p_{T} (qq initial states) (GeV)", 50, 0, 1000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=20; var.max=400; }
  plotvars.emplace_back("leadingptak4jet_pt_initialqg", "Leading p_{T} ak4-jet, p_{T} (qg initial states) (GeV)", 50, 0, 1000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=20; var.max=400; }
  plotvars.emplace_back("leadingptak4jet_pt_initialgg", "Leading p_{T} ak4-jet, p_{T} (gg initial states) (GeV)", 50, 0, 1000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=20; var.max=400; }
  plotvars.emplace_back("leadingptak4jet_pt_metcut", "Leading p_{T} ak4-jet, p_{T} (E_{T}^{miss}>150 GeV) (GeV)", 40, 0, 1000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=16; var.max=400; }
  plotvars.emplace_back("leadingptak4jet_pt_metcut_initialqq", "Leading p_{T} ak4-jet, p_{T} (E_{T}^{miss}>150 GeV, qq initial states) (GeV)", 40, 0, 1000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=16; var.max=400; }
  plotvars.emplace_back("leadingptak4jet_pt_metcut_initialqg", "Leading p_{T} ak4-jet, p_{T} (E_{T}^{miss}>150 GeV, qg initial states) (GeV)", 40, 0, 1000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=16; var.max=400; }
  plotvars.emplace_back("leadingptak4jet_pt_metcut_initialgg", "Leading p_{T} ak4-jet, p_{T} (E_{T}^{miss}>150 GeV, gg initial states) (GeV)", 40, 0, 1000);
  if (isDYLLJets){ auto& var=plotvars.back(); var.nbins=16; var.max=400; }

  plotvars.emplace_back("leadingptak4jet_axis2", "Leading p_{T} ak4-jet, axis_{2}", 10, 0, 0.2);
  plotvars.emplace_back("leadingptak4jet_multiplicity", "Leading p_{T} ak4-jet, N_{PFcands}", 50, 0, 50);
  plotvars.emplace_back("leadingptak4jet_ptD", "Leading p_{T} ak4-jet, p_{T}^{D}", 20, 0, 1);

  plotvars.emplace_back("ntftops", "N_{tops}", 5, 0, 5);
  plotvars.emplace_back("ntftops_initialqq", "N_{tops} (qq initial states)", 5, 0, 5);
  plotvars.emplace_back("ntftops_initialqg", "N_{tops} (qg initial states)", 5, 0, 5);
  plotvars.emplace_back("ntftops_initialgg", "N_{tops} (gg initial states)", 5, 0, 5);
  plotvars.emplace_back("ntftops_metcut", "N_{tops} (E_{T}^{miss}>150 GeV)", 5, 0, 5);
  plotvars.emplace_back("ntftops_metcut_initialqq", "N_{tops} (E_{T}^{miss}>150 GeV, qq initial states)", 5, 0, 5);
  plotvars.emplace_back("ntftops_metcut_initialqg", "N_{tops} (E_{T}^{miss}>150 GeV, qg initial states)", 5, 0, 5);
  plotvars.emplace_back("ntftops_metcut_initialgg", "N_{tops} (E_{T}^{miss}>150 GeV, gg initial states)", 5, 0, 5);

  plotvars.emplace_back("highestdisctftop_mass", "m_{3j} with highest top disc. (GeV)", 50, 100, 250);
  plotvars.emplace_back("highestdisctftop_mass_initialqq", "m_{3j} with highest top disc. (qq initial states) (GeV)", 50, 100, 250);
  plotvars.emplace_back("highestdisctftop_mass_initialqg", "m_{3j} with highest top disc. (qg initial states) (GeV)", 50, 100, 250);
  plotvars.emplace_back("highestdisctftop_mass_initialgg", "m_{3j} with highest top disc. (gg initial states) (GeV)", 50, 100, 250);
  plotvars.emplace_back("highestdisctftop_mass_metcut", "m_{3j} with highest top disc. (E_{T}^{miss}>150 GeV) (GeV)", 50, 100, 250);
  plotvars.emplace_back("highestdisctftop_mass_metcut_initialqq", "m_{3j} with highest top disc. (E_{T}^{miss}>150 GeV, qq initial states) (GeV)", 50, 100, 250);
  plotvars.emplace_back("highestdisctftop_mass_metcut_initialqg", "m_{3j} with highest top disc. (E_{T}^{miss}>150 GeV, qg initial states) (GeV)", 50, 100, 250);
  plotvars.emplace_back("highestdisctftop_mass_metcut_initialgg", "m_{3j} with highest top disc. (E_{T}^{miss}>150 GeV, gg initial states) (GeV)", 50, 100, 250);

  plotvars.emplace_back("nelectrons", "N_{e}", 6, 0, 6);
  plotvars.emplace_back("nelectrons_initialqq", "N_{e} (qq initial states)", 6, 0, 6);
  plotvars.emplace_back("nelectrons_initialqg", "N_{e} (qg initial states)", 6, 0, 6);
  plotvars.emplace_back("nelectrons_initialgg", "N_{e} (gg initial states)", 6, 0, 6);
  plotvars.emplace_back("nelectrons_metcut", "N_{e} (E_{T}^{miss}>150 GeV)", 6, 0, 6);
  plotvars.emplace_back("nelectrons_metcut_initialqq", "N_{e} (E_{T}^{miss}>150 GeV, qq initial states)", 6, 0, 6);
  plotvars.emplace_back("nelectrons_metcut_initialqg", "N_{e} (E_{T}^{miss}>150 GeV, qg initial states)", 6, 0, 6);
  plotvars.emplace_back("nelectrons_metcut_initialgg", "N_{e} (E_{T}^{miss}>150 GeV, gg initial states)", 6, 0, 6);

  plotvars.emplace_back("nmuons", "N_{#mu}", 6, 0, 6);
  plotvars.emplace_back("nmuons_initialqq", "N_{#mu} (qq initial states)", 6, 0, 6);
  plotvars.emplace_back("nmuons_initialqg", "N_{#mu} (qg initial states)", 6, 0, 6);
  plotvars.emplace_back("nmuons_initialgg", "N_{#mu} (gg initial states)", 6, 0, 6);
  plotvars.emplace_back("nmuons_metcut", "N_{#mu} (E_{T}^{miss}>150 GeV)", 6, 0, 6);
  plotvars.emplace_back("nmuons_metcut_initialqq", "N_{#mu} (E_{T}^{miss}>150 GeV, qq initial states)", 6, 0, 6);
  plotvars.emplace_back("nmuons_metcut_initialqg", "N_{#mu} (E_{T}^{miss}>150 GeV, qg initial states)", 6, 0, 6);
  plotvars.emplace_back("nmuons_metcut_initialgg", "N_{#mu} (E_{T}^{miss}>150 GeV, gg initial states)", 6, 0, 6);

  /*
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
  */
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

    size_t const ngenparts = genparticles_pz->size();
    std::vector<MELAParticle> genparticles; genparticles.reserve(ngenparts);
    for (size_t ip=0; ip<ngenparts; ip++){
      TLorentzVector v; v.SetXYZM(genparticles_px->at(ip), genparticles_py->at(ip), genparticles_pz->at(ip), genparticles_mass->at(ip));
      genparticles.emplace_back(genparticles_id->at(ip), v);
    }
    bool isQQ=false;
    bool isQG=false;
    bool isGG=false;
    isGG = (ngenparts>=2 && PDGHelpers::isAGluon(genparticles.at(0).id) && PDGHelpers::isAGluon(genparticles.at(1).id));
    isQQ = (ngenparts>=2 && !PDGHelpers::isAGluon(genparticles.at(0).id) && !PDGHelpers::isAGluon(genparticles.at(1).id));
    isQG = (ngenparts>=2 && !(isGG || isQQ));

    size_t const nelectrons = electrons_pt->size();
    std::vector<MELAParticle> electrons; electrons.reserve(nelectrons);
    for (size_t ip=0; ip<nelectrons; ip++){
      TLorentzVector v; v.SetPtEtaPhiM(electrons_pt->at(ip), electrons_eta->at(ip), electrons_phi->at(ip), electrons_mass->at(ip));
      electrons.emplace_back(0, v);
    }

    size_t const nmuons = muons_pt->size();
    std::vector<MELAParticle> muons; muons.reserve(nmuons);
    for (size_t ip=0; ip<nmuons; ip++){
      TLorentzVector v; v.SetPtEtaPhiM(muons_pt->at(ip), muons_eta->at(ip), muons_phi->at(ip), muons_mass->at(ip));
      muons.emplace_back(0, v);
    }

    size_t const nak4jets = ak4jets_pt->size();
    std::vector<MELAParticle> ak4jets; ak4jets.reserve(nak4jets);
    for (size_t ip=0; ip<nak4jets; ip++){
      TLorentzVector v; v.SetPtEtaPhiM(ak4jets_pt->at(ip), ak4jets_eta->at(ip), ak4jets_phi->at(ip), ak4jets_mass->at(ip));
      ak4jets.emplace_back(0, v);
    }

    size_t const nak8jets = ak8jets_pt->size();
    std::vector<MELAParticle> ak8jets; ak8jets.reserve(nak8jets);
    for (size_t ip=0; ip<nak8jets; ip++){
      TLorentzVector v; v.SetPtEtaPhiM(ak8jets_pt->at(ip), ak8jets_eta->at(ip), ak8jets_phi->at(ip), ak8jets_mass->at(ip));
      ak8jets.emplace_back(0, v);
    }

    size_t const ntftops = tftops_pt->size();
    std::vector<MELAParticle> tftops; tftops.reserve(ntftops);
    for (size_t ip=0; ip<ntftops; ip++){
      TLorentzVector v; v.SetPtEtaPhiM(tftops_pt->at(ip), tftops_eta->at(ip), tftops_phi->at(ip), tftops_mass->at(ip));
      tftops.emplace_back(0, v);
    }

    for (size_t ivar=0; ivar<nvars; ivar++){
      auto& var=plotvars.at(ivar);
      var.reset();

      if (var.name.Contains("metcut") && pfmet<150.f) continue;
      if (var.name.Contains("initialqq") && !isQQ) continue;
      if (var.name.Contains("initialqg") && !isQG) continue;
      if (var.name.Contains("initialgg") && !isGG) continue;

      if (var.name.BeginsWith("recomet")) var.setVal(pfmet, genweights);
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
            if (var.name.BeginsWith("leadingptak4jet_axis2")) v=ak4jets_axis2->at(ip);
            else if (var.name.BeginsWith("leadingptak4jet_ptD")) v=ak4jets_ptDistribution->at(ip);
            else if (var.name.BeginsWith("leadingptak4jet_pt")) v=jet.pt();
            else if (var.name.BeginsWith("leadingptak4jet_multiplicity")) v=ak4jets_npfcands->at(ip);
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
