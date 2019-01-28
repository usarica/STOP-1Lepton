#include <cassert>
#include "BtagHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace MELAStreamHelpers;



TString BtagHelpers::getBtagSFFileName(BtagWPType type, bool isFastSim){
  TString res = BTAGSFPATH;
  if (isFastSim) res += "run2_fastsim/";
  else if (theDataYear >= 2016) res += "run2_25ns/";
  else if (theDataYear == 2015) res += "run2_50ns/";

  // From cmstas/CORE/Config.cc
  TString fn_btagSF_DeepCSV, fn_btagSF_FS_DeepCSV, fn_btagSF_CSVv2, fn_btagSF_FS_CSVv2;
  if (theDataYear == 2016 && theDataVersion == kCMSSW_8_0_X){
    fn_btagSF_DeepCSV = "DeepCSV_Moriond17_B_H.csv";
    fn_btagSF_FS_DeepCSV = "fastsim_deepcsv_ttbar_26_1_2017.csv";
    fn_btagSF_CSVv2 = "CSVv2_Moriond17_B_H.csv";
    fn_btagSF_FS_CSVv2 = "fastsim_csvv2_ttbar_26_1_2017.csv";
  }
  else if (theDataYear == 2016 && theDataVersion == kCMSSW_9_4_X){
    fn_btagSF_DeepCSV = "DeepCSV_Moriond17_B_H.csv";
    fn_btagSF_FS_DeepCSV = "fastsim_deepcsv_ttbar_26_1_2017.csv";
    fn_btagSF_CSVv2 = "CSVv2_Moriond17_B_H.csv";
    fn_btagSF_FS_CSVv2 = "fastsim_csvv2_ttbar_26_1_2017.csv";
  }
  else if (theDataYear == 2017){
    fn_btagSF_DeepCSV = "DeepCSV_94XSF_V3_B_F.csv";
    fn_btagSF_FS_DeepCSV = "fastsim_deepcsv_ttbar_26_1_2017.csv";
    fn_btagSF_CSVv2 = "CSVv2_94XSF_V2_B_F.csv";
    fn_btagSF_FS_CSVv2 = "fastsim_csvv2_ttbar_26_1_2017.csv";
  }
  else if (theDataYear == 2018){
    fn_btagSF_DeepCSV = "DeepCSV_94XSF_V3_B_F.csv";
    fn_btagSF_FS_DeepCSV = "fastsim_deepcsv_ttbar_26_1_2017.csv";
    fn_btagSF_CSVv2 = "CSVv2_94XSF_V2_B_F.csv";
    fn_btagSF_FS_CSVv2 = "fastsim_csvv2_ttbar_26_1_2017.csv";
  }
  else{
    MELAerr << "BtagHelpers::getBtagSFFileName: Cannot determine the b tag SF file name for year " << theDataYear << " and version " << theDataVersion << ". Aborting..." << endl;
    assert(0);
  }
  switch (type){
  case kCSVv2_Loose:
  case kCSVv2_Medium:
  case kCSVv2_Tight:
    res += (!isFastSim ? fn_btagSF_CSVv2 : fn_btagSF_FS_CSVv2);
    break;
  case kDeepCSV_Loose:
  case kDeepCSV_Medium:
  case kDeepCSV_Tight:
    res += (!isFastSim ? fn_btagSF_DeepCSV : fn_btagSF_FS_DeepCSV);
    break;
  default:
    MELAerr << "BtagHelpers::getBtagSFFileName: No implementation for b tag WP " << type << ". Aborting..." << endl;
    assert(0);
    break;
  }

  return res;
}

TString BtagHelpers::getBtagEffFileName(BtagWPType /*type*/, bool isFastSim){
  TString res = BTAGSFPATH;
  if (isFastSim) res += "run2_fastsim/";
  else if (theDataYear >= 2016) res += "run2_25ns/";
  else if (theDataYear == 2015) res += "run2_50ns/";

  if (isFastSim){
    res += "btageff__SMS-T1bbbb-T1qqqq_25ns_Moriond17.root";
  }
  else{
    if (theDataYear == 2016) res += "btageff__ttbar_powheg_pythia8_25ns_Moriond17_deepCSV.root";
    else if (theDataYear >= 2017) res += "btageff__ttbar_amc_94X_deepCSV.root";
  }

  return res;
}

std::vector<TString> BtagHelpers::getBtagEffHistogramNames(BtagWPType type, bool /*isFastSim*/){
  const TString s_loose_btag_eff_b = "h2_BTaggingEff_csv_loose_Eff_b";
  const TString s_loose_btag_eff_c = "h2_BTaggingEff_csv_loose_Eff_c";
  const TString s_loose_btag_eff_udsg = "h2_BTaggingEff_csv_loose_Eff_udsg";
  const TString s_medium_btag_eff_b = "h2_BTaggingEff_csv_med_Eff_b";
  const TString s_medium_btag_eff_c = "h2_BTaggingEff_csv_med_Eff_c";
  const TString s_medium_btag_eff_udsg = "h2_BTaggingEff_csv_med_Eff_udsg";
  const TString s_tight_btag_eff_b = "h2_BTaggingEff_csv_tight_Eff_b";
  const TString s_tight_btag_eff_c = "h2_BTaggingEff_csv_tight_Eff_c";
  const TString s_tight_btag_eff_udsg = "h2_BTaggingEff_csv_tight_Eff_udsg";

  std::vector<TString> res;
  switch (type){
  case kCSVv2_Loose:
  case kDeepCSV_Loose:
    res = std::vector<TString>{ s_loose_btag_eff_b, s_loose_btag_eff_c, s_loose_btag_eff_udsg };
    break;
  case kCSVv2_Medium:
  case kDeepCSV_Medium:
    res = std::vector<TString>{ s_medium_btag_eff_b, s_medium_btag_eff_c, s_medium_btag_eff_udsg };
    break;
  case kCSVv2_Tight:
  case kDeepCSV_Tight:
    res = std::vector<TString>{ s_tight_btag_eff_b, s_tight_btag_eff_c, s_tight_btag_eff_udsg };
    break;
  default:
    MELAerr << "BtagHelpers::getBtagEffHistogramNames: No implementation for b tag histogram names for WP type " << type << ". Aborting..." << endl;
    assert(0);
  }
  return res;
}

float BtagHelpers::getBtagWP(BtagWPType type){
  // From cmstas/CORE/Config.cc
  float WP_DEEPCSV_TIGHT, WP_DEEPCSV_MEDIUM, WP_DEEPCSV_LOOSE, WP_CSVv2_TIGHT, WP_CSVv2_MEDIUM, WP_CSVv2_LOOSE;
  if (theDataYear == 2016 && theDataVersion == kCMSSW_8_0_X){
    WP_DEEPCSV_TIGHT = 0.8958;
    WP_DEEPCSV_MEDIUM = 0.6324;
    WP_DEEPCSV_LOOSE = 0.2219;
    WP_CSVv2_TIGHT = 0.9535;
    WP_CSVv2_MEDIUM = 0.8484;
    WP_CSVv2_LOOSE = 0.5426;
  }
  else if (theDataYear == 2016 && theDataVersion == kCMSSW_9_4_X){
    WP_DEEPCSV_TIGHT = 0.8958;
    WP_DEEPCSV_MEDIUM = 0.6324;
    WP_DEEPCSV_LOOSE = 0.2219;
    WP_CSVv2_TIGHT = 0.9535;
    WP_CSVv2_MEDIUM = 0.8484;
    WP_CSVv2_LOOSE = 0.5426;
  }
  else if (theDataYear == 2017){
    // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    WP_DEEPCSV_TIGHT = 0.8001;
    WP_DEEPCSV_MEDIUM = 0.4941;
    WP_DEEPCSV_LOOSE = 0.1522;
    WP_CSVv2_TIGHT = 0.9693;
    WP_CSVv2_MEDIUM = 0.8838;
    WP_CSVv2_LOOSE = 0.5803;
  }
  else if (theDataYear == 2018){
    // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    WP_DEEPCSV_TIGHT = 0.7527;
    WP_DEEPCSV_MEDIUM = 0.4184;
    WP_DEEPCSV_LOOSE = 0.1241;
    // FIXME: Not really supported yet, so keep 2017 values
    WP_CSVv2_TIGHT = 0.9693;
    WP_CSVv2_MEDIUM = 0.8838;
    WP_CSVv2_LOOSE = 0.5803;
  }
  else{
    MELAerr << "BtagHelpers::getBtagWP: Cannot determine the b tag SF file name for year " << theDataYear << " and version " << theDataVersion << ". Aborting..." << endl;
    assert(0);
  }
  switch (type){
  case kCSVv2_Loose:
    return WP_CSVv2_LOOSE;
  case kCSVv2_Medium:
    return WP_CSVv2_MEDIUM;
  case kCSVv2_Tight:
    return WP_CSVv2_TIGHT;
  case kDeepCSV_Loose:
    return WP_DEEPCSV_LOOSE;
  case kDeepCSV_Medium:
    return WP_DEEPCSV_MEDIUM;
  case kDeepCSV_Tight:
    return WP_DEEPCSV_TIGHT;
  default:
    MELAerr << "BtagHelpers::getBtagWP: No implementation for b tag WP " << type << ". Aborting..." << endl;
    assert(0);
    return 0; // Just return to avoid warnings
  }
}
