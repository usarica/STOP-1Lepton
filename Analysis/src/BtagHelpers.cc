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
    fn_btagSF_DeepCSV = "DeepCSV_Moriond17_B_H.csv";              // to be updated
    fn_btagSF_FS_DeepCSV = "fastsim_deepcsv_ttbar_26_1_2017.csv"; // to be updated
    fn_btagSF_CSVv2 = "CSVv2_Moriond17_B_H.csv";               // to be updated
    fn_btagSF_FS_CSVv2 = "fastsim_csvv2_ttbar_26_1_2017.csv";  // to be updated
  }
  else if (theDataYear == 2017){
    fn_btagSF_DeepCSV = "DeepCSV_94XSF_V3_B_F.csv";
    fn_btagSF_FS_DeepCSV = "fastsim_deepcsv_ttbar_26_1_2017.csv"; // to be updated
    fn_btagSF_CSVv2 = "CSVv2_94XSF_V2_B_F.csv";
    fn_btagSF_FS_CSVv2 = "fastsim_csvv2_ttbar_26_1_2017.csv"; // to be updated
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
    WP_DEEPCSV_TIGHT = 0.8001;
    WP_DEEPCSV_MEDIUM = 0.4941;
    WP_DEEPCSV_LOOSE = 0.1522;
    WP_CSVv2_TIGHT = 0.9693;
    WP_CSVv2_MEDIUM = 0.8838;
    WP_CSVv2_LOOSE = 0.5803;
  }
  else if (theDataYear == 2018){
    WP_DEEPCSV_TIGHT = 0.8001;
    WP_DEEPCSV_MEDIUM = 0.4941;
    WP_DEEPCSV_LOOSE = 0.1522;
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
