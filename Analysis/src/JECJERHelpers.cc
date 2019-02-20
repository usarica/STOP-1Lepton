#include <cassert>
#include <unordered_map>
#include "StdExtensions.h"
#include "JECJERHelpers.h"
#include "HostHelpersCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace MELAStreamHelpers;


TString JECJERHelpers::getJECJERTypeName(JECJERType type){
  if (type==kAK4) return "AK4PFchs";
  else if (type==kAK8 && theDataVersion == kCMSSW_8_0_X) return "AK8PFchs";
  else if (type==kAK8) return "AK8PFPuppi";
  else{ MELAerr << "JECJERHelpers::getJECJERTypeName: JECJERType " << type << " is not implemented!" << endl; assert(0); return ""; }
}

TString JECJERHelpers::getJECFilePath(JECJERType /*type*/, bool isMC, bool isFastSim){
  unordered_map<TString, TString> eraMap;
  if (theDataYear == 2016 && theDataVersion == kCMSSW_8_0_X){ // Does not have AK8 corrections
    eraMap["2016B"] = eraMap["2016C"] = eraMap["2016D"] = "Summer16_23Sep2016BCDV4_DATA";
    eraMap["2016E"] = eraMap["2016F"] = "Summer16_23Sep2016EFV4_DATA";
    eraMap["2016G"] = "Summer16_23Sep2016GV4_DATA";
    eraMap["2016H"] = "Summer16_23Sep2016HV4_DATA";
    eraMap["MC_noFS"] = "Summer16_23Sep2016V4_MC";
    eraMap["MC_FS"] = "Spring16_FastSimV1";
  }
  else if (theDataYear == 2016 && theDataVersion == kCMSSW_9_4_X){
    eraMap["2016B"] = eraMap["2016C"] = eraMap["2016D"] = "Summer16_07Aug2017BCD_V11_DATA";
    eraMap["2016E"] = eraMap["2016F"] = "Summer16_07Aug2017EF_V11_DATA";
    eraMap["2016G"] = eraMap["2016H"] = "Summer16_07Aug2017GH_V11_DATA";
    eraMap["MC_noFS"] = "Summer16_07Aug2017_V11_MC";
    eraMap["MC_FS"] = "Spring16_FastSimV1";
  }
  else if (theDataYear == 2017 && theDataVersion == kCMSSW_9_4_X){
    eraMap["2017B"] = "Fall17_17Nov2017B_V32_DATA";
    eraMap["2017C"] = "Fall17_17Nov2017C_V32_DATA";
    eraMap["2017D"] = eraMap["2017E"] = "Fall17_17Nov2017DE_V32_DATA";
    eraMap["2017F"] = "Fall17_17Nov2017F_V32_DATA";
    eraMap["MC_noFS"] = eraMap["MC_FS"] = "Fall17_17Nov2017_V32_MC";
  }
  else if (theDataYear == 2018 && theDataVersion == kCMSSW_10_X){ // FIXME: 2018 to be updated!
    eraMap["2018A"] = eraMap["2018B"] = eraMap["2018C"] = eraMap["2018D"] = "Fall17_17Nov2017F_V32_DATA";
    eraMap["MC_noFS"] = eraMap["MC_FS"] = "Fall17_17Nov2017_V32_MC";
  }
  else{
    MELAerr << "JECJERHelpers::getJECFilePath: Data year " << theDataYear << " and data version " << theDataVersion << " are not recognized. Aborting..." << endl;
    assert(0);
  }
  if (isMC) return (isFastSim ? eraMap["MC_FS"] : eraMap["MC_noFS"]);
  else{
    unordered_map<TString, TString>::const_iterator it = eraMap.find(theDataPeriod);
    if (it==eraMap.cend()){
      if (SampleHelpers::testDataPeriodIsLikeData()){
        MELAerr << "JECJERHelpers::getJECFilePath: Era map does not contain the data period " << theDataPeriod << " for year " << theDataYear << " and data version " << theDataVersion << ". Aborting..." << endl;
        assert(0);
      }
      return "";
    }
    else return it->second;
  }
}
std::vector<TString> JECJERHelpers::getJECFileNames(JECJERType type, bool isMC, bool isFastSim){
  std::vector<TString> res;
  if (type == kAK8 && theDataYear == 2016 && theDataVersion == kCMSSW_8_0_X) return res;
  TString jecVer = getJECFilePath(type, isMC, isFastSim);
  if (jecVer=="") return res;
  TString jetType = getJECJERTypeName(type);

  TString sprep = JECDATADIR;
  if (theDataYear >= 2016) sprep += "run2_25ns/";
  else if (theDataYear == 2015) sprep += "run2_50ns/";

  res.push_back(sprep+jecVer+"/"+jecVer+"_L1FastJet_"+jetType+".txt");
  res.push_back(sprep+jecVer+"/"+jecVer+"_L2Relative_"+jetType+".txt");
  res.push_back(sprep+jecVer+"/"+jecVer+"_L3Absolute_"+jetType+".txt");
  if (!(isMC && isFastSim) || type == kAK8) res.push_back(sprep+jecVer+"/"+jecVer+"_L2L3Residual_"+jetType+".txt");
  return res;
}
TString JECJERHelpers::getJECUncertaintyFileName(JECJERType type, bool isMC, bool isFastSim){
  TString res;
  if (type == kAK8 && theDataYear == 2016 && theDataVersion == kCMSSW_8_0_X) return res;
  TString jecVer = getJECFilePath(type, isMC, isFastSim);
  if (jecVer=="") return res;
  TString jetType = getJECJERTypeName(type);

  TString sprep = JECDATADIR;
  if (theDataYear >= 2016) sprep += "run2_25ns/";
  else if (theDataYear == 2015) sprep += "run2_50ns/";

  res = sprep+jecVer+"/"+jecVer+"_Uncertainty_"+jetType+".txt";
  return res;
}

TString JECJERHelpers::getJERFilePath(JECJERType /*type*/, bool isMC, bool /*isFastSim*/){
  TString res;
  TString dataOrMC;
  if (isMC) dataOrMC="MC";
  else dataOrMC="DATA";
  switch (theDataYear){
  case 2016:
    res = "Summer16_25nsV1_";
    break;
  case 2017:
  case 2018: // FIXME: Update for the 2018 recipe needed
    res = "Fall17_V3_";
    break;
  default:
    MELAerr << "JECJERHelpers::getJERFilePath: Year " << theDataYear << " is unknown." << endl;
    assert(0);
    return "";
  }
  res += dataOrMC;
  return res;
}
TString JECJERHelpers::getJERPtFileName(JECJERType type, bool isMC, bool isFastSim){
  TString res = JERDATADIR;
  TString fext = getJERFilePath(type, isMC, isFastSim);
  res += fext + "/" + fext + "_PtResolution_" + getJECJERTypeName(type) + ".txt";
  if (!HostHelpers::FileExists(res.Data())){
    MELAerr << "JECJERHelpers::getJERPtFileName: File " << res << " does not exist! Aborting..." << endl;
    assert(0);
  }
  return res;
}
TString JECJERHelpers::getJERPhiFileName(JECJERType type, bool isMC, bool isFastSim){
  TString res = JERDATADIR;
  TString fext = getJERFilePath(type, isMC, isFastSim);
  res += fext + "/" + fext + "_PhiResolution_" + getJECJERTypeName(type) + ".txt";
  if (!HostHelpers::FileExists(res.Data())){
    MELAerr << "JECJERHelpers::getJERPhiFileName: File " << res << " does not exist! Aborting..." << endl;
    assert(0);
  }
  return res;
}
TString JECJERHelpers::getJERSFFileName(JECJERType type, bool isFastSim){
  TString res = JERDATADIR;
  TString fext = getJERFilePath(type, true, isFastSim);
  res += fext + "/" + fext + "_SF_" + getJECJERTypeName(type) + ".txt";
  if (!HostHelpers::FileExists(res.Data())){
    MELAerr << "JECJERHelpers::getJERSFFileName: File " << res << " does not exist! Aborting..." << endl;
    assert(0);
  }
  return res;
}
