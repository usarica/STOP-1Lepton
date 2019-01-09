#include <cassert>
#include "JECJERHelpers.h"
#include "HostHelpersCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace MELAStreamHelpers;


TString JECJERHelpers::getJECJERTypeName(JECJERType type){
  if (type==kAK4) return "AK4PFchs";
  else if (type==kAK8) return "AK8PFchs";
  else{ MELAerr << "JECJERHelpers::getJECJERTypeName: JECJERType " << type << " is not implemented!" << endl; assert(0); return ""; }
}
TString JECJERHelpers::getJERFilePath(JECJERType type, bool /*isFastSim*/){
  switch (theDataYear){
  case 2016:
    return "Summer16_25nsV1";
  case 2017:
  case 2018: // FIXME: Update for the 2018 recipe needed
    return "Fall17_V3";
  default:
    MELAerr << "JECJERHelpers::getJERFilePath: Year " << theDataYear << " is unknown." << endl;
    assert(0);
    return "";
  }
}
TString JECJERHelpers::getJERPtFileName(JECJERType type, bool isFastSim){
  TString res = JERDATADIR;
  TString fext = getJERFilePath(type, isFastSim);
  res += fext + "/" + fext + "_PtResolution_" + getJECJERTypeName(type) + ".txt";
  assert(HostHelpers::FileExists(res.Data()));
  return res;
}
TString JECJERHelpers::getJERPhiFileName(JECJERType type, bool isFastSim){
  TString res = JERDATADIR;
  TString fext = getJERFilePath(type, isFastSim);
  res += fext + "/" + fext + "_PhiResolution_" + getJECJERTypeName(type) + ".txt";
  assert(HostHelpers::FileExists(res.Data()));
  return res;
}
TString JECJERHelpers::getJERSFFileName(JECJERType type, bool isFastSim){
  TString res = JERDATADIR;
  TString fext = getJERFilePath(type, isFastSim);
  res += fext + "/" + fext + "_SF_" + getJECJERTypeName(type) + ".txt";
  assert(HostHelpers::FileExists(res.Data()));
  return res;
}
