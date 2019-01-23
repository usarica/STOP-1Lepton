#ifndef JECJERHELPERS_H
#define JECJERHELPERS_H

#include <vector>
#include "Samples.h"


const TString JECDATADIR = CMSTASCOREPKGPATH + "Tools/jetcorr/data/";
const TString JERDATADIR = STOP1LPKGPATH + "../../cms-jet/JRDatabase/textFiles/";


namespace JECJERHelpers{
  enum JECJERType{
    kAK4,
    kAK8
  };
  enum JECJERSystematicType{
    sNominal,
    sJECUp, sJECDn,
    sJERUp, sJERDn
  };

  TString getJECJERTypeName(JECJERType type);

  TString getJECFilePath(JECJERType type, bool isMC, bool isFastSim);
  std::vector<TString> getJECFileNames(JECJERType type, bool isMC, bool isFastSim);
  TString getJECUncertaintyFileName(JECJERType type, bool isMC, bool isFastSim);

  TString getJERFilePath(JECJERType type, bool isMC, bool isFastSim);
  TString getJERPtFileName(JECJERType type, bool isMC, bool isFastSim);
  TString getJERPhiFileName(JECJERType type, bool isMC, bool isFastSim);
  TString getJERSFFileName(JECJERType type, bool isFastSim);

}


#endif
