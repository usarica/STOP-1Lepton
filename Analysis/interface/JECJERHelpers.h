#ifndef JECJERHELPERS_H
#define JECJERHELPERS_H

#include <vector>
#include "Samples.h"


const TString JERDATADIR = STOP1LPKGPATH + "../../cms-jet/JRDatabase/textFiles/";


namespace JECJERHelpers{
  enum JECJERType{
    kAK4,
    kAK8
  };

  TString getJECJERTypeName(JECJERType type);
  TString getJERFilePath(JECJERType type, bool isFastSim);
  TString getJERPtFileName(JECJERType type, bool isFastSim);
  TString getJERPhiFileName(JECJERType type, bool isFastSim);
  TString getJERSFFileName(JECJERType type, bool isFastSim);

}


#endif
