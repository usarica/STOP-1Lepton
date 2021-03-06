#ifndef BTAGHELPERS_H
#define BTAGHELPERS_H

#include <vector>
#include "Samples.h"


const TString BTAGSFPATH = CMSTASCOREPKGPATH + "Tools/btagsf/data/";


namespace BtagHelpers{
  enum BtagWPType{
    kCSVv2_Loose,
    kCSVv2_Medium,
    kCSVv2_Tight,
    kDeepCSV_Loose,
    kDeepCSV_Medium,
    kDeepCSV_Tight
  };

  TString getBtagSFFileName(BtagWPType type, bool isFastSim);
  TString getBtagEffFileName(BtagWPType type, bool isFastSim);
  std::vector<TString> getBtagEffHistogramNames(BtagWPType type, bool isFastSim);
  float getBtagWP(BtagWPType type);

}


#endif
