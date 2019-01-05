#ifndef TFTOPTAGGERHELPERS_H
#define TFTOPTAGGERHELPERS_H

#include <TopTagger/TopTagger/interface/TopTagger.h>
#include "Samples.h"
#include "AK4JetObject.h"
#include "TFTopObject.h"


const TString TFTOPTAGGERPKGDATAPATH = STOP1LPKGDATAPATH + "/TFTopTagger/";


namespace TFTopTaggerHelpers{
  extern TopTagger tagger;
  extern bool TopTaggerCfg_firsttime;

  bool setTFTopTaggerConfigFile();

  std::vector<TFTopObject*> getTopsFromResolvedJets(std::vector<AK4JetObject*> const& jets);

}


#endif
