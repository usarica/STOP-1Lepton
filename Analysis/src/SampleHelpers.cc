#include "SampleHelpers.h"


namespace SampleHelpers{
  const DatasetInfoExtractor datasetInfoExtractor=setupDatasetInfoExtractor();
}

DatasetInfoExtractor SampleHelpers::setupDatasetInfoExtractor(){
  DatasetInfoExtractor res;
  res.loadFromFile(TString(CMSTASCOREPKGPATH + "Tools/datasetinfo/scale1fbs.txt").Data());
  return res;
}
float SampleHelpers::getDatasetXsec(const TString& strsample){
  return 0.;
}
