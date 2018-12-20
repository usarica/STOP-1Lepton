#ifndef SAMPLEHELPERS_H
#define SAMPLEHELPERS_H

#include "Samples.h"
#include "DatasetInfoExtractor.h"


namespace SampleHelpers{
  extern const DatasetInfoExtractor datasetInfoExtractor;

  DatasetInfoExtractor setupDatasetInfoExtractor();
  float getDatasetXsec(const TString& strsample);
}


#endif
