#ifndef SAMPLEHELPERS_H
#define SAMPLEHELPERS_H

#include "Samples.h"
#include "DatasetInfoExtractor.h"
#include "FrameworkOptionParser.h"


namespace SampleHelpers{
  extern const DatasetInfoExtractor datasetInfoExtractor;

  DatasetInfoExtractor setupDatasetInfoExtractor();
  float getDatasetXsec(const TString& strsample, const TString& strtag);

  TString getDatasetDirectoryName(FrameworkOptionParser const& opts);
  void setupUsingOptions(FrameworkOptionParser const& opts);

}


#endif
