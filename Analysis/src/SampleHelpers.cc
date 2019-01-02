#include <cassert>
#include "SampleHelpers.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;



namespace SampleHelpers{
  const DatasetInfoExtractor datasetInfoExtractor=setupDatasetInfoExtractor();
}

DatasetInfoExtractor SampleHelpers::setupDatasetInfoExtractor(){
  DatasetInfoExtractor res;
  res.loadFromFile(TString(CMSTASCOREPKGPATH + "Tools/datasetinfo/scale1fbs.txt").Data());
  return res;
}
float SampleHelpers::getDatasetXsec(const TString& strsample, const TString& strtag){ return datasetInfoExtractor.getXsecFromFile(strsample.Data(), strtag.Data()); }

TString SampleHelpers::getDatasetDirectoryName(std::string sname, std::string stag){
  if (sname.find('/')==0) sname.replace(0, 1, "");
  bool replaceAllSlashes=true;
  do{
    replaceAllSlashes = replaceString<std::string, const char*>(sname, "/", "_");
  }
  while (replaceAllSlashes);
  if (stag!="") return Form("%s_%s", sname.c_str(), stag.c_str());
  else return TString(sname.c_str());
}
TString SampleHelpers::getDatasetDirectoryName(TString sname, TString stag){ return SampleHelpers::getDatasetDirectoryName(std::string(sname.Data()), std::string(stag.Data())); }
TString SampleHelpers::getDatasetDirectoryName(FrameworkOptionParser const& opts){
  SampleHelpers::setupUsingOptions(opts);
  TString sname_tag = SampleHelpers::getDatasetDirectoryName(opts.sampleName(), opts.sampleTag());
  if (theInputDirectory=="" || (opts.sampleName().find('/')!=std::string::npos && theInputDirectory=="./")){
    MELAerr << "SampleHelpers::getDatasetDirectoryName: The main input directory is not set up!" << endl;
    assert(0);
  }
  return Form("%s/%s", theInputDirectory.Data(), sname_tag.Data());
}

void SampleHelpers::setupUsingOptions(FrameworkOptionParser const& opts){
  setDataPeriod(opts.dataPeriod().c_str());
  setDataVersion(opts.dataVersion().c_str());
  setInputDirectory(opts.inputDir().c_str());
}
