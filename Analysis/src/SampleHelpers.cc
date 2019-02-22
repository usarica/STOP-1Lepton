#include <cassert>
#include "SampleHelpers.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


namespace SampleHelpers{
  const DatasetInfoExtractor datasetInfoExtractor=setupDatasetInfoExtractor();
  FrameworkOptionParser const* prevOpts = nullptr;
  FrameworkOptionParser const* currentOpts = nullptr;

  // These functions are hidden from the user!
  void setDataPeriod(TString s);
  void setDataVersion(TString s);
  void setInputDirectory(TString s);
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

void SampleHelpers::setDataPeriod(TString s){
  theDataPeriod = s;
  theDataYear = -1;
  if (theDataPeriod.Contains("2016")) theDataYear = 2016;
  else if (theDataPeriod.Contains("2017")) theDataYear = 2017;
  else if (theDataPeriod.Contains("2018")) theDataYear = 2018;
  else assert(0);
}
void SampleHelpers::setDataVersion(TString s){
  s.ToLower();
  if (s.Contains("80x")) theDataVersion = kCMSSW_8_0_X;
  else if (s.Contains("94x")) theDataVersion = kCMSSW_9_4_X;
  else if (s.Contains("10x")) theDataVersion = kCMSSW_10_X;
  else assert(0);
}
void SampleHelpers::setInputDirectory(TString s){ theInputDirectory=s; }
void SampleHelpers::setupUsingOptions(FrameworkOptionParser const& opts, bool doForce){
  if (currentOpts != &opts || doForce){
    setDataPeriod(opts.dataPeriod().c_str());
    setDataVersion(opts.dataVersion().c_str());
    setInputDirectory(opts.inputDir().c_str());

    // Refresh pointers to the previous and current options
    prevOpts = currentOpts;
    currentOpts = &opts;
  }
}
FrameworkOptionParser const* SampleHelpers::getCurrentOptions(){ return currentOpts; }
void SampleHelpers::revertToPreviousOptions(){
  if (!prevOpts) return;
  setupUsingOptions(*prevOpts);
}
