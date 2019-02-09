#ifndef FRAMEWORKOPTIONPARSER_H
#define FRAMEWORKOPTIONPARSER_H

#include <string>
#include <vector>
#include "LHEWeightHandler.h"


class FrameworkOptionParser{
protected:
  std::vector<std::string> rawOptions;

  std::vector<std::string> theInputFileNames;
  std::string indir;
  std::string outdir;
  std::string theCondorSite;
  std::string condorOutdir;
  std::string sample;
  std::string sampletag;
  std::string outputName;
  std::string theDataPeriod;
  std::string theDataVersion;
  int maxEvents;
  int recordEveryN;
  bool isMCflag;
  bool isFastSimflag;

  LHEWeightHandler::ExceptionalCases exceptionalCases;

  bool findTagFromDatasetFile();

public:
  FrameworkOptionParser(int argc, char** argv);
  FrameworkOptionParser(std::string opts);
  ~FrameworkOptionParser(){}

  void analyze();
  void interpretOption(const std::string& wish, std::string value);
  void printOptionsHelp();

  std::vector<std::string> const& inputFileNames() const{ return theInputFileNames; }
  std::string const& inputDir() const{ return indir; }
  std::string const& outputDir() const{ return outdir; }
  std::string const& condorSite() const{ return theCondorSite; }
  std::string const& condorOutputDir() const{ return condorOutdir; }
  std::string const& outputFilename() const{ return outputName; }
  std::string const& sampleName() const{ return sample; }
  std::string const& sampleTag() const{ return sampletag; }
  std::string const& dataPeriod() const{ return theDataPeriod; }
  std::string const& dataVersion() const{ return theDataVersion; }
  int const& maxEventsToProcess() const{ return maxEvents; }
  int const& recordEveryNEvents() const{ return recordEveryN; }
  bool isFastSim() const{ return isFastSimflag; }
  bool isMC() const{ return isMCflag; }
  bool isData() const{ return !this->isMC(); }
  LHEWeightHandler::ExceptionalCases const& getLHEExceptionalCases() const{ return exceptionalCases; }

};

#endif
