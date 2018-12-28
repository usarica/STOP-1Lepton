#ifndef FRAMEWORKOPTIONPARSER_H
#define FRAMEWORKOPTIONPARSER_H

#include <string>
#include <vector>
#include "LHEWeightHandler.h"


class FrameworkOptionParser{
protected:
  std::vector<std::string> rawOptions;

  std::string indir;
  std::string outdir;
  std::string sample;
  std::string sampletag;
  std::string outputName;
  std::string theDataPeriod;
  int maxEvents;
  bool isMCflag;

  LHEWeightHandler::ExceptionalCases exceptionalCases;

  bool findTagFromDatasetFile();

public:
  FrameworkOptionParser(int argc, char** argv);
  FrameworkOptionParser(std::string opts);
  ~FrameworkOptionParser(){}

  void analyze();
  void interpretOption(const std::string& wish, std::string value);
  void printOptionsHelp();

  std::string const& inputDir() const{ return indir; }
  std::string const& outputDir() const{ return outdir; }
  std::string const& outputFilename() const{ return outputName; }
  std::string const& sampleName() const{ return sample; }
  std::string const& sampleTag() const{ return sampletag; }
  std::string const& dataPeriod() const{ return theDataPeriod; }
  int const& maxEventsToProcess() const{ return maxEvents; }
  bool isMC() const{ return isMCflag; }
  bool isData() const{ return !this->isMC(); }
  LHEWeightHandler::ExceptionalCases const& getLHEExceptionalCases() const{ return exceptionalCases; }

};

#endif
