#ifndef OPTIONPARSER_H
#define OPTIONPARSER_H

#include <string>
#include <vector>


class OptionParser{
protected:
  std::vector<std::string> rawOptions;
  std::vector<std::string> filename;

  std::string indir;
  std::string outdir;
  std::string coutput;
  int maxEvents;

public:
  OptionParser(int argc, char** argv);
  ~OptionParser(){}

  void analyze();
  void interpretOption(const std::string& wish, std::string value);
  void printOptionsHelp();

  std::string const& inputDir() const{ return indir; }
  std::string const& outputDir() const{ return outdir; }
  std::string const& outputFilename() const{ return coutput; }
  std::vector<std::string> const& inputfiles() const{ return filename; }
  int const& maxEventsToProcess() const{ return maxEvents; }

};

#endif
