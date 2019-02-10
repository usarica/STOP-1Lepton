#ifndef FRAMEWORKJOBSPLITTER_H
#define FRAMEWORKJOBSPLITTER_H

#include <string>
#include <vector>


class FrameworkJobSplitter{
protected:
  std::string rawopts;
  int nFilesPerJob;
  std::vector<std::string> theSplitOptions;

  void splitOptions();

public:
  FrameworkJobSplitter(std::string opts, int nFilesPerJob_);
  ~FrameworkJobSplitter(){}

  std::vector<std::string> const& getSplitOptions(){ return theSplitOptions; }

};

#endif
