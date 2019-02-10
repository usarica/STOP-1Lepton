#include <cassert>
#include "FrameworkJobSplitter.h"
#include "FrameworkOptionParser.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


FrameworkJobSplitter::FrameworkJobSplitter(std::string opts, int nFilesPerJob_) :
  rawopts(opts),
  nFilesPerJob(nFilesPerJob_)
{
  splitOptions();
}

void FrameworkJobSplitter::splitOptions(){
  if (nFilesPerJob==0) return;
  else if (nFilesPerJob<0){
    theSplitOptions.reserve(1);
    theSplitOptions.push_back(rawopts);
    return;
  }

  FrameworkOptionParser optparser(rawopts);
  std::string const& stroutput = optparser.outputFilename();
  std::vector<std::string> const& files = optparser.inputFileNames();
  std::vector<std::vector<std::string>> filegroups;
  {
    unsigned int ifile=0;
    for (auto const& file:files){
      if ((int) ifile==nFilesPerJob) ifile=0;
      if (ifile==0) filegroups.push_back(std::vector<std::string>());
      filegroups.back().push_back(file);
      ifile++;
    }
  }

  theSplitOptions.reserve(filegroups.size());
  for (unsigned int ifg=0; ifg<filegroups.size(); ifg++){
    auto const& fg = filegroups.at(ifg);
    std::string rawopts_copy = rawopts;
    std::string stroutput_copy = stroutput;
    HelperFunctions::replaceString<std::string, const std::string>(stroutput_copy, ".root", std::string("_")+std::to_string(ifg)+".root");
    HelperFunctions::replaceString<std::string, const std::string>(rawopts_copy, stroutput, stroutput_copy);
    std::string strappend = " inputfiles=";
    for (unsigned int ifile=0; ifile<fg.size(); ifile++){ strappend += fg.at(ifile); if (ifile!=fg.size()-1) strappend += ","; }
    theSplitOptions.push_back(rawopts_copy+strappend);
  }
}
