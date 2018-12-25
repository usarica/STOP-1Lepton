#include "FrameworkOptionParser.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"
#include "Samples.h"
#include "SampleHelpers.h"
#include "FrameworkTag.h"
#include <algorithm>


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


FrameworkOptionParser::FrameworkOptionParser(int argc, char** argv) :
  indir("./"),
  outdir("./"),
  sample(""),
  sampletag(""),
  outputName("tmp.root"),
  theDataPeriod(""),
  maxEvents(-1),
  isMCflag(true)
{
  if (argc>0) MELAout << "Executing " << argv[0] << " with " << (argc>1 ? "options" : "no options.");
  for (int a=1; a<argc; a++){
    string tmpArg(argv[a]);
    rawOptions.push_back(tmpArg);
    MELAout << ' ' << argv[a];
  }
  if (argc>0) MELAout << endl;
  analyze();
}
FrameworkOptionParser::FrameworkOptionParser(std::string opts) :
  indir("./"),
  outdir("./"),
  sample(""),
  sampletag(""),
  outputName("tmp.root"),
  theDataPeriod(""),
  maxEvents(-1),
  isMCflag(true)
{
  splitOptionRecursive(opts, rawOptions, ' ');
  analyze();
}

void FrameworkOptionParser::analyze(){
  bool hasInvalidOption=false;
  bool redefinedOutputFile=false;
  char rawdelimiter = '=';
  for (auto const& rawopt:rawOptions){
    string wish, value;
    splitOption(rawopt, wish, value, rawdelimiter);
    std::transform(wish.begin(), wish.end(), wish.begin(), ::tolower);
    interpretOption(wish, value);
    if (wish=="outfile") redefinedOutputFile=true;
  }

  if (sample==""){ MELAerr << "You have to specify the input sample." << endl; hasInvalidOption |= true; }
  else if (sample.find(".root")!=string::npos){
    MELAerr << "Sample name " << sample << " cannot have the \\.root\\ exptension!" << endl;
    hasInvalidOption |= true;
  }
  if (sampletag==""){
    MELAerr << "You have to specify the input sample tag." << endl;
    if (sample.find('/')==0){ // Attempt to find the latest tag for this sample
      hasInvalidOption |= !findTagFromDatasetFile();
    }
    else hasInvalidOption |= true;
  }
  if (theDataPeriod==""){ MELAerr << "You have to specify the dataset period." << endl; hasInvalidOption |= true; }

  // Check for any invalid options and print an error
  //

  // Warnings-only
  if (!redefinedOutputFile) MELAout << "WARNING: No output file specified. Defaulting to " << outputName << "." << endl;

  // Print help if needed and abort at this point, nowhere later
  if (hasInvalidOption) printOptionsHelp();

  // Append extra "/" if they do not exist.
  unsigned int tlen=(unsigned int) indir.length();
  if (tlen>1 && indir[tlen-1]!='/') indir.append("/");
  tlen=(unsigned int) outdir.length();
  if (tlen>1 && outdir[tlen-1]!='/') outdir.append("/");

}

void FrameworkOptionParser::interpretOption(const std::string& wish, std::string value){
  if (wish.empty()){
    if (value.find("help")!=string::npos) printOptionsHelp();
    else{
      MELAerr << "FrameworkOptionParser::interpretOption: Option " << value << " is not recognized!" << endl;
      printOptionsHelp();
    }
  }

  else if (wish=="indir") indir = value;
  else if (wish=="outdir") outdir = value;
  else if (wish=="sample") sample = value;
  else if (wish=="sampletag" || wish=="tag") sampletag = value;
  else if (wish=="outfile") outputName = value;
  else if (wish=="period" || wish=="dataperiod" || wish=="year") theDataPeriod = value;

  else if (wish=="maxevents") maxEvents = (int) atoi(value.c_str());

  else if (wish=="ismc" || wish=="isdata"){ HelperFunctions::castStringToValue(value, isMCflag); if (wish=="isdata") isMCflag = !isMCflag; }

  else MELAerr << "Unknown specified argument: " << value << " with specifier " << wish << endl;
}

bool FrameworkOptionParser::findTagFromDatasetFile(){
  MELAout << "FrameworkOptionParser::findTagFromDatasetFile: Attempting to find the latest tag for sample " << sample << "..." << endl;
  std::unordered_map<std::string, DatasetInfoExtractor::datasetInfo> const& dsmap = SampleHelpers::datasetInfoExtractor.get_dslist();
  FrameworkTag latest_tag;
  for (std::unordered_map<std::string, DatasetInfoExtractor::datasetInfo>::const_iterator it=dsmap.cbegin(); it!=dsmap.cend(); it++){
    std::string strkey = it->first;
    if (strkey.find(sample)!=string::npos){
      replaceString<std::string, const char*>(strkey, sample.c_str(), "");
      FrameworkTag tmp_tag(strkey);
      if (tmp_tag>=latest_tag) latest_tag=tmp_tag;
    }
  }
  if (latest_tag != FrameworkTag()){ sampletag=latest_tag.getRawTag(); return true; }
  else return false;
}

void FrameworkOptionParser::printOptionsHelp(){
  MELAout << endl;
  MELAout << "The options implemented for the LHEAnalyzer (format: specifier=value):\n\n";

  MELAout << "- indir: Location of input files. Default=\"./\"\n\n";
  MELAout << "- outdir: Location of the output file. Default=\"./\"\n\n";
  MELAout << "- sample: Sample name. Example: sample=/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM. Default=\"\"\n\n";
  MELAout << "- sampletag/tag: Sample tag. Example: sampletag=CMS4_V09-04-19. Default=\"\"\n\n";
  MELAout << "- outfile: Output file name. Default=\"tmp.root\"\n\n";
  MELAout << "- period/dataperiod/year: The data period (2016, 2017, 2018 etc.). Default=\"\"\n\n";
  MELAout << "- maxevents: Maximum number of events to process. Default=-1 (all events)\n\n";
  MELAout << "- ismc/isdata: Specify whether the sample is from simulation or real data. Default=true (ismc=true)\n\n";

  MELAout << endl;
  assert(0);
}
