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
  theCondorSite(""),
  condorOutdir(""),
  sample(""),
  sampletag(""),
  outputName("tmp.root"),
  theDataPeriod(""),
  theDataVersion(""),
  maxEvents(-1),
  recordEveryN(-1),
  isMCflag(true),
  isFastSimflag(false),
  exceptionalCases()
{
  if (argc>0){
    MELAout << "Executing " << argv[0] << " with " << (argc>1 ? "options" : "no options.");
    if (argc>1) rawOptions.reserve(argc-1);
    for (int a=1; a<argc; a++){
      rawOptions.emplace_back(argv[a]);
      MELAout << ' ' << argv[a];
    }
    MELAout << endl;
  }
  analyze();
}
FrameworkOptionParser::FrameworkOptionParser(std::string opts) :
  indir("./"),
  outdir("./"),
  theCondorSite(""),
  condorOutdir(""),
  sample(""),
  sampletag(""),
  outputName("tmp.root"),
  theDataPeriod(""),
  theDataVersion(""),
  maxEvents(-1),
  recordEveryN(-1),
  isMCflag(true),
  isFastSimflag(false),
  exceptionalCases()
{
  splitOptionRecursive(opts, rawOptions, ' ');
  analyze();
}

void FrameworkOptionParser::analyze(){
  bool hasInvalidOption=false;
  bool redefinedOutputFile=false;
  const char rawdelimiter = '=';
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
  if (theDataPeriod==""){ MELAerr << "You have to specify the data set period." << endl; hasInvalidOption |= true; }
  else{
    if (this->isData()){
      int try_year=-1;
      bool has_exception = false;
      try{ try_year = stoi(theDataPeriod.c_str()); }
      catch (std::invalid_argument& e){ has_exception = true; }
      if (!has_exception && try_year>0){
        if (theDataPeriod == Form("%i", try_year)){ // Check if the data period string contains just the year
          MELAout << "The data period " << theDataPeriod << " contains only the year of the data file. Searching the file name for the data period..." << endl;
          const char test_chars[]="ABCDEFGHIJK";
          const unsigned int n_test_chars = strlen(test_chars);
          for (unsigned int ic=0; ic<n_test_chars; ic++){
            string test_data_period = Form("%i%c", try_year, test_chars[ic]);
            if (sample.find(test_data_period.c_str())!=string::npos){
              MELAout << "\t- The data period " << test_data_period << " is found." << endl;
              theDataPeriod = test_data_period;
              break;
            }
          }
        }
      }
    }
    if (theDataVersion==""){
      MELAout << "You have to specify the data set version, but attempting to find out from the data period." << endl;
      if (theDataPeriod.find("2018")!=std::string::npos) theDataVersion="10x";
      else if (theDataPeriod.find("2017")!=std::string::npos) theDataVersion="94x";
      else if (theDataPeriod.find("2016")!=std::string::npos){
        if (sample.find("2018")!=std::string::npos || sample.find("2019")!=std::string::npos || sample.find("94X")!=std::string::npos || sample.find("94x")!=std::string::npos) theDataVersion="94x"; // 94X samples are made in 2018 and 2019
        else theDataVersion="80x";
      }
      else{
        MELAerr << "\t- Could not determine the data set version." << endl;
        hasInvalidOption |= true;
      }
    }
  }

  // Check Condor specifications
  if ((theCondorSite!="" && condorOutdir=="") || (theCondorSite=="" && condorOutdir!="")){
    MELAerr << "Either the transfer site or the target directory is not defined for HTCondor." << endl;
    hasInvalidOption |= true;
  }
  if (theCondorSite!="" && condorOutdir!="" && outdir.find('/')==0){
    MELAerr << "Local output directory has to be relative when a remote transfer site is specified." << endl;
    hasInvalidOption |= true;
  }

  // Check for any invalid options and print an error
  if (isData() && isFastSim()){ MELAerr << "FastSim option is only usable in the MC!" << endl; isFastSimflag=false; }

  // Warnings-only
  if (!redefinedOutputFile) MELAout << "WARNING: No output file specified. Defaulting to " << outputName << "." << endl;

  if (theInputFileNames.empty()){
    TString dirname = SampleHelpers::getDatasetDirectoryName(*this);
    MELAout << "Attempting to find the input file names in " << dirname << "..." << endl;
    std::vector<TString> fnamelist = SampleHelpers::lsdir(dirname);
    if (!fnamelist.empty()){
      theInputFileNames.reserve(fnamelist.size());
      for (TString const& fname:fnamelist){ if (fname.Contains(".root")) theInputFileNames.emplace_back(fname.Data()); }
    }
    if (theInputFileNames.empty()){
      MELAerr << "\t- No ROOT files found!" << endl;
      hasInvalidOption |= true;
    }
    else MELAout << "\t- Found the following files: " << theInputFileNames << endl;
  }

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
  else if (wish=="inputfiles") splitOptionRecursive(value, theInputFileNames, ',', true);
  else if (wish=="outdir") outdir = value;
  else if (wish=="condorsite") theCondorSite = value;
  else if (wish=="condoroutdir") condorOutdir = value;
  else if (wish=="sample") sample = value;
  else if (wish=="sampletag" || wish=="tag") sampletag = value;
  else if (wish=="outfile") outputName = value;
  else if (wish=="period" || wish=="dataperiod" || wish=="year") theDataPeriod = value;
  else if (wish=="version" || wish=="dataversion" || wish=="release") theDataVersion = value;

  else if (wish=="maxevents") maxEvents = (int) atoi(value.c_str());
  else if (wish=="recordeveryn") recordEveryN = (int) atoi(value.c_str());

  else if (wish=="ismc" || wish=="isdata"){ HelperFunctions::castStringToValue(value, isMCflag); if (wish=="isdata") isMCflag = !isMCflag; }
  else if (wish=="isfastsim") HelperFunctions::castStringToValue(value, isFastSimflag);

  else if (wish=="specialpdf_nnpdf30_nlo_nf_4_pdfas_madgraph_1000offset_powhegstyle_case1") HelperFunctions::castStringToValue(value, exceptionalCases.specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1);
  else if (wish=="specialpdf_nnpdf31_nnlo_as_0118_nf_4") HelperFunctions::castStringToValue(value, exceptionalCases.specialPDF_NNPDF31_NNLO_as_0118_nf_4);
  else if (wish=="specialpdf_nnpdf31_nnlo_as_0118_madgraph_1000offset_case1") HelperFunctions::castStringToValue(value, exceptionalCases.specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1);

  else MELAerr << "Unknown specified argument: " << value << " with specifier " << wish << endl;
}

bool FrameworkOptionParser::findTagFromDatasetFile(){
  if (sample=="") return false;
  MELAout << "FrameworkOptionParser::findTagFromDatasetFile: Attempting to find the latest tag for sample " << sample << "..." << endl;
  std::unordered_map<std::string, DatasetInfoExtractor::datasetInfo> const& dsmap = SampleHelpers::datasetInfoExtractor.get_dslist();
  FrameworkTag latest_tag;
  if (this->isMC()){ // Search the DS info file for the MC
    for (std::unordered_map<std::string, DatasetInfoExtractor::datasetInfo>::const_iterator it=dsmap.cbegin(); it!=dsmap.cend(); it++){
      std::string strkey = it->first;
      if (strkey.find(sample)!=string::npos){
        replaceString<std::string, const char*>(strkey, sample.c_str(), "");
        FrameworkTag tmp_tag(strkey);
        if (tmp_tag>=latest_tag) latest_tag=tmp_tag;
      }
    }
  }
  else{ // Search the available folders for the data
    TString partial_sample = SampleHelpers::getDatasetDirectoryName(sample, "");
    //MELAout << "Partial sample dir: " << partial_sample << endl;
    std::vector<TString> dirlist = SampleHelpers::lsdir(indir);
    for (auto const& dir:dirlist){
      TString tmptag = dir;
      if (tmptag.Contains(partial_sample)) replaceString<TString, const char*>(tmptag, (partial_sample+"_").Data(), "");
      else continue;
      //MELAout << "Checking dir: " << dir << ", tag: " << tmptag << endl;
      FrameworkTag tmp_tag(tmptag.Data());
      if (tmp_tag>=latest_tag) latest_tag=tmp_tag;
    }
  }
  if (latest_tag != FrameworkTag()){
    sampletag=latest_tag.getRawTag();
    MELAout << "FrameworkOptionParser::findTagFromDatasetFile: Tag " << sampletag << " is found!" << endl;
    return true;
  }
  else return false;
}

void FrameworkOptionParser::printOptionsHelp(){
  MELAout << endl;
  MELAout << "The options implemented for the FrameworkOptionParser (format: specifier=value):\n\n";

  MELAout << "- indir: Location of input files. Default=\"./\"\n\n";
  MELAout << "- inputfiles: Comma-separated list of input file names. Default=\"\", corresponding to all files in the input directory.\n\n";
  MELAout << "- condorsite: Condor transsfer site. Default=\"\"\n\n";
  MELAout << "- condoroutdir: Transfer directory at the output site. Default=\"\"\n\n";
  MELAout << "- outdir: Local output directory. Default=\"./\"\n\t- When a remote transfer site is specified, this directory has to be relative.\n\n";
  MELAout << "- outfile: Output file name. Default=\"tmp.root\"\n\n";
  MELAout << "- sample: Sample name. Example: sample=/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM. Default=\"\"\n\n";
  MELAout << "- sampletag/tag: Sample tag. Example: sampletag=CMS4_V09-04-19. Default=\"\"\n\n";
  MELAout << "- period/dataperiod/year: The data period (2016, 2017, 2018 etc.). Default=\"\"\n\n";
  MELAout << "- version/dataversion/release: The data version (80x, 94x, 10x etc.). Default depends on the data period.\n\n";
  MELAout << "- maxevents: Maximum number of events to process. Default=-1 (all events)\n\n";
  MELAout << "- recordeveryn: Maximum number of events to loop before recording results to the base tree. Default=-1 (all events in the sub-tree)\n\n";
  MELAout << "- ismc/isdata: Specify whether the sample is from simulation or real data. Default=true (ismc=true)\n\n";
  MELAout << "- isfastsim: Specify whether the simulation sample is using FastSim. Default=false\n\n";
  MELAout << "- specialpdf_nnpdf30_nlo_nf_4_pdfas_madgraph_1000offset_powhegstyle_case1: Set the MC LHE weights flag for specialPDF_NNPDF30_nlo_nf_4_pdfas_Madgraph_1000offset_POWHEGStyle_Case1. Default=false\n\n";
  MELAout << "- specialpdf_nnpdf31_nnlo_as_0118_nf_4: Set the MC LHE weights flag for specialPDF_NNPDF31_NNLO_as_0118_nf_4. Default=false\n\n";
  MELAout << "- specialpdf_nnpdf31_nnlo_as_0118_madgraph_1000offset_case1: Set the MC LHE weights flag for specialPDF_NNPDF31_NNLO_as_0118_Madgraph_1000offset_Case1. Default=false\n\n";

  MELAout << endl;
  assert(0);
}
