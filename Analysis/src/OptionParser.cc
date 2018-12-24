#include "OptionParser.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"
#include "Samples.h"
#include "SampleHelpers.h"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


OptionParser::OptionParser(int argc, char** argv) :
  indir("./"),
  outdir("./"),
  coutput("tmp.root"),
  maxEvents(-1)
{
  for (int a=0; a<argc; a++){
    string tmpArg(argv[a]);
    rawOptions.push_back(tmpArg);
  }
  analyze();
}

void OptionParser::analyze(){
  bool hasInvalidOption=false;
  bool redefinedOutputFile=false;
  char rawdelimiter = '=';
  for (unsigned int opt=1; opt<rawOptions.size(); opt++){
    string wish, value;
    splitOption(rawOptions.at(opt), wish, value, rawdelimiter);
    interpretOption(wish, value);
    if (wish=="outfile") redefinedOutputFile=true;
  }

  if (filename.empty()){ MELAerr << "You have to specify the input files." << endl; if (!hasInvalidOption) hasInvalidOption=true; }
  else{
    for (unsigned int f=0; f<filename.size(); f++){
      if (filename.at(f).find(".root")!=string::npos){
        MELAerr << "Inconsistent file name " << filename.at(f) << "!" << endl;
        if (!hasInvalidOption) hasInvalidOption=true;
      }
    }
  }

  // Check for any invalid options and print an error
  //

  // Warnings-only
  if (!redefinedOutputFile) MELAout << "WARNING: No output file specified. Defaulting to " << coutput << "." << endl;

  // Print help if needed and abort at this point, nowhere later
  if (hasInvalidOption) printOptionsHelp();

  // Append extra "/" if they do not exist.
  unsigned int tlen=(unsigned int) indir.length();
  if (tlen>1 && indir[tlen-1]!='/') indir.append("/");
  tlen=(unsigned int) outdir.length();
  if (tlen>1 && outdir[tlen-1]!='/') outdir.append("/");

}

void OptionParser::interpretOption(const std::string& wish, std::string value){
  if (wish.empty()){
    if (value.find(".lhe")!=string::npos || value.find(".root")!=string::npos) filename.push_back(value);
    else if (value.find("help")!= string::npos) printOptionsHelp();
    else MELAerr << "Unknown unspecified argument: " << value << endl;
  }

  else if (wish=="indir") indir = value;
  else if (wish=="outdir") outdir = value;
  else if (wish=="outfile") coutput = value;

  else if (wish=="maxevents" || wish=="maxEvents") maxEvents = (int) atoi(value.c_str());

  else MELAerr << "Unknown specified argument: " << value << " with specifier " << wish << endl;
}

void OptionParser::printOptionsHelp(){
  MELAout << endl;
  MELAout << "The options implemented for the LHEAnalyzer (format: specifier=value):\n\n";

  MELAout << "- No option specifier: Input files with extension .lhe or .root. Multiple input files can be passed as different arguments.\n\n";
  MELAout << "- indir: Location of input files. Default=\"./\"\n\n";
  MELAout << "- outdir: Location of the output file. Default=\"./\"\n\n";
  MELAout << "- outfile: Output file name. Default=\"tmp.root\"\n\n";
  MELAout << "- maxevents / maxEvents: Maximum number of events to process. Default=-1 (all events)\n\n";

  MELAout << endl;
  assert(0);
}
