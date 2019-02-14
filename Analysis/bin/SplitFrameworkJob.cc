#include <algorithm>
#include "FrameworkJobSplitter.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


int main(int argc, char** argv){
  using namespace std;
  using namespace HelperFunctions;
  using namespace MELAStreamHelpers;

  if (argc<3){
    MELAerr << "SplitFrameworkJob: You need the executable name (ie. SplitFrameworkJob), the option string in double quotes, and the number of splittings, which cannot be 0." << endl;
    MELAerr << "\t- You passed " << argc << " arguments, which were ";
    for (int ix=0; ix<argc; ix++) MELAerr << argv[ix] << ' ';
    MELAerr << endl;
    return 1;
  }

  int nfiles=0;
  std::string outfile="";
  {
    std::vector<std::string> rawOptions; rawOptions.reserve(argc-2);
    for (int a=2; a<argc; a++) rawOptions.emplace_back(argv[a]);
    const char rawdelimiter = '=';
    for (auto const& rawopt:rawOptions){
      std::string wish, value;
      splitOption(rawopt, wish, value, rawdelimiter);
      std::transform(wish.begin(), wish.end(), wish.begin(), ::tolower);
      if (wish=="nfiles") nfiles=stoi(value.c_str());
      else if (wish=="outfile") outfile=value;
    }
  }

  FrameworkJobSplitter splitter(argv[1], nfiles);
  // It is important that the file is opened right before the final printout to avoid printing other stuff.
  if (outfile!="") MELAout.open(outfile.c_str(), std::ios_base::app);
  for (auto const& stropt:splitter.getSplitOptions()) MELAout << "\"" << stropt << "\"" << endl;
  if (outfile!="") MELAout.close();
  return 0;
}