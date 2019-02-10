#include "FrameworkJobSplitter.h"
#include "MELAStreamHelpers.hh"


int main(int argc, char** argv){
  using namespace std;
  using namespace MELAStreamHelpers;

  if (argc!=3){
    MELAerr << "SplitFrameworkJob: You need the executable name (ie. SplitFrameworkJob), the option string in double quotes, and the number of splittings, which cannot be 0." << endl;
    return 1;
  }

  FrameworkJobSplitter splitter(argv[1], atoi(argv[2]));
  for (auto const& stropt:splitter.getSplitOptions()) MELAout << "\"" << stropt << "\"" << endl;
  return 0;
}