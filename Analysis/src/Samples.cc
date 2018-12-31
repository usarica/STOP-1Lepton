#include <cassert>
#include "Samples.h"

namespace SampleHelpers{
  int theDataYear=2016;
  TString theDataPeriod="2016"; // Initialize the extern here to 2016
  TString theInputDirectory=""; // Initialize the extern here to empty string
}

void SampleHelpers::setDataPeriod(const TString s){
  theDataPeriod = s;
  theDataYear = -1;
  if (theDataPeriod.Contains("2016")) theDataYear = 2016;
  else if (theDataPeriod.Contains("2017")) theDataYear = 2017;
  else if (theDataPeriod.Contains("2018")) theDataYear = 2018;
  else assert(0);
}
void SampleHelpers::setInputDirectory(const TString s){ theInputDirectory=s; }
