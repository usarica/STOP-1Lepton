#include <cassert>
#include "Samples.h"

namespace SampleHelpers{
  int theDataYear=2016;
  DataVersion theDataVersion=kCMSSW_8_0_X;
  TString theDataPeriod="2016"; // Initialize the extern here to 2016
  TString theInputDirectory=""; // Initialize the extern here to empty string
}

void SampleHelpers::setDataPeriod(TString s){
  theDataPeriod = s;
  theDataYear = -1;
  if (theDataPeriod.Contains("2016")) theDataYear = 2016;
  else if (theDataPeriod.Contains("2017")) theDataYear = 2017;
  else if (theDataPeriod.Contains("2018")) theDataYear = 2018;
  else assert(0);
}
void SampleHelpers::setDataVersion(TString s){
  s.ToLower();
  if (s.Contains("80x")) theDataVersion = kCMSSW_8_0_X;
  else if (s.Contains("94x")) theDataVersion = kCMSSW_9_4_X;
  else if (s.Contains("10x")) theDataVersion = kCMSSW_10_X;
  else assert(0);
}
void SampleHelpers::setInputDirectory(TString s){ theInputDirectory=s; }
