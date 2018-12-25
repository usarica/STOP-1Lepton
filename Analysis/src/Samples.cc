#include "Samples.h"

namespace SampleHelpers{
  TString theDataPeriod="2016"; // Initialize the extern here to 2016
  TString theInputDirectory=""; // Initialize the extern here to empty string
}

void SampleHelpers::setDataPeriod(const TString s){ theDataPeriod = s; }
void SampleHelpers::setInputDirectory(const TString s){ theInputDirectory=s; }
