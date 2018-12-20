#include "Samples.h"

namespace SampleHelpers{
  TString theDataPeriod="2016"; // Initialize the extern here to 2016
}

void SampleHelpers::setDataPeriod(const TString theDataPeriod_){ theDataPeriod = theDataPeriod_; }
