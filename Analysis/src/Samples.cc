#include <cassert>
#include <string>
#include <stdexcept>
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

bool SampleHelpers::testDataPeriodIsLikeData(){
  int try_year=-1;
  bool has_exception = false;
  try{ try_year = std::stoi(theDataPeriod.Data()); }
  catch (std::invalid_argument& e){ has_exception = true; }
  if (!has_exception && try_year>0){
    // Check if the data period string contains just the year
    if (theDataPeriod == Form("%i", try_year)) return false;
    else{
      const char test_chars[]="ABCDEFGHIJK";
      const unsigned int n_test_chars = strlen(test_chars);
      for (unsigned int ic=0; ic<n_test_chars; ic++){
        TString test_data_period = Form("%i%c", try_year, test_chars[ic]);
        if (theDataPeriod.Contains(test_data_period)) return true;
      }
      return false;
    }
  }
  else return false;
}

