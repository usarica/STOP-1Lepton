#include <cmath>
#include <cassert>
#include "CategorizationHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


namespace CategorizationHelpers{
  CategorizationScheme globalCategorizationScheme=Moriond19_Analysis;
}

using namespace CategorizationHelpers;


TString CategorizationHelpers::getCategoryName(CategorizationHelpers::Category category){
  switch(category){
  case Inclusive:
    return "Inclusive";
  default:
    return "";
  }
}

TString CategorizationHelpers::getCategoryLabel(CategorizationHelpers::Category category){
  switch (category){
  case Inclusive:
    return "Inclusive";
  default:
    return "";
  }
}

void CategorizationHelpers::setGlobalCategorizationScheme(CategorizationHelpers::CategorizationScheme scheme){ globalCategorizationScheme=scheme; }

CategorizationHelpers::Category CategorizationHelpers::getCategory(SimpleEntry const& sel_vars){
  if (globalCategorizationScheme==Moriond19_Analysis){

  }
  else if (globalCategorizationScheme==Moriond19_ttZCR){

  }
  return Inclusive;
}
std::vector<CategorizationHelpers::Category> CategorizationHelpers::getAllowedCategories(CategorizationHelpers::CategorizationScheme scheme){
  std::vector<CategorizationHelpers::Category> res;
  res.push_back(Inclusive);
  /*
  switch (scheme){
  case Moriond19_ttZCR:
  case Moriond19_Analysis:
    res.push_back(JJVBFTagged);
    res.push_back(Untagged);
  default:
    break;
  }
  */
  return res;
}
bool CategorizationHelpers::testCategoryAgainstGlobalScheme(CategorizationHelpers::Category theCategory){
  std::vector<Category> cats=getAllowedCategories(globalCategorizationScheme);
  for (Category const& cat:cats){ if (cat==theCategory) return true; }
  return false;
}
