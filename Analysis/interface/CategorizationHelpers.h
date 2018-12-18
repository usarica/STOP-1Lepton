#ifndef CATEGORIZATIONHELPERS_H
#define CATEGORIZATIONHELPERS_H

#include <vector>
#include "TString.h"
#include "SimpleEntry.h"


namespace CategorizationHelpers{
  enum Category{
    Inclusive,
    nCategories
  };
  enum CategorizationScheme{
    Moriond19_Analysis,
    Moriond19_ttZCR,
    nCatSchemes
  };

  extern CategorizationScheme globalCategorizationScheme;
  void setGlobalCategorizationScheme(CategorizationHelpers::CategorizationScheme scheme);
  std::vector<CategorizationHelpers::Category> getAllowedCategories(CategorizationHelpers::CategorizationScheme scheme);
  bool testCategoryAgainstGlobalScheme(CategorizationHelpers::Category theCategory);

  TString getCategoryName(CategorizationHelpers::Category category);
  TString getCategoryLabel(CategorizationHelpers::Category category);

  CategorizationHelpers::Category getCategory(SimpleEntry const& sel_vars);

}


#endif
