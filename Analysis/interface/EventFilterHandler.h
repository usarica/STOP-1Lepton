#ifndef EVENTFILTERHANDLER_H
#define EVENTFILTERHANDLER_H

#include <vector>
#include <unordered_map>
#include "IvyBase.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"
#include "SystematicVariations.h"


class EventFilterHandler : public IvyBase{
protected:
  bool passEventFilterFlag;
  std::unordered_map<TString, bool> product_HLTpaths;

  void clear();

public:
  // Constructors
  EventFilterHandler();

  // Destructors
  ~EventFilterHandler(){ clear(); }

  bool constructFilter();
  std::unordered_map<TString, bool> const& getHLTProduct() const{ return product_HLTpaths; }
  bool passEventFilters() const{ return passEventFilterFlag; }

  void bookBranches(BaseTree* tree);
  static std::vector<TString> getMETFilterFlags(FrameworkTree* fwktree);
  static std::unordered_map<TString, std::vector<TString>> getHLTPaths(FrameworkTree* fwktree);

  // Special event filters for various specific issues
  bool test2017_2018HTFilter(
    std::vector<AK4JetObject*> const* ak4jets,
    SystematicsHelpers::SystematicVariationTypes syst
  ) const;
  bool test2018HEMFilter(
    std::vector<ElectronObject*> const* electrons,
    std::vector<PhotonObject*> const* photons,
    std::vector<AK4JetObject*> const* ak4jets,
    std::vector<AK8JetObject*> const* ak8jets,
    SystematicsHelpers::SystematicVariationTypes syst
  ) const;

};


#endif
