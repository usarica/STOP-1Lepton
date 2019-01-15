#ifndef EVENTFILTERHANDLER_H
#define EVENTFILTERHANDLER_H

#include <vector>
#include <unordered_map>
#include "IvyBase.h"


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

};


#endif
