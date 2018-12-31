#ifndef FRAMEWORKTREE_H
#define FRAMEWORKTREE_H

#include "SampleHelpers.h"
#include "BaseEDMInputTree.h"
#include "GoodEventFilter.h"
#include "FrameworkOptionParser.h"
#include "FrameworkTag.h"


// Forward declarations
class FrameworkSet;


class FrameworkTree : public BaseEDMInputTree{
protected:
  FrameworkOptionParser options;
  FrameworkTag tag;

  FrameworkSet* associatedSet;
  RunNumber_t* RunNumberRef;
  Lumisection_t* LumisectionRef;
  EventNumber_t* EventNumberRef;

  void autoBookBranches();

public:
  FrameworkTree(FrameworkOptionParser const& opts, const TString fname, const TString treename=CMS4_EVENTS_TREE_NAME);
  ~FrameworkTree(){}

  void setAssociatedSet(FrameworkSet* inSet){ associatedSet = inSet; }
  FrameworkSet* getAssociatedSet(){ return associatedSet; }
  FrameworkSet const* getAssociatedSet() const{ return associatedSet; }

  bool isMC() const{ return options.isMC(); }
  bool isValidEvent() const;
  FrameworkOptionParser const& getOptions() const{ return options; }
  FrameworkTag const& getTag() const{ return tag; }

  TString constructSampleIdentifier();

};


#endif
