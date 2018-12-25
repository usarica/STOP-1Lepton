#ifndef ANALYSISTREE_H
#define ANALYSISTREE_H

#include "SampleHelpers.h"
#include "BaseTree.h"
#include "GoodEventFilter.h"
#include "FrameworkOptionParser.h"
#include "FrameworkTag.h"


// Forward declarations
class AnalysisSet;


class AnalysisTree : public BaseTree{
protected:
  FrameworkOptionParser options;
  FrameworkTag tag;

  AnalysisSet* associatedSet;
  RunNumber_t* RunNumberRef;
  Lumisection_t* LumisectionRef;
  EventNumber_t* EventNumberRef;

  void autoBookBranches();

public:
  AnalysisTree(FrameworkOptionParser const& opts, const TString fname, const TString treename=CMS4_EVENTS_TREE_NAME);
  ~AnalysisTree(){}

  void setAssociatedSet(AnalysisSet* inSet){ associatedSet = inSet; }
  AnalysisSet* getAssociatedSet(){ return associatedSet; }
  AnalysisSet const* getAssociatedSet() const{ return associatedSet; }

  bool isMC() const{ return options.isMC(); }
  bool isValidEvent() const;
  FrameworkOptionParser const& getOptions() const{ return options; }
  FrameworkTag const& getTag() const{ return tag; }

  TString constructSampleIdentifier();

};


#endif
