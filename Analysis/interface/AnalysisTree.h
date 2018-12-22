#ifndef ANALYSISTREE_H
#define ANALYSISTREE_H

#include "Samples.h"
#include "BaseTree.h"
#include "GoodEventFilter.h"


// Forward declarations
class AnalysisSet;


class AnalysisTree : public BaseTree{
protected:
  AnalysisSet* associatedSet;
  RunNumber_t* RunNumberRef;
  Lumisection_t* LumisectionRef;
  EventNumber_t* EventNumberRef;

  void autoBookBranches();

public:
  const bool isMC;

  AnalysisTree(TString strsample, bool isMC_, const TString treename=CMS4_EVENTS_TREE_NAME);
  ~AnalysisTree(){}

  void setAssociatedSet(AnalysisSet* inSet){ associatedSet = inSet; }
  AnalysisSet* getAssociatedSet(){ return associatedSet; }
  AnalysisSet const* getAssociatedSet() const{ return associatedSet; }

  bool isValidEvent() const;

  // Static functions
  static TString constructSampleIdentifier(TString strsample);

};


#endif
