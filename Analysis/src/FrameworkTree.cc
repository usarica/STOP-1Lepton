#include "SampleHelpersCore.h"
#include "FrameworkTree.h"


FrameworkTree::FrameworkTree(FrameworkOptionParser const& opts, const TString fname, const TString treename) :
  BaseTree(SampleHelpers::getDatasetDirectoryName(opts)+'/'+fname, treename, "", ""),
  options(opts),
  tag(opts.sampleTag()),
  associatedSet(nullptr), RunNumberRef(nullptr), LumisectionRef(nullptr), EventNumberRef(nullptr)
{
  sampleIdentifier = FrameworkTree::constructSampleIdentifier();
  if (this->isValid()) autoBookBranches();
}

TString FrameworkTree::constructSampleIdentifier(){
  TString res="";
  std::vector<TString> splitstr; char delimiter='/';
  HelperFunctions::splitOptionRecursive(options.sampleName().c_str(), splitstr, delimiter);
  for (TString const& strtmp:splitstr){ if (strtmp!=""){ res = strtmp; break; } }
  return res;
}

void FrameworkTree::autoBookBranches(){
  if (!this->isMC()){
    this->bookBranch<RunNumber_t>("uint_eventMaker_evtrun_CMS3.obj", 0); this->getValRef("uint_eventMaker_evtrun_CMS3.obj", RunNumberRef);
    this->bookBranch<Lumisection_t>("uint_eventMaker_evtlumiBlock_CMS3.obj", 0); this->getValRef("uint_eventMaker_evtlumiBlock_CMS3.obj", LumisectionRef);
    this->bookBranch<EventNumber_t>("ull_eventMaker_evtevent_CMS3.obj", 0); this->getValRef("ull_eventMaker_evtevent_CMS3.obj", EventNumberRef);
  }
  else{
    // Calculation of xsec * BR is done outside
  }
}

bool FrameworkTree::isValidEvent() const{
  if (!BaseTree::isValidEvent()) return false;
  if (this->isMC()) return true;
  if (!RunNumberRef || !LumisectionRef) return false;
  return GoodEventFilter::testEvent(*RunNumberRef, *LumisectionRef);
}
