#include "FrameworkVariables.hh"
#include "SampleHelpersCore.h"
#include "FrameworkTree.h"


FrameworkTree::FrameworkTree(FrameworkOptionParser const& opts, const TString fname, const TString treename) :
  BaseEDMInputTree(SampleHelpers::getDatasetDirectoryName(opts)+'/'+fname, treename, "", ""),
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
    this->bookBranch<RunNumber_t>(_event_RunNumber_, 0); this->getValRef(_event_RunNumber_, RunNumberRef);
    this->bookBranch<Lumisection_t>(_event_Lumisection_, 0); this->getValRef(_event_Lumisection_, LumisectionRef);
    this->bookBranch<EventNumber_t>(_event_EventNumber_, 0); this->getValRef(_event_EventNumber_, EventNumberRef);
  }
  // Calculation of xsec * BR is done outside in this framework
}

bool FrameworkTree::isValidEvent() const{
  if (!BaseEDMInputTree::isValidEvent()) return false;
  if (this->isMC()) return true;
  if (!RunNumberRef || !LumisectionRef) return false;
  return GoodEventFilter::testEvent(*RunNumberRef, *LumisectionRef);
}
