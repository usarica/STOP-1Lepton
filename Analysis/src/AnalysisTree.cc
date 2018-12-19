#include "SampleHelpersCore.h"
#include "AnalysisTree.h"


AnalysisTree::AnalysisTree(TString strsample, bool isMC_, const TString treename) :
  BaseTree(strsample, treename, "", ""),
  associatedSet(nullptr), RunNumberRef(nullptr), LumisectionRef(nullptr),
  isMC(isMC_),
  sampleIdentifier(AnalysisTree::constructSampleIdentifier(strsample))
{
  if (this->isValid()) autoBookBranches();
}

TString AnalysisTree::constructSampleIdentifier(TString strsample){
  TString res="";
  std::vector<TString> splitstr; char delimiter='/';
  HelperFunctions::splitOptionRecursive(strsample, splitstr, delimiter);
  for (std::vector<TString>::reverse_iterator rit=splitstr.rbegin(); rit!=splitstr.rend(); rit++){
    const TString& strtmp = *rit;
    if (strtmp!=""){ res = strtmp; break; }
  }
  return res;
}

void AnalysisTree::autoBookBranches(){
  if (!isMC){
    this->bookBranch<RunNumber_t>("uint_eventMaker_evtrun_CMS3.obj", 0); this->getValRef("uint_eventMaker_evtrun_CMS3.obj", RunNumberRef);
    this->bookBranch<Lumisection_t>("uint_eventMaker_evtlumiBlock_CMS3.obj", 0); this->getValRef("uint_eventMaker_evtlumiBlock_CMS3.obj", LumisectionRef);
  }
  else{
    // Calculation of xsec * BR is done outside
  }
}

bool AnalysisTree::isValidEvent() const{
  if (isMC) return true;
  if (!RunNumberRef || !LumisectionRef) return false;
  return GoodEventFilter::testEvent(*RunNumberRef, *LumisectionRef);
}
