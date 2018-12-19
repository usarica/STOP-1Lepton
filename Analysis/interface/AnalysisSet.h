#ifndef ANALYSISSET_H
#define ANALYSISSET_H

#include <iostream>
#include <vector>
#include "AnalysisTree.h"


class AnalysisSet{
public:
  enum NormScheme{
    NormScheme_None,
    NormScheme_OneOverNgen,
    NormScheme_OneOverNgen_RenormBySumOneOverNgen,
    NormScheme_OneOverNgen_RelRenormToSumNgen,
    NormScheme_NgenOverNgenWPU,
    NormScheme_XsecOnly,
    NormScheme_XsecOverNgen,
    NormScheme_XsecOverNgen_RenormBySumXsecOverNgen,
    NormScheme_XsecOverNgen_RelRenormToSumNgen
  };

protected:
  std::vector<AnalysisTree*> treeList;
  std::unordered_map<AnalysisTree*, float> permanentWeights;

public:
  AnalysisSet(const TString& strname, const bool isMC, const TString treename);
  AnalysisSet(const std::vector<TString>& strlist, const bool isMC, const TString treename);
  ~AnalysisSet();

  bool addAnalysisTree(const TString& strname, const bool isMC, const TString treename);
  bool addAnalysisTreeList(const std::vector<TString>& strlist, const bool isMC, const TString treename);
  bool dissociateAnalysisTree(AnalysisTree*& tree);
  bool associateAnalysisTree(AnalysisTree*& tree);

  AnalysisTree* getAnalysisTree(TString sampleid) const;
  const std::vector<AnalysisTree*>& getAnalysisTreeList() const;
  std::vector<AnalysisTree*>& getAnalysisTreeList();

  AnalysisTree* getSelectedEvent(const int evid);
  AnalysisTree* getFailedEvent(const int evid);

};


#endif
