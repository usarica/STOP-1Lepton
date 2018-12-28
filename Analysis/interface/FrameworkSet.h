#ifndef FRAMEWORKSET_H
#define FRAMEWORKSET_H

#include <iostream>
#include <vector>
#include "FrameworkTree.h"


class FrameworkSet{
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
  std::vector<FrameworkTree*> treeList;
  std::unordered_map<FrameworkTree*, float> permanentWeights;
  std::unordered_map<TString, std::vector<std::string>> sampleHeaders;

  void addSampleHeader(FrameworkTree* tree);

public:
  FrameworkSet(FrameworkOptionParser const& opts, const TString treename);
  ~FrameworkSet();

  bool addFrameworkTree(FrameworkOptionParser const& opts, const TString& fname, const TString treename);
  bool addFrameworkTreeList(FrameworkOptionParser const& opts, const TString treename);
  bool dissociateFrameworkTree(FrameworkTree*& tree);
  bool associateFrameworkTree(FrameworkTree*& tree);

  std::vector<FrameworkTree*> getFrameworkTree(TString sampleid) const;
  const std::vector<FrameworkTree*>& getFrameworkTreeList() const;
  std::vector<FrameworkTree*>& getFrameworkTreeList();

  FrameworkTree* getSelectedEvent(const int evid);
  FrameworkTree* getFailedEvent(const int evid);

  std::vector<std::string> getSampleHeader(FrameworkTree* tree);

};


#endif
