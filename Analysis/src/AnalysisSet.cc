#include "SampleHelpersCore.h"
#include "AnalysisSet.h"
#include <algorithm>


using namespace std;


AnalysisSet::AnalysisSet(const TString& strname, const bool isMC, const TString treename){
  addAnalysisTree(strname, isMC, treename);
}
AnalysisSet::AnalysisSet(const std::vector<TString>& strlist, const bool isMC, const TString treename){
  addAnalysisTreeList(strlist, isMC, treename);
}
AnalysisSet::~AnalysisSet(){
  for (auto& tree:treeList) delete tree;
  treeList.clear();
}
bool AnalysisSet::addAnalysisTree(const TString& strname, const bool isMC, const TString treename){
  AnalysisTree* tree = new AnalysisTree(strname, isMC, treename);
  if (tree->isValid()) treeList.push_back(tree);
  else{ delete tree; tree=nullptr; }
  if (!tree) cerr << "AnalysisSet::addAnalysisTree(" << strname << ") is invalid!" << endl;
  else cout << "AnalysisSet::addAnalysisTree(" << strname << ") is successful!" << endl;
  if (tree) tree->setAssociatedSet(this);
  return (tree!=nullptr);
}
bool AnalysisSet::addAnalysisTreeList(const std::vector<TString>& strlist, const bool isMC, const TString treename){
  bool res=true;
  for (auto const& s:strlist) res &= addAnalysisTree(s, isMC, treename);
  return res;
}
bool AnalysisSet::dissociateAnalysisTree(AnalysisTree*& tree){
  if (!tree) return false;
  auto it = std::find(treeList.begin(), treeList.end(), tree);
  if (it!=treeList.end()){ treeList.erase(it); tree->setAssociatedSet(nullptr); return true; }
  return false;
}
bool AnalysisSet::associateAnalysisTree(AnalysisTree*& tree){
  if (!tree) return false;
  AnalysisSet* theSet = tree->getAssociatedSet();
  if (theSet==this) return true;
  else if (theSet){
    theSet->dissociateAnalysisTree(tree);
    treeList.push_back(tree);
    tree->setAssociatedSet(this);
    return true;
  }
  return false;
}


AnalysisTree* AnalysisSet::getAnalysisTree(TString sampleid) const{
  for (auto tree:treeList){ if (tree->sampleIdentifier==sampleid) return tree; }
  AnalysisTree* res=nullptr;
  return res;
}
const std::vector<AnalysisTree*>& AnalysisSet::getAnalysisTreeList() const{ return treeList; }
std::vector<AnalysisTree*>& AnalysisSet::getAnalysisTreeList(){ return treeList; }

AnalysisTree* AnalysisSet::getSelectedEvent(const int evid){
  int ev = evid;
  AnalysisTree* isfound=nullptr;
  for (auto& tree:treeList){
    int nevts = tree->getSelectedNEvents();
    if (ev<nevts){ tree->getSelectedEvent(ev); isfound = tree; break; }
    else ev -= nevts;
  }
  return isfound;
}
AnalysisTree* AnalysisSet::getFailedEvent(const int evid){
  int ev = evid;
  AnalysisTree* isfound=nullptr;
  for (auto& tree:treeList){
    int nevts = tree->getFailedNEvents();
    if (ev<nevts){ tree->getFailedEvent(ev); isfound = tree; break; }
    else ev -= nevts;
  }
  return isfound;
}


