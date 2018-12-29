#include <algorithm>
#include "SampleHelpers.h"
#include "FrameworkSet.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


FrameworkSet::FrameworkSet(FrameworkOptionParser const& opts, const TString treename){
  addFrameworkTreeList(opts, treename);
}
FrameworkSet::~FrameworkSet(){
  for (auto& tree:treeList) delete tree;
  treeList.clear();
}

bool FrameworkSet::addFrameworkTree(FrameworkOptionParser const& opts, const TString& fname, const TString treename){
  FrameworkTree* tree = new FrameworkTree(opts, fname, treename);
  if (tree->isValid()){
    tree->silenceUnused();
    treeList.push_back(tree);
  }
  else{ delete tree; tree=nullptr; }
  if (!tree) MELAerr << "FrameworkSet::addFrameworkTree: Tree from file " << fname << " is invalid!" << endl;
  else MELAout << "FrameworkSet::addFrameworkTree(" << fname << "::" << treename << ") is successful!" << endl;
  if (tree) tree->setAssociatedSet(this);
  return (tree!=nullptr);
}
bool FrameworkSet::addFrameworkTreeList(FrameworkOptionParser const& opts, const TString treename){
  bool res=true;
  TString dirname = SampleHelpers::getDatasetDirectoryName(opts);
  std::vector<TString> fnamelist = SampleHelpers::lsdir(dirname);
  for (auto const& s:fnamelist){ if (s.Contains(".root")) res &= addFrameworkTree(opts, s, treename); }
  return res;
}

bool FrameworkSet::dissociateFrameworkTree(FrameworkTree*& tree){
  if (!tree) return false;
  auto it = std::find(treeList.begin(), treeList.end(), tree);
  if (it!=treeList.end()){ treeList.erase(it); tree->setAssociatedSet(nullptr); return true; }
  return false;
}
bool FrameworkSet::associateFrameworkTree(FrameworkTree*& tree){
  if (!tree) return false;
  FrameworkSet* theSet = tree->getAssociatedSet();
  if (theSet==this) return true;
  else if (theSet){
    theSet->dissociateFrameworkTree(tree);
    treeList.push_back(tree);
    tree->setAssociatedSet(this);
    return true;
  }
  return false;
}

std::vector<FrameworkTree*> FrameworkSet::getFrameworkTree(TString sampleid) const{
  std::vector<FrameworkTree*> res;
  for (auto tree:treeList){ if (tree->sampleIdentifier==sampleid) res.push_back(tree); }
  return res;
}
const std::vector<FrameworkTree*>& FrameworkSet::getFrameworkTreeList() const{ return treeList; }
std::vector<FrameworkTree*>& FrameworkSet::getFrameworkTreeList(){ return treeList; }

FrameworkTree* FrameworkSet::getSelectedEvent(const int evid){
  int ev = evid;
  FrameworkTree* isfound=nullptr;
  for (auto& tree:treeList){
    int nevts = tree->getSelectedNEvents();
    if (ev<nevts){ tree->getSelectedEvent(ev); isfound = tree; break; }
    else ev -= nevts;
  }
  return isfound;
}
FrameworkTree* FrameworkSet::getFailedEvent(const int evid){
  int ev = evid;
  FrameworkTree* isfound=nullptr;
  for (auto& tree:treeList){
    int nevts = tree->getFailedNEvents();
    if (ev<nevts){ tree->getFailedEvent(ev); isfound = tree; break; }
    else ev -= nevts;
  }
  return isfound;
}


