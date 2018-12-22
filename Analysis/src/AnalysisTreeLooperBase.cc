#include <algorithm>
#include <utility>
#include <iterator>
#include <cassert>
#include "AnalysisTreeLooperBase.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


AnalysisTreeLooperBase::AnalysisTreeLooperBase() : IvyBase(), sampleIdOpt(AnalysisTreeLooperBase::kNoStorage), maxNEvents(-1) { setExternalProductList(); setExternalProductTree(); }
AnalysisTreeLooperBase::AnalysisTreeLooperBase(AnalysisTree* inTree) : IvyBase(), sampleIdOpt(AnalysisTreeLooperBase::kNoStorage), maxNEvents(-1) { this->addTree(inTree); setExternalProductList(); setExternalProductTree(); }
AnalysisTreeLooperBase::AnalysisTreeLooperBase(std::vector<AnalysisTree*> const& inTreeList) :
  IvyBase(),
  sampleIdOpt(AnalysisTreeLooperBase::kNoStorage),
  treeList(inTreeList),
  maxNEvents(-1)
{
  setExternalProductList();
  setExternalProductTree();
}
AnalysisTreeLooperBase::AnalysisTreeLooperBase(AnalysisSet const* inTreeSet) :
  IvyBase(),
  sampleIdOpt(AnalysisTreeLooperBase::kNoStorage),
  treeList(inTreeSet->getAnalysisTreeList()),
  maxNEvents(-1)
{
  setExternalProductList();
  setExternalProductTree();
}
AnalysisTreeLooperBase::~AnalysisTreeLooperBase(){}

void AnalysisTreeLooperBase::addTree(AnalysisTree* tree){ this->treeList.push_back(tree); }


void AnalysisTreeLooperBase::addExternalIvyObject(TString objname, IvyBase* obj){
  if (!obj) return;
  if (obj == (IvyBase*)this){
    if (verbosity>=TVar::ERROR) MELAerr << "AnalysisTreeLooperBase::addExternalIvyObject: " << objname << " is the same as the AnalysisTreeLooperBase! This object can never be added, so ignoring the association!" << endl;
    return;
  }
  if (externalIvyObjects.find(objname)!=externalIvyObjects.end() && verbosity>=TVar::ERROR) MELAerr << "AnalysisTreeLooperBase::addExternalIvyObject: " << objname << " already exists but will override it regardless." << endl;
  externalIvyObjects[objname] = obj;
}
void AnalysisTreeLooperBase::addExternalFunction(TString fcnname, void(*fcn)(AnalysisTreeLooperBase*, SimpleEntry&)){
  if (!fcn) return;
  if (externalFunctions.find(fcnname)!=externalFunctions.end()) MELAerr << "AnalysisTreeLooperBase::addExternalFunction: " << fcnname << " already exists but will override it regardless." << endl;
  externalFunctions[fcnname] = fcn;
}


void AnalysisTreeLooperBase::setExternalProductList(std::vector<SimpleEntry>* extProductListRef){
  if (extProductListRef) this->productListRef=extProductListRef;
  else this->productListRef=&(this->productList);
}

void AnalysisTreeLooperBase::setExternalProductTree(BaseTree* extTree){
  this->productTree=extTree;
  this->productListRef=&(this->productList); // To make sure product list collects some events before flushing
}

void AnalysisTreeLooperBase::setMaximumEvents(int n){ maxNEvents=n; }

void AnalysisTreeLooperBase::setSampleIdStorageOption(AnalysisTreeLooperBase::SampleIdStorageType opt){ sampleIdOpt=opt; }

void AnalysisTreeLooperBase::addProduct(SimpleEntry& product, unsigned int* ev_rec){
  this->productListRef->push_back(product);
  if (ev_rec) (*ev_rec)++;
}

void AnalysisTreeLooperBase::recordProductsToTree(){
  if (!this->productTree) return;
  BaseTree::writeSimpleEntries(this->productListRef->begin(), this->productListRef->end(), this->productTree);
  this->clearProducts();
}

void AnalysisTreeLooperBase::loop(bool loopSelected, bool loopFailed, bool keepProducts){
  const TString strVarRunNumber = "uint_eventMaker_evtrun_CMS3.obj";
  const TString strVarEventNumber = "ull_eventMaker_evtevent_CMS3.obj";
  // Loop over the trees
  unsigned int ev_acc=0;
  unsigned int ev_rec=0;
  const bool storeSampleIdByRunAndEventNumber = (sampleIdOpt==kStoreByRunAndEventNumber);
  const bool storeSampleIdByHashVal = (sampleIdOpt==kStoreByHashVal);
  vector<unsigned int> loopRecSelList, loopTotalSelList, loopRecFailList, loopTotalFailList;
  vector<unsigned int>::iterator it_loopRecSelList, it_loopTotalSelList, it_loopRecFailList, it_loopTotalFailList;
  if (verbosity>=TVar::INFO && !treeList.empty()){
    loopRecSelList.assign(treeList.size(), 0); it_loopRecSelList=loopRecSelList.begin();
    loopTotalSelList.assign(treeList.size(), 0); it_loopTotalSelList=loopTotalSelList.begin();
    loopRecFailList.assign(treeList.size(), 0); it_loopRecFailList=loopRecFailList.begin();
    loopTotalFailList.assign(treeList.size(), 0); it_loopTotalFailList=loopTotalFailList.begin();
  }
  if (storeSampleIdByRunAndEventNumber){ // Check if RunNumber and EventNumber variables are consumed
    bool doAbort=false;
    if (valuints.find(strVarRunNumber)==valuints.cend()){
      MELAerr << "AnalysisTreeLooperBase::loop: RunNumber is not a consumed variable!" << endl;
      doAbort=true;
    }
    if (valulonglongs.find(strVarEventNumber)==valulonglongs.cend()){
      MELAerr << "AnalysisTreeLooperBase::loop: EventNumber is not a consumed variable!" << endl;
      doAbort=true;
    }
    assert(!doAbort);
  }
  for (AnalysisTree*& tree:treeList){
    // Skip the tree if it cannot be linked
    if (!(this->wrapTree(tree))) continue;

    // Wrap external ivy objects to the current tree
    bool externalObjectsWrapped = true;
    for (auto it_ivy=externalIvyObjects.begin(); it_ivy!=externalIvyObjects.end(); it_ivy++) externalObjectsWrapped &= it_ivy->second->wrapTree(tree);
    if (!externalObjectsWrapped) continue;

    float wgtExternal = 1;
    //AnalysisSet const* associatedSet = tree->getAssociatedSet();
    //if (associatedSet) wgtExternal *= associatedSet->getPermanentWeight(tree);
    if (wgtExternal==0.){
      MELAerr << "AnalysisTreeLooperBase::loop: External weights are 0 for the " << tree->sampleIdentifier << " sample. Skipping..." << endl;
      continue;
    }

    size_t sampleId=0;
    if (storeSampleIdByHashVal){ std::hash<TString> tmphash; sampleId=tmphash(tree->sampleIdentifier); }

    // Loop over selected events
    if (loopSelected){
      MELAout << "AnalysisTreeLooperBase::loop: Looping over " << tree->sampleIdentifier << " selected events" << endl;
      int ev=0;
      const int nevents = tree->getSelectedNEvents();
      while (tree->getSelectedEvent(ev)){
        if (maxNEvents>=0 && (int) ev_rec==maxNEvents) break;
        SimpleEntry product;
        if (tree->isValidEvent()){
          if (this->runEvent(tree, wgtExternal, product)){
            if (keepProducts){
              if (storeSampleIdByHashVal){
                product.setNamedVal("SampleId", sampleId);
                product.setNamedVal("EventNumber", ev_acc);
              }
              else if (storeSampleIdByRunAndEventNumber){
                product.setNamedVal("RunNumber", *(valuints[strVarRunNumber]));
                product.setNamedVal("EventNumber", *(valulonglongs[strVarEventNumber]));
              }
              this->addProduct(product, &ev_rec);
              if (verbosity>=TVar::INFO) (*it_loopRecSelList)++;
            }
          }
        }
        HelperFunctions::progressbar(ev, nevents);
        ev++; ev_acc++;
        if (verbosity>=TVar::INFO) (*it_loopTotalSelList)++;
      }
    }
    // Loop over failed events
    if (loopFailed){
      MELAout << "AnalysisTreeLooperBase::loop: Looping over " << tree->sampleIdentifier << " failed events" << endl;
      int ev=0;
      const int nevents = tree->getFailedNEvents();
      while (tree->getFailedEvent(ev)){
        if (maxNEvents>=0 && (int) ev_rec==maxNEvents) break;
        SimpleEntry product;
        if (tree->isValidEvent()){
          if (this->runEvent(tree, wgtExternal, product)){
            if (keepProducts){
              if (storeSampleIdByHashVal){
                product.setNamedVal("SampleId", sampleId);
                product.setNamedVal("EventNumber", ev_acc);
              }
              else if (storeSampleIdByRunAndEventNumber){
                product.setNamedVal("RunNumber", *(valuints[strVarRunNumber]));
                product.setNamedVal("EventNumber", *(valulonglongs[strVarEventNumber]));
              }
              this->addProduct(product, &ev_rec);
              if (verbosity>=TVar::INFO) (*it_loopRecFailList)++;
            }
          }
        }
        HelperFunctions::progressbar(ev, nevents);
        ev++; ev_acc++;
        if (verbosity>=TVar::INFO) (*it_loopTotalFailList)++;
      }
    }

    // Record products to external tree
    this->recordProductsToTree();

    if (verbosity>=TVar::INFO){
      it_loopRecSelList++;
      it_loopRecFailList++;
      it_loopTotalSelList++;
      it_loopTotalFailList++;
    }
  } // End loop over the trees
  MELAout << "AnalysisTreeLooperBase::loop: Total number of products: " << ev_rec << " / " << ev_acc << endl;
  if (verbosity>=TVar::INFO){
    for (unsigned int it=0; it<treeList.size(); it++){
      MELAout << "\t- AnalysisTreeLooperBase::loop: Total number of selected | failed products in tree " << it << ": "
        << loopRecSelList.at(it) << " / " << loopTotalSelList.at(it)
        << " | "
        << loopRecFailList.at(it) << " / " << loopTotalFailList.at(it)
        << endl;
    }
  }
}

std::vector<SimpleEntry> const& AnalysisTreeLooperBase::getProducts() const{ return *productListRef; }

void AnalysisTreeLooperBase::moveProducts(std::vector<SimpleEntry>& targetColl){
  MELAout << "AnalysisTreeLooperBase::moveProducts: Moving " << productListRef->size() << " products into a list of initial size " << targetColl.size() << endl;
  std::move(productListRef->begin(), productListRef->end(), std::back_inserter(targetColl));
  clearProducts();
  MELAout << "AnalysisTreeLooperBase::moveProducts: Target list final size: " << targetColl.size() << endl;
}

void AnalysisTreeLooperBase::clearProducts(){ std::vector<SimpleEntry> emptyList; std::swap(emptyList, *productListRef); }
