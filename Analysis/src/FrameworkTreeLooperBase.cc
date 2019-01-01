#include <algorithm>
#include <utility>
#include <iterator>
#include <cassert>
#include "FrameworkVariables.hh"
#include "FrameworkTreeLooperBase.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


FrameworkTreeLooperBase::FrameworkTreeLooperBase() : IvyBase(), sampleIdOpt(FrameworkTreeLooperBase::kNoStorage), maxNEvents(-1) { setExternalProductList(); setExternalProductTree(); }
FrameworkTreeLooperBase::FrameworkTreeLooperBase(FrameworkTree* inTree) : IvyBase(), sampleIdOpt(FrameworkTreeLooperBase::kNoStorage), maxNEvents(-1) { this->addTree(inTree); setExternalProductList(); setExternalProductTree(); }
FrameworkTreeLooperBase::FrameworkTreeLooperBase(std::vector<FrameworkTree*> const& inTreeList) :
  IvyBase(),
  sampleIdOpt(FrameworkTreeLooperBase::kNoStorage),
  treeList(inTreeList),
  maxNEvents(-1)
{
  setExternalProductList();
  setExternalProductTree();
}
FrameworkTreeLooperBase::FrameworkTreeLooperBase(FrameworkSet const* inTreeSet) :
  IvyBase(),
  sampleIdOpt(FrameworkTreeLooperBase::kNoStorage),
  treeList(inTreeSet->getFrameworkTreeList()),
  maxNEvents(-1)
{
  setExternalProductList();
  setExternalProductTree();
}
FrameworkTreeLooperBase::~FrameworkTreeLooperBase(){}

void FrameworkTreeLooperBase::addTree(FrameworkTree* tree){ this->treeList.push_back(tree); }


void FrameworkTreeLooperBase::addExternalIvyObject(TString objname, IvyBase* obj){
  if (!obj) return;
  if (obj == (IvyBase*)this){
    if (verbosity>=TVar::ERROR) MELAerr << "FrameworkTreeLooperBase::addExternalIvyObject: " << objname << " is the same as the FrameworkTreeLooperBase! This object can never be added, so ignoring the association!" << endl;
    return;
  }
  if (externalIvyObjects.find(objname)!=externalIvyObjects.end() && verbosity>=TVar::ERROR) MELAerr << "FrameworkTreeLooperBase::addExternalIvyObject: " << objname << " already exists but will override it regardless." << endl;
  externalIvyObjects[objname] = obj;
}
void FrameworkTreeLooperBase::addExternalScaleFactorHandler(TString objname, ScaleFactorHandlerBase* obj){
  if (!obj) return;
  if (externalScaleFactorHandlers.find(objname)!=externalScaleFactorHandlers.end() && verbosity>=TVar::ERROR) MELAerr << "FrameworkTreeLooperBase::addExternalScaleFactorHandler: " << objname << " already exists but will override it regardless." << endl;
  externalScaleFactorHandlers[objname] = obj;
}
void FrameworkTreeLooperBase::addExternalFunction(TString fcnname, void(*fcn)(FrameworkTreeLooperBase*, SimpleEntry&)){
  if (!fcn) return;
  if (externalFunctions.find(fcnname)!=externalFunctions.end()) MELAerr << "FrameworkTreeLooperBase::addExternalFunction: " << fcnname << " already exists but will override it regardless." << endl;
  externalFunctions[fcnname] = fcn;
}


void FrameworkTreeLooperBase::setExternalProductList(std::vector<SimpleEntry>* extProductListRef){
  if (extProductListRef) this->productListRef=extProductListRef;
  else this->productListRef=&(this->productList);
}

void FrameworkTreeLooperBase::setExternalProductTree(BaseTree* extTree){
  this->productTree=extTree;
  this->productListRef=&(this->productList); // To make sure product list collects some events before flushing
}

void FrameworkTreeLooperBase::setMaximumEvents(int n){ maxNEvents=n; }

void FrameworkTreeLooperBase::setSampleIdStorageOption(FrameworkTreeLooperBase::SampleIdStorageType opt){ sampleIdOpt=opt; }

void FrameworkTreeLooperBase::addProduct(SimpleEntry& product, unsigned int* ev_rec){
  this->productListRef->push_back(product);
  if (ev_rec) (*ev_rec)++;
}

void FrameworkTreeLooperBase::recordProductsToTree(){
  if (!this->productTree) return;
  BaseTree::writeSimpleEntries(this->productListRef->begin(), this->productListRef->end(), this->productTree);
  this->clearProducts();
}

void FrameworkTreeLooperBase::loop(bool loopSelected, bool loopFailed, bool keepProducts){
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
    if (valuints.find(_event_RunNumber_)==valuints.cend()){
      MELAerr << "FrameworkTreeLooperBase::loop: RunNumber is not a consumed variable!" << endl;
      doAbort=true;
    }
    if (valulonglongs.find(_event_EventNumber_)==valulonglongs.cend()){
      MELAerr << "FrameworkTreeLooperBase::loop: EventNumber is not a consumed variable!" << endl;
      doAbort=true;
    }
    assert(!doAbort);
  }
  for (FrameworkTree*& tree:treeList){
    // Skip the tree if it cannot be linked
    if (!(this->wrapTree(tree))) continue;

    // Setup tree-specific dataset details
    SampleHelpers::setupUsingOptions(tree->getOptions());

    // Wrap external ivy objects to the current tree
    bool externalObjectsWrapped = true;
    for (auto it_ext=externalIvyObjects.begin(); it_ext!=externalIvyObjects.end(); it_ext++) externalObjectsWrapped &= it_ext->second->wrapTree(tree);
    if (!externalObjectsWrapped) continue;

    bool externalSFHandlersSetup = true;
    for (auto it_ext=externalScaleFactorHandlers.begin(); it_ext!=externalScaleFactorHandlers.end(); it_ext++) externalSFHandlersSetup &= it_ext->second->setup();
    if (!externalSFHandlersSetup) continue;

    float wgtExternal = 1;
    //FrameworkSet const* associatedSet = tree->getAssociatedSet();
    //if (associatedSet) wgtExternal *= associatedSet->getPermanentWeight(tree);
    if (wgtExternal==0.){
      MELAerr << "FrameworkTreeLooperBase::loop: External weights are 0 for the " << tree->sampleIdentifier << " sample. Skipping..." << endl;
      continue;
    }

    size_t sampleId=0;
    if (storeSampleIdByHashVal){ std::hash<TString> tmphash; sampleId=tmphash(tree->sampleIdentifier); }

    // Loop over selected events
    if (loopSelected){
      MELAout << "FrameworkTreeLooperBase::loop: Looping over " << tree->sampleIdentifier << " selected events" << endl;
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
                product.setNamedVal("RunNumber", *(valuints[_event_RunNumber_]));
                product.setNamedVal("EventNumber", *(valulonglongs[_event_EventNumber_]));
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
      MELAout << "FrameworkTreeLooperBase::loop: Looping over " << tree->sampleIdentifier << " failed events" << endl;
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
                product.setNamedVal("RunNumber", *(valuints[_event_RunNumber_]));
                product.setNamedVal("EventNumber", *(valulonglongs[_event_EventNumber_]));
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
  MELAout << "FrameworkTreeLooperBase::loop: Total number of products: " << ev_rec << " / " << ev_acc << endl;
  if (verbosity>=TVar::INFO){
    for (unsigned int it=0; it<treeList.size(); it++){
      MELAout << "\t- FrameworkTreeLooperBase::loop: Total number of selected | failed products in tree " << it << ": "
        << loopRecSelList.at(it) << " / " << loopTotalSelList.at(it)
        << " | "
        << loopRecFailList.at(it) << " / " << loopTotalFailList.at(it)
        << endl;
    }
  }
}

std::vector<SimpleEntry> const& FrameworkTreeLooperBase::getProducts() const{ return *productListRef; }

void FrameworkTreeLooperBase::moveProducts(std::vector<SimpleEntry>& targetColl){
  MELAout << "FrameworkTreeLooperBase::moveProducts: Moving " << productListRef->size() << " products into a list of initial size " << targetColl.size() << endl;
  std::move(productListRef->begin(), productListRef->end(), std::back_inserter(targetColl));
  clearProducts();
  MELAout << "FrameworkTreeLooperBase::moveProducts: Target list final size: " << targetColl.size() << endl;
}

void FrameworkTreeLooperBase::clearProducts(){ std::vector<SimpleEntry> emptyList; std::swap(emptyList, *productListRef); }
