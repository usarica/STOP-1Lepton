#include <cassert>
#include <algorithm>
#include <utility>
#include <iterator>
#include <chrono>
#include "FrameworkVariables.hh"
#include "FrameworkTreeLooperBase.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


FrameworkTreeLooperBase::FrameworkTreeLooperBase() : IvyBase(), sampleIdOpt(FrameworkTreeLooperBase::kNoStorage), maxNEvents(-1), recordEveryNEvents(-1){ setExternalProductList(); setExternalProductTree(); }
FrameworkTreeLooperBase::FrameworkTreeLooperBase(FrameworkTree* inTree) : IvyBase(), sampleIdOpt(FrameworkTreeLooperBase::kNoStorage), maxNEvents(-1), recordEveryNEvents(-1) { this->addTree(inTree); setExternalProductList(); setExternalProductTree(); }
FrameworkTreeLooperBase::FrameworkTreeLooperBase(std::vector<FrameworkTree*> const& inTreeList) :
  IvyBase(),
  sampleIdOpt(FrameworkTreeLooperBase::kNoStorage),
  treeList(inTreeList),
  maxNEvents(-1),
  recordEveryNEvents(-1)
{
  setExternalProductList();
  setExternalProductTree();
}
FrameworkTreeLooperBase::FrameworkTreeLooperBase(FrameworkSet const* inTreeSet) :
  IvyBase(),
  sampleIdOpt(FrameworkTreeLooperBase::kNoStorage),
  treeList(inTreeSet->getFrameworkTreeList()),
  maxNEvents(-1),
  recordEveryNEvents(-1)
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

void FrameworkTreeLooperBase::setRecordEveryNEvents(int n){ recordEveryNEvents=n; }

void FrameworkTreeLooperBase::setSampleIdStorageOption(FrameworkTreeLooperBase::SampleIdStorageType opt){
  sampleIdOpt=opt;
  if (sampleIdOpt==kStoreByRunAndEventNumber){
    this->addConsumed<RunNumber_t>(_event_RunNumber_);
    this->addConsumed<Lumisection_t>(_event_Lumisection_);
    this->addConsumed<EventNumber_t>(_event_EventNumber_);
  }
}

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

  auto time_start = std::chrono::steady_clock::now();
  for (FrameworkTree*& tree:treeList){
    if (storeSampleIdByRunAndEventNumber){
      tree->bookEDMBranch<RunNumber_t>(_event_RunNumber_, 0);
      tree->bookEDMBranch<Lumisection_t>(_event_Lumisection_, 0);
      tree->bookEDMBranch<EventNumber_t>(_event_EventNumber_, 0);
    }

    // Skip if maximum events are already reached
    if (maxNEvents>=0 && (int) ev_rec==maxNEvents) break;

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
      auto time_begin = std::chrono::steady_clock::now();
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
                // Acquire run and event numbers, and the lumisection
                RunNumber_t theRunNumberVal=0;
                Lumisection_t theLumisectionVal=0;
                EventNumber_t theEventNumberVal=0;
                bool allSampleIdVariablesPresent = (
                  this->getConsumedValue(_event_RunNumber_, theRunNumberVal)
                  &&
                  this->getConsumedValue(_event_Lumisection_, theLumisectionVal)
                  &&
                  this->getConsumedValue(_event_EventNumber_, theEventNumberVal)
                  );
                assert(allSampleIdVariablesPresent);

                product.setNamedVal("RunNumber", theRunNumberVal);
                product.setNamedVal("LumiSection", theLumisectionVal);
                product.setNamedVal("EventNumber", theEventNumberVal);
              }
              this->addProduct(product, &ev_rec);
              if (verbosity>=TVar::INFO) (*it_loopRecSelList)++;
            }
          }
        }
        HelperFunctions::progressbar(ev, nevents);
        ev++; ev_acc++;
        if (verbosity>=TVar::INFO) (*it_loopTotalSelList)++;

        // Record products to external tree if recordEveryNEvents>0 is specified
        if (recordEveryNEvents>0 && ev_rec%recordEveryNEvents==0) this->recordProductsToTree();
      }
      auto time_end = std::chrono::steady_clock::now();
      float duration = (std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_begin)).count(); // In miliseconds
      MELAout << "\t- Average rate of processing: " << float(ev)/duration << " kHz" << endl;
    }
    // Loop over failed events
    if (loopFailed){
      MELAout << "FrameworkTreeLooperBase::loop: Looping over " << tree->sampleIdentifier << " failed events" << endl;
      int ev=0;
      const int nevents = tree->getFailedNEvents();
      auto time_begin = std::chrono::steady_clock::now();
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
                // Acquire run and event numbers, and the lumisection
                RunNumber_t theRunNumberVal=0;
                Lumisection_t theLumisectionVal=0;
                EventNumber_t theEventNumberVal=0;
                bool allSampleIdVariablesPresent = (
                  this->getConsumedValue(_event_RunNumber_, theRunNumberVal)
                  &&
                  this->getConsumedValue(_event_Lumisection_, theLumisectionVal)
                  &&
                  this->getConsumedValue(_event_EventNumber_, theEventNumberVal)
                  );
                assert(allSampleIdVariablesPresent);

                product.setNamedVal("RunNumber", theRunNumberVal);
                product.setNamedVal("LumiSection", theLumisectionVal);
                product.setNamedVal("EventNumber", theEventNumberVal);
              }
              this->addProduct(product, &ev_rec);
              if (verbosity>=TVar::INFO) (*it_loopRecFailList)++;
            }
          }
        }
        HelperFunctions::progressbar(ev, nevents);
        ev++; ev_acc++;
        if (verbosity>=TVar::INFO) (*it_loopTotalFailList)++;

        // Record products to external tree if recordEveryNEvents>0 is specified
        if (recordEveryNEvents>0 && ev_rec%recordEveryNEvents==0) this->recordProductsToTree();
      }
      auto time_end = std::chrono::steady_clock::now();
      float duration = (std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_begin)).count(); // In miliseconds
      MELAout << "\t- Average rate of processing: " << float(ev)/duration << " kHz" << endl;
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
  {
    auto time_end = std::chrono::steady_clock::now();
    float duration = (std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start)).count(); // In miliseconds
    MELAout << "\t- Average rate of processing: " << float(ev_acc)/duration << " kHz" << endl;
  }
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
