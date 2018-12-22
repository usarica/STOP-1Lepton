#include <algorithm>
#include <utility>
#include <iterator>
#include <cassert>
#include "BaseTreeLooper.h"
#include "BaseTreeLooper.hpp"


using namespace std;


BaseTreeLooper::BaseTreeLooper() : sampleIdOpt(BaseTreeLooper::kNoStorage), maxNEvents(-1), verbose(false) { setExternalProductList(); setExternalProductTree(); }
BaseTreeLooper::BaseTreeLooper(AnalysisTree* inTree) : sampleIdOpt(BaseTreeLooper::kNoStorage), maxNEvents(-1), verbose(false) { this->addTree(inTree); setExternalProductList(); setExternalProductTree(); }
BaseTreeLooper::BaseTreeLooper(std::vector<AnalysisTree*> const& inTreeList) :
  sampleIdOpt(BaseTreeLooper::kNoStorage),
  treeList(inTreeList),
  maxNEvents(-1),
  verbose(false)
{
  setExternalProductList();
  setExternalProductTree();
}
BaseTreeLooper::BaseTreeLooper(AnalysisSet const* inTreeSet) :
  sampleIdOpt(BaseTreeLooper::kNoStorage),
  treeList(inTreeSet->getAnalysisTreeList()),
  maxNEvents(-1),
  verbose(false)
{
  setExternalProductList();
  setExternalProductTree();
}
BaseTreeLooper::~BaseTreeLooper(){}

void BaseTreeLooper::setVerbosity(bool flag){ verbose = flag; }

void BaseTreeLooper::addTree(AnalysisTree* tree){ this->treeList.push_back(tree); }


void BaseTreeLooper::addExternalFunction(TString fcnname, void(*fcn)(BaseTreeLooper*, SimpleEntry&)){
  if (!fcn) return;
  if (externalFunctions.find(fcnname)!=externalFunctions.end()) MELAerr << "BaseTreeLooper::addExternalFunction: " << fcnname << " already exists but will override it regardless." << endl;
  externalFunctions[fcnname] = fcn;
}


void BaseTreeLooper::setExternalProductList(std::vector<SimpleEntry>* extProductListRef){
  if (extProductListRef) this->productListRef=extProductListRef;
  else this->productListRef=&(this->productList);
}

void BaseTreeLooper::setExternalProductTree(BaseTree* extTree){
  this->productTree=extTree;
  this->productListRef=&(this->productList); // To make sure product list collects some events before flushing
}

void BaseTreeLooper::setMaximumEvents(int n){ maxNEvents=n; }

void BaseTreeLooper::setSampleIdStorageOption(BaseTreeLooper::SampleIdStorageType opt){ sampleIdOpt=opt; }

void BaseTreeLooper::addProduct(SimpleEntry& product, unsigned int* ev_rec){
  this->productListRef->push_back(product);
  if (ev_rec) (*ev_rec)++;
}

void BaseTreeLooper::recordProductsToTree(){
  if (!this->productTree) return;
  BaseTree::writeSimpleEntries(this->productListRef->begin(), this->productListRef->end(), this->productTree);
  this->clearProducts();
}

bool BaseTreeLooper::linkConsumes(AnalysisTree* tree){
  bool process = tree->isValid();
  if (process){
    process &= this->linkConsumed<short>(tree);
    process &= this->linkConsumed<unsigned int>(tree);
    process &= this->linkConsumed<int>(tree);
    process &= this->linkConsumed<unsigned long>(tree);
    process &= this->linkConsumed<long>(tree);
    process &= this->linkConsumed<long long>(tree);
    process &= this->linkConsumed<float>(tree);
    process &= this->linkConsumed<double>(tree);
    process &= this->linkConsumed<std::vector<short>>(tree);
    process &= this->linkConsumed<std::vector<unsigned int>>(tree);
    process &= this->linkConsumed<std::vector<int>>(tree);
    process &= this->linkConsumed<std::vector<unsigned long>>(tree);
    process &= this->linkConsumed<std::vector<long>>(tree);
    process &= this->linkConsumed<std::vector<long long>>(tree);
    process &= this->linkConsumed<std::vector<float>>(tree);
    process &= this->linkConsumed<std::vector<double>>(tree);
    // Silence unused branches
    tree->silenceUnused();
  }
  if (!process) MELAerr << "BaseTreeLooper::linkConsumes: Linking failed for some reason for tree " << tree->sampleIdentifier << endl;
  return process;
}

void BaseTreeLooper::loop(bool loopSelected, bool loopFailed, bool keepProducts){
  const TString strVarRunNumber = "uint_eventMaker_evtrun_CMS3.obj";
  const TString strVarEventNumber = "ull_eventMaker_evtevent_CMS3.obj";
  // Loop over the trees
  unsigned int ev_acc=0;
  unsigned int ev_rec=0;
  const bool storeSampleIdByRunAndEventNumber = (sampleIdOpt==kStoreByRunAndEventNumber);
  const bool storeSampleIdByHashVal = (sampleIdOpt==kStoreByHashVal);
  vector<unsigned int> loopRecSelList, loopTotalSelList, loopRecFailList, loopTotalFailList;
  vector<unsigned int>::iterator it_loopRecSelList, it_loopTotalSelList, it_loopRecFailList, it_loopTotalFailList;
  if (verbose && !treeList.empty()){
    loopRecSelList.assign(treeList.size(), 0); it_loopRecSelList=loopRecSelList.begin();
    loopTotalSelList.assign(treeList.size(), 0); it_loopTotalSelList=loopTotalSelList.begin();
    loopRecFailList.assign(treeList.size(), 0); it_loopRecFailList=loopRecFailList.begin();
    loopTotalFailList.assign(treeList.size(), 0); it_loopTotalFailList=loopTotalFailList.begin();
  }
  if (storeSampleIdByRunAndEventNumber){ // Check if RunNumber and EventNumber variables are consumed
    bool doAbort=false;
    if (valuints.find(strVarRunNumber)==valuints.cend()){
      MELAerr << "BaseTreeLooper::loop: RunNumber is not a consumed variable!" << endl;
      doAbort=true;
    }
    if (valulonglongs.find(strVarEventNumber)==valulonglongs.cend()){
      MELAerr << "BaseTreeLooper::loop: EventNumber is not a consumed variable!" << endl;
      doAbort=true;
    }
    assert(!doAbort);
  }
  for (AnalysisTree*& tree:treeList){
    // Skip the tree if it cannot be linked
    if (!(this->linkConsumes(tree))) continue;

    float wgtExternal = 1;
    //AnalysisSet const* associatedSet = tree->getAssociatedSet();
    //if (associatedSet) wgtExternal *= associatedSet->getPermanentWeight(tree);
    if (wgtExternal==0.){
      MELAerr << "BaseTreeLooper::loop: External weights are 0 for the " << tree->sampleIdentifier << " sample. Skipping..." << endl;
      continue;
    }

    size_t sampleId=0;
    if (storeSampleIdByHashVal){ std::hash<TString> tmphash; sampleId=tmphash(tree->sampleIdentifier); }

    // Loop over selected events
    if (loopSelected){
      MELAout << "BaseTreeLooper::loop: Looping over " << tree->sampleIdentifier << " selected events" << endl;
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
              if (verbose) (*it_loopRecSelList)++;
            }
          }
        }
        HelperFunctions::progressbar(ev, nevents);
        ev++; ev_acc++;
        if (verbose) (*it_loopTotalSelList)++;
      }
    }
    // Loop over failed events
    if (loopFailed){
      MELAout << "BaseTreeLooper::loop: Looping over " << tree->sampleIdentifier << " failed events" << endl;
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
              if (verbose) (*it_loopRecFailList)++;
            }
          }
        }
        HelperFunctions::progressbar(ev, nevents);
        ev++; ev_acc++;
        if (verbose) (*it_loopTotalFailList)++;
      }
    }

    // Record products to external tree
    this->recordProductsToTree();

    if (verbose){
      it_loopRecSelList++;
      it_loopRecFailList++;
      it_loopTotalSelList++;
      it_loopTotalFailList++;
    }
  } // End loop over the trees
  MELAout << "BaseTreeLooper::loop: Total number of products: " << ev_rec << " / " << ev_acc << endl;
  if (verbose){
    for (unsigned int it=0; it<treeList.size(); it++){
      MELAout << "\t- BaseTreeLooper::loop: Total number of selected | failed products in tree " << it << ": "
        << loopRecSelList.at(it) << " / " << loopTotalSelList.at(it)
        << " | "
        << loopRecFailList.at(it) << " / " << loopTotalFailList.at(it)
        << endl;
    }
  }
}

std::vector<SimpleEntry> const& BaseTreeLooper::getProducts() const{ return *productListRef; }

void BaseTreeLooper::moveProducts(std::vector<SimpleEntry>& targetColl){
  MELAout << "BaseTreeLooper::moveProducts: Moving " << productListRef->size() << " products into a list of initial size " << targetColl.size() << endl;
  std::move(productListRef->begin(), productListRef->end(), std::back_inserter(targetColl));
  clearProducts();
  MELAout << "BaseTreeLooper::moveProducts: Target list final size: " << targetColl.size() << endl;
}

void BaseTreeLooper::clearProducts(){ std::vector<SimpleEntry> emptyList; std::swap(emptyList, *productListRef); }
