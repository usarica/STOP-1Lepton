#ifndef BASETREELOOPER_H
#define BASETREELOOPER_H

#include "AnalysisSet.h"
#include "HelperFunctions.h"
#include "CMSLorentzVector.h"


class BaseTreeLooper{
public:
  enum SampleIdStorageType{
    kNoStorage,
    kStoreByRunAndEventNumber,
    kStoreByHashVal
  };

protected:
  SampleIdStorageType sampleIdOpt; // When not kNoStorage, stores a sample identifier and original tree entry index in each output event

  // List of trees to loop over
  std::vector<AnalysisTree*> treeList;

  // Max. events to process
  int maxNEvents;

  // Verbosity flag
  bool verbose;

  // Consumes
  std::unordered_map<TString, short*> valshorts;
  std::unordered_map<TString, unsigned int*> valuints;
  std::unordered_map<TString, int*> valints;
  std::unordered_map<TString, unsigned long*> valulongs;
  std::unordered_map<TString, long*> vallongs;
  std::unordered_map<TString, unsigned long long*> valulonglongs;
  std::unordered_map<TString, long long*> vallonglongs;
  std::unordered_map<TString, float*> valfloats;
  std::unordered_map<TString, double*> valdoubles;
  std::unordered_map<TString, CMSLorentzVector*> valCMSLorentzVectors;

  std::unordered_map<TString, std::vector<short>*> valVshorts;
  std::unordered_map<TString, std::vector<unsigned int>*> valVuints;
  std::unordered_map<TString, std::vector<int>*> valVints;
  std::unordered_map<TString, std::vector<unsigned long>*> valVulongs;
  std::unordered_map<TString, std::vector<long>*> valVlongs;
  std::unordered_map<TString, std::vector<unsigned long long>*> valVulonglongs;
  std::unordered_map<TString, std::vector<long long>*> valVlonglongs;
  std::unordered_map<TString, std::vector<float>*> valVfloats;
  std::unordered_map<TString, std::vector<double>*> valVdoubles;
  std::unordered_map<TString, std::vector<CMSLorentzVector>*> valVCMSLorentzVectors;

  template<typename T> bool linkConsumed(AnalysisTree* tree);
  bool linkConsumes(AnalysisTree* tree);

  // Get consumed map
  template<typename T> void getConsumedMap(std::unordered_map<TString, T*>*& theMap);
  template<typename T> void getConsumedMap(std::unordered_map<TString, T*> const*& theMap) const;

  // Get consumed
  template<typename T> bool getConsumed(TString name, T const*& val) const;

  // External dependencies
  std::unordered_map<TString, void(*)(BaseTreeLooper*, SimpleEntry&)> externalFunctions;

  // List of products
  std::vector<SimpleEntry> productList;
  std::vector<SimpleEntry>* productListRef;
  BaseTree* productTree;
  void addProduct(SimpleEntry& product, unsigned int* ev_rec=nullptr);

  // Flush product list into tree
  void recordProductsToTree();

  // Abstract function to loop over a single event
  virtual bool runEvent(AnalysisTree* tree, float const& externalWgt, SimpleEntry& product)=0;

public:
  // Constructors
  BaseTreeLooper();
  BaseTreeLooper(AnalysisTree* inTree);
  BaseTreeLooper(std::vector<AnalysisTree*> const& inTreeList);
  BaseTreeLooper(AnalysisSet const* inTreeSet);
  void addTree(AnalysisTree* tree);

  // Destructors
  virtual ~BaseTreeLooper();

  // Set verbosity
  void setVerbosity(bool flag);

  // Add the necessary objects
  template<typename T> void addConsumed(TString name);
  void addExternalFunction(TString fcnname, void(*fcn)(BaseTreeLooper*, SimpleEntry&));
  void setExternalProductList(std::vector<SimpleEntry>* extProductListRef=nullptr);
  void setExternalProductTree(BaseTree* extTree=nullptr);

  // Max. events
  void setMaximumEvents(int n);

  // Sample id storage option
  // POWHEG can be stored by mH, but might be better to use hash in others
  void setSampleIdStorageOption(SampleIdStorageType opt);

  // Function to loop over the tree list
  virtual void loop(bool loopSelected, bool loopFailed, bool keepProducts);

  // Get the products
  std::vector<SimpleEntry> const& getProducts() const;
  // Move the products
  void moveProducts(std::vector<SimpleEntry>& targetColl);
  // Clear the products
  void clearProducts();

};


#endif
