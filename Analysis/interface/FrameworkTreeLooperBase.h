#ifndef FRAMEWORKTREELOOPERBASE_H
#define FRAMEWORKTREELOOPERBASE_H

#include "IvyBase.h"
#include "FrameworkSet.h"
#include "HelperFunctions.h"


class FrameworkTreeLooperBase : public IvyBase{
public:
  enum SampleIdStorageType{
    kNoStorage,
    kStoreByRunAndEventNumber,
    kStoreByHashVal
  };

protected:
  SampleIdStorageType sampleIdOpt; // When not kNoStorage, stores a sample identifier and original tree entry index in each output event

  // List of trees to loop over
  std::vector<FrameworkTree*> treeList;

  // Max. events to process
  int maxNEvents;

  // External dependencies
  std::unordered_map<TString, IvyBase*> externalIvyObjects;
  std::unordered_map<TString, void(*)(FrameworkTreeLooperBase*, SimpleEntry&)> externalFunctions;

  // List of products
  std::vector<SimpleEntry> productList;
  std::vector<SimpleEntry>* productListRef;
  BaseTree* productTree;
  void addProduct(SimpleEntry& product, unsigned int* ev_rec=nullptr);

  // Flush product list into tree
  void recordProductsToTree();

  // Abstract function to loop over a single event
  virtual bool runEvent(FrameworkTree* tree, float const& externalWgt, SimpleEntry& product)=0;

public:
  // Constructors
  FrameworkTreeLooperBase();
  FrameworkTreeLooperBase(FrameworkTree* inTree);
  FrameworkTreeLooperBase(std::vector<FrameworkTree*> const& inTreeList);
  FrameworkTreeLooperBase(FrameworkSet const* inTreeSet);

  void addTree(FrameworkTree* tree);

  // Destructors
  virtual ~FrameworkTreeLooperBase();

  // Add the necessary objects
  void addExternalIvyObject(TString objname, IvyBase* obj);
  void addExternalFunction(TString fcnname, void(*fcn)(FrameworkTreeLooperBase*, SimpleEntry&));
  void setExternalProductList(std::vector<SimpleEntry>* extProductListRef=nullptr);
  void setExternalProductTree(BaseTree* extTree=nullptr);

  // Max. events
  void setMaximumEvents(int n);

  // Sample id storage option
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
