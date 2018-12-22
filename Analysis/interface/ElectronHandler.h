#ifndef ELECTRONHANDLER_H
#define ELECTRONHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "ElectronObject.h"


class ElectronHandler : public IvyBase{
public:
  typedef ElectronObject ProductType_t;

protected:
  std::vector<ProductType_t*> productList;

  void clear(){ for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

public:
  // Constructors
  ElectronHandler() : IvyBase() {}

  // Destructors
  ~ElectronHandler(){ clear(); }

  bool constructElectrons();

  bool wrapTree(BaseTree* tree);

  static void bookBranches(BaseTree* tree);

};


#endif
