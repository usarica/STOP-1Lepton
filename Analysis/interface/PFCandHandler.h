#ifndef PFCANDHANDLER_H
#define PFCANDHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "PFCandObject.h"


class PFCandHandler : public IvyBase{
public:
  typedef PFCandObject ProductType_t;

protected:
  std::vector<ProductType_t*> productList;

  void clear(){ for (ProductType_t*& prod:productList) delete prod; productList.clear(); }

public:
  // Constructors
  PFCandHandler();

  // Destructors
  ~PFCandHandler(){ clear(); }

  bool constructPFCands();
  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  static void bookBranches(BaseTree* tree);

};


#endif
