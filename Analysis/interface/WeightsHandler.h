#ifndef WEIGHTSHANDLER_H
#define WEIGHTSHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "WeightsObject.h"
#include "LHEWeightHandler.h"


class WeightsHandler : public IvyBase{
public:
  typedef WeightsObject ProductType_t;

protected:
  ProductType_t* product;
  LHEWeightHandler* weightHandler_DefaultPDF;
  LHEWeightHandler* weightHandler_2016;

  void clear(){ delete product; product=nullptr; }

public:
  // Constructors
  WeightsHandler();

  // Destructors
  ~WeightsHandler(){ clear(); delete weightHandler_DefaultPDF; delete weightHandler_2016; }

  bool constructWeights();
  ProductType_t const* getProduct() const{ return product; }

  void bookBranches(BaseTree* tree);

};


#endif
