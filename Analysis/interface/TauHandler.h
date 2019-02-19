#ifndef TAUHANDLER_H
#define TAUHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "MuonObject.h"
#include "ElectronObject.h"
#include "TauObject.h"


class TauHandler : public IvyBase{
public:
  typedef TauObject ProductType_t;

protected:
  std::vector<ProductType_t*> productList;

  std::vector<MuonObject*> const* registeredMuons;
  std::vector<ElectronObject*> const* registeredElectrons;

  void clear(){ for (ProductType_t*& prod:productList) delete prod; productList.clear(); } // Do not clear registered objects here

public:
  // Constructors
  TauHandler();

  // Destructors
  ~TauHandler(){ clear(); }

  bool constructTaus();
  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  void registerParticles(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons){ registeredMuons = muons; registeredElectrons = electrons; }

  static void bookBranches(BaseTree* tree);

};


#endif
