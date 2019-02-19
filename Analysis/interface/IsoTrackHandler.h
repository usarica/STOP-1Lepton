#ifndef ISOTRACKHANDLER_H
#define ISOTRACKHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "MuonObject.h"
#include "ElectronObject.h"
#include "IsoTrackObject.h"


class IsoTrackHandler : public IvyBase{
public:
  typedef IsoTrackObject ProductType_t;

protected:
  std::vector<ProductType_t*> productList;

  std::vector<MuonObject*> const* registeredMuons;
  std::vector<ElectronObject*> const* registeredElectrons;
  // No photons since photons don't have tracks!

  void clear(){ for (ProductType_t*& prod:productList) delete prod; productList.clear(); } // Do not clear registered objects here

public:
  // Constructors
  IsoTrackHandler();

  // Destructors
  ~IsoTrackHandler(){ clear(); }

  bool constructIsoTracks();
  std::vector<ProductType_t*> const& getProducts() const{ return productList; }

  void registerParticles(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons){ registeredMuons = muons; registeredElectrons = electrons; }

  static void bookBranches(BaseTree* tree);

};


#endif
