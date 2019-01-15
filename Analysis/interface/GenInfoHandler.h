#ifndef GENINFOHANDLER_H
#define GENINFOHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "GenParticleObject.h"


class GenInfoHandler : public IvyBase{
protected:
  bool doEventInfo;
  bool doParticleInfo;
  std::vector<GenParticleObject*> genparticles;
  GenEventInfo* geninfo;

  void clear();

public:
  // Constructors
  GenInfoHandler();

  // Destructors
  ~GenInfoHandler(){ clear(); }

  std::vector<GenParticleObject*>& getGenParticles(){ return genparticles; }
  std::vector<GenParticleObject*> const& getGenParticles() const{ return genparticles; }

  GenEventInfo* getGenInfo(){ return geninfo; }
  GenEventInfo const* getGenInfo() const{ return geninfo; }

  bool constructGenInfo();

  void setEventInfoFlag(bool flag){ doEventInfo = flag; }
  void setParticleInfoFlag(bool flag){ doParticleInfo = flag; }

  void bookBranches(BaseTree* tree);

};


#endif
