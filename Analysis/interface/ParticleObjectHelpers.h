#ifndef PARTICLEOBJECTHELPERS_H
#define PARTICLEOBJECTHELPERS_H

#include <vector>
#include <algorithm>
#include "ParticleObject.h"


namespace ParticleObjectHelpers{
  template<typename T> bool objHasGreaterPt(T const& earlier, T const& later);
  template<typename T> bool ptrHasGreaterPt(T const* earlier, T const* later);

  template<typename T> void sortByGreaterPt(std::vector<T>& vec);
  template<typename T> void sortByGreaterPt(std::vector<T*>& vec);

}

template<typename T> bool ParticleObjectHelpers::objHasGreaterPt(T const& earlier, T const& later){ return (earlier.getFinalMomentum().Pt() > later.getFinalMomentum().Pt()); }
template<typename T> bool ParticleObjectHelpers::ptrHasGreaterPt(T const* earlier, T const* later){ return (earlier && ((later && earlier->getFinalMomentum().Pt() > later->getFinalMomentum().Pt()) || !later)); }

template<typename T> void ParticleObjectHelpers::sortByGreaterPt(std::vector<T>& vec){ std::sort(vec.begin(), vec.end(), ParticleObjectHelpers::objHasGreaterPt<T>); }
template<typename T> void ParticleObjectHelpers::sortByGreaterPt(std::vector<T*>& vec){ std::sort(vec.begin(), vec.end(), ParticleObjectHelpers::ptrHasGreaterPt<T>); }


#endif
