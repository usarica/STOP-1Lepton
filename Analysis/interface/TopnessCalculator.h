#ifndef TOPNESSCALCULATOR_H
#define TOPNESSCALCULATOR_H

#include <vector>
#include "MELAParticle.h"
#include "TFitter.h"


namespace TopnessCalculator{
  double topnessFunction(
    double const& pwx_, double const& pwy_, double const& pwz_, double const& pnz_,
    double const& plx_, double const& ply_, double const& plz_, double const& ple_,
    double const& pb1x_, double const& pb1y_, double const& pb1z_, double const& pb1e_,
    double const& pb2x_, double const& pb2y_, double const& pb2z_, double const& pb2e_,
    double const& pmx_, double const& pmy_, double const& pmz_, double const& pme_
  );
  void minuitFcn_Topness(int&, double*, double& result, double* par, int){
    result = topnessFunction(
      par[0], par[1], par[2], par[3],
      par[4], par[5], par[6], par[7],
      par[8], par[9], par[10], par[11],
      par[12], par[13], par[14], par[15],
      par[16], par[17], par[18], par[19]
    );
  }

  double topnessModFunction(
    double const& pwx_, double const& pwy_, double const& pwz_, double const& pnz_,
    double const& plx_, double const& ply_, double const& plz_, double const& ple_,
    double const& pb1x_, double const& pb1y_, double const& pb1z_, double const& pb1e_,
    double const& pb2x_, double const& pb2y_, double const& pb2z_, double const& pb2e_,
    double const& pmx_, double const& pmy_, double const& pmz_, double const& pme_
  );
  void minuitFcn_ModTopness(int&, double*, double& result, double* par, int){
    result = topnessModFunction(
      par[0], par[1], par[2], par[3],
      par[4], par[5], par[6], par[7],
      par[8], par[9], par[10], par[11],
      par[12], par[13], par[14], par[15],
      par[16], par[17], par[18], par[19]
    );
  }

  float topnessMinimization(bool isModified, float const& METx, float const& METy, MELAParticle const* const& lep, MELAParticle const* const& bjet1, MELAParticle const* const& bjet2);
  float CalcTopness(bool isModified, float const& MET, float const& METphi, MELAParticle const* const& lep, std::vector<MELAParticle const*> const& bjets, std::vector<MELAParticle const*> const& addjets);

}


#endif
