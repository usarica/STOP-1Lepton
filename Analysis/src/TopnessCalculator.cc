#include "TopnessCalculator.h"


double TopnessCalculator::topnessFunction(
  double const& pwx_, double const& pwy_, double const& pwz_, double const& pnz_,
  double const& plx_, double const& ply_, double const& plz_, double const& ple_,
  double const& pb1x_, double const& pb1y_, double const& pb1z_, double const& pb1e_,
  double const& pb2x_, double const& pb2y_, double const& pb2z_, double const& pb2e_,
  double const& pmx_, double const& pmy_, double const& pmz_, double const& pme_
){
  constexpr double mW = 81.;
  constexpr double mT = 172.;
  constexpr double aW = 5.;
  constexpr double aT = 15.;
  constexpr double aCM = 1000.;

  // construct the lorentz vectors
  TLorentzVector vW(pwx_, pwy_, pwz_, (sqrt((mW*mW)+(pwx_*pwx_)+(pwy_*pwy_)+(pwz_*pwz_))));
  TLorentzVector vL(plx_, ply_, plz_, ple_);
  TLorentzVector vB1(pb1x_, pb1y_, pb1z_, pb1e_);
  TLorentzVector vB2(pb2x_, pb2y_, pb2z_, pb2e_);
  TLorentzVector vMET(pmx_, pmy_, pmz_, pme_);
  TLorentzVector vN((pmx_-pwx_), (pmy_-pwy_), pnz_, (sqrt(pow((pmx_-pwx_), 2)+pow((pmy_-pwy_), 2)+pow(pnz_, 2))));
  // construct the w-term (lep)
  double tWL = (pow(((mW*mW) - ((vL+vN).M2())), 2)) / (pow(aW, 4));
  // construct the tL-term [seen lepton]
  double tTL = (pow(((mT*mT) - ((vL+vB1+vN).M2())), 2)) / (pow(aT, 4));
  // construct the tM-term [miss lepton]
  double tTM = (pow(((mT*mT) - ((vB2+vW).M2())), 2)) / (pow(aT, 4));
  // construct the CM-term
  double tCM = (pow(((4*(mT*mT)) - ((vL+vN+vW+vB1+vB2).M2())), 2)) / (pow(aCM, 4));
  // calculate Topness
  double Topness = tWL + tTL + tTM + tCM;
  return Topness;
}
void TopnessCalculator::minuitFcn_Topness(int&, double*, double& result, double* par, int){
  result = topnessFunction(
    par[0], par[1], par[2], par[3],
    par[4], par[5], par[6], par[7],
    par[8], par[9], par[10], par[11],
    par[12], par[13], par[14], par[15],
    par[16], par[17], par[18], par[19]
  );
}

double TopnessCalculator::topnessModFunction(
  double const& pwx_, double const& pwy_, double const& pwz_, double const& pnz_,
  double const& plx_, double const& ply_, double const& plz_, double const& ple_,
  double const& pb1x_, double const& pb1y_, double const& pb1z_, double const& pb1e_,
  double const& pb2x_, double const& pb2y_, double const& pb2z_, double const& pb2e_,
  double const& pmx_, double const& pmy_, double const& pmz_, double const& pme_
){
  constexpr double mW = 81.;
  constexpr double mT = 172.;
  constexpr double aW = 5.;
  constexpr double aT = 15.;

  // construct the lorentz vectors
  TLorentzVector vW(pwx_, pwy_, pwz_, (sqrt((mW*mW)+(pwx_*pwx_)+(pwy_*pwy_)+(pwz_*pwz_))));
  TLorentzVector vL(plx_, ply_, plz_, ple_);
  TLorentzVector vB1(pb1x_, pb1y_, pb1z_, pb1e_);
  TLorentzVector vB2(pb2x_, pb2y_, pb2z_, pb2e_);
  TLorentzVector vMET(pmx_, pmy_, pmz_, pme_);
  TLorentzVector vN((pmx_-pwx_), (pmy_-pwy_), pnz_, (sqrt(pow((pmx_-pwx_), 2)+pow((pmy_-pwy_), 2)+pow(pnz_, 2))));
  // construct the w-term (lep)
  double tWL = (pow(((mW*mW) - ((vL+vN).M2())), 2)) / (pow(aW, 4));
  // construct the tM-term [miss lepton]
  double tTM = (pow(((mT*mT) - ((vB2+vW).M2())), 2)) / (pow(aT, 4));
  // calculate Topness
  double Topness = tWL + tTM;
  return Topness;
}
void TopnessCalculator::minuitFcn_ModTopness(int&, double*, double& result, double* par, int){
  result = topnessModFunction(
    par[0], par[1], par[2], par[3],
    par[4], par[5], par[6], par[7],
    par[8], par[9], par[10], par[11],
    par[12], par[13], par[14], par[15],
    par[16], par[17], par[18], par[19]
  );
}

float TopnessCalculator::topnessMinimization(bool isModified, float const& METx, float const& METy, MELAParticle const* const& lep, MELAParticle const* const& bjet1, MELAParticle const* const& bjet2){
  TFitter minimizer(4);
  double p1 = -1;
  minimizer.ExecuteCommand("SET PRINTOUT", &p1, 1);
  if (isModified) minimizer.SetFCN(minuitFcn_ModTopness);
  else minimizer.SetFCN(minuitFcn_Topness);
  // get variables for Topness
  float iLpx = lep->x();
  float iLpy = lep->y();
  float iLpz = lep->z();
  float iLpe = lep->t();
  float iB1px = bjet1->x();
  float iB1py = bjet1->y();
  float iB1pz = bjet1->z();
  float iB1pe = bjet1->t();
  float iB2px = bjet2->x();
  float iB2py = bjet2->y();
  float iB2pz = bjet2->z();
  float iB2pe = bjet2->t();
  float iMpx = METx;
  float iMpy = METy;
  float iMpz = 0;
  float iMpe = sqrt(pow(METx, 2)+pow(METy, 2));

  // Define parameters [param number, param name, init val, estimated distance to min, bla, bla] // 300,3000,-3000,3000
  minimizer.SetParameter(0, "pwx", 0, 500, -3000, 3000);
  minimizer.SetParameter(1, "pwy", 0, 500, -3000, 3000);
  minimizer.SetParameter(2, "pwz", 0, 500, -3000, 3000);
  minimizer.SetParameter(3, "pnz", 0, 500, -3000, 3000);
  // fixed parameters
  minimizer.SetParameter(4, "plx", iLpx, 0, iLpx-0.001, iLpx+0.001);
  minimizer.SetParameter(5, "ply", iLpy, 0, iLpy-0.001, iLpy+0.001);
  minimizer.SetParameter(6, "plz", iLpz, 0, iLpz-0.001, iLpz+0.001);
  minimizer.SetParameter(7, "ple", iLpe, 0, iLpe-0.001, iLpe+0.001);
  minimizer.SetParameter(8, "pb1x", iB1px, 0, iB1px-0.001, iB1px+0.001);
  minimizer.SetParameter(9, "pb1y", iB1py, 0, iB1py-0.001, iB1py+0.001);
  minimizer.SetParameter(10, "pb1z", iB1pz, 0, iB1pz-0.001, iB1pz+0.001);
  minimizer.SetParameter(11, "pb1e", iB1pe, 0, iB1pe-0.001, iB1pe+0.001);
  minimizer.SetParameter(12, "pb2x", iB2px, 0, iB2px-0.001, iB2px+0.001);
  minimizer.SetParameter(13, "pb2y", iB2py, 0, iB2py-0.001, iB2py+0.001);
  minimizer.SetParameter(14, "pb2z", iB2pz, 0, iB2pz-0.001, iB2pz+0.001);
  minimizer.SetParameter(15, "pb2e", iB2pe, 0, iB2pe-0.001, iB2pe+0.001);
  minimizer.SetParameter(16, "pmx", iMpx, 0, iMpx-0.001, iMpx+0.001);
  minimizer.SetParameter(17, "pmy", iMpy, 0, iMpy-0.001, iMpy+0.001);
  minimizer.SetParameter(18, "pmz", iMpz, 0, iMpz-0.001, iMpz+0.001);
  minimizer.SetParameter(19, "pme", iMpe, 0, iMpe-0.001, iMpe+0.001);
  minimizer.FixParameter(4);
  minimizer.FixParameter(5);
  minimizer.FixParameter(6);
  minimizer.FixParameter(7);
  minimizer.FixParameter(8);
  minimizer.FixParameter(9);
  minimizer.FixParameter(10);
  minimizer.FixParameter(11);
  minimizer.FixParameter(12);
  minimizer.FixParameter(13);
  minimizer.FixParameter(14);
  minimizer.FixParameter(15);
  minimizer.FixParameter(16);
  minimizer.FixParameter(17);
  minimizer.FixParameter(18);
  minimizer.FixParameter(19);

  minimizer.ExecuteCommand("SIMPLEX", 0, 0);
  //Get the best fit values
  float pwx_fit = minimizer.GetParameter(0);
  float pwy_fit = minimizer.GetParameter(1);
  float pwz_fit = minimizer.GetParameter(2);
  float pnz_fit = minimizer.GetParameter(3);

  // get the function value at best fit
  if (isModified) return topnessModFunction(
    pwx_fit, pwy_fit, pwz_fit, pnz_fit,
    iLpx, iLpy, iLpz, iLpe,
    iB1px, iB1py, iB1pz, iB1pe,
    iB2px, iB2py, iB2pz, iB2pe,
    iMpx, iMpy, iMpz, iMpe
  );
  else return topnessFunction(
    pwx_fit, pwy_fit, pwz_fit, pnz_fit,
    iLpx, iLpy, iLpz, iLpe,
    iB1px, iB1py, iB1pz, iB1pe,
    iB2px, iB2py, iB2pz, iB2pe,
    iMpx, iMpy, iMpz, iMpe
  );
}
float TopnessCalculator::CalcTopness(bool isModified, float const& MET, float const& METphi, MELAParticle const* const& lep, std::vector<MELAParticle const*> const& bjets, std::vector<MELAParticle const*> const& addjets){
  float metx = MET*TMath::Cos(METphi);
  float mety = MET*TMath::Sin(METphi);

  float topness = 1e9;
  if ((bjets.size()+addjets.size())<2) return -1;
  std::vector<MELAParticle const*> alljets = bjets; for (MELAParticle const* jet:addjets) alljets.push_back(jet);

  for (auto it_b1=alljets.cbegin(); it_b1!=alljets.cend(); it_b1++){
    auto const& b1=*it_b1;
    for (auto it_b2=it_b1; it_b2!=alljets.cend(); it_b2++){
      auto const& b2=*it_b2;
      float tmptop = std::min(topnessMinimization(isModified, metx, mety, lep, b1, b2), topnessMinimization(isModified, metx, mety, lep, b2, b1));
      if (tmptop<topness) topness = tmptop;
    }
  }
  return topness;
}
