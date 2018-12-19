#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include "Samples.h"
#include "GoodEventFilter.h"
#include "MELAStreamHelpers.hh"


namespace GoodEventFilter{
  std::vector<GoodRunLumisection> const list_GoodRunLumisection = getGoodRunLumisectionList();
}

using namespace std;
using namespace MELAStreamHelpers;


std::vector<GoodRunLumisection> GoodEventFilter::getGoodRunLumisectionList(){
  std::vector<GoodRunLumisection> res;
  RunNumber_t irun;
  Lumisection_t ilumi, jlumi;

  const TString filename = STOP1LPKGDATAPATH + "Good_Runs_Lumisections.txt";
  ifstream fin; fin.open(filename.Data());
  if (fin.good()){
    while (!fin.eof()){
      fin >> irun >> ilumi >> jlumi;
      res.emplace_back(irun, ilumi, jlumi);
    }
  }

  return res;
}
bool GoodEventFilter::testEvent(RunNumber_t const& irun, Lumisection_t const& ilumi){
  for (GoodRunLumisection const& grl:list_GoodRunLumisection){
    if (grl.testEvent(irun, ilumi)) return true;
  }
  return false;
}
void GoodEventFilter::printGoodRunLumisectionList(){
  for (GoodRunLumisection const& grl:list_GoodRunLumisection) MELAout
    << "Run " << grl.getRunNumber()
    << ", lumi. sections [" << grl.getLumisectionBegin() << ", " << grl.getLumisectionEnd() << "]"
    << endl;
}
