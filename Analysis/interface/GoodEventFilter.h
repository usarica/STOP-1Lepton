#ifndef GOODEVENTFILTER_H
#define GOODEVENTFILTER_H

#include <vector>
#include "FrameworkVariables.hh"


class GoodRunLumisection{
protected:
  RunNumber_t RunNumber;
  Lumisection_t Lumisection_begin;
  Lumisection_t Lumisection_end;

public:
  GoodRunLumisection() : RunNumber(0), Lumisection_begin(0), Lumisection_end(0){}
  GoodRunLumisection(RunNumber_t RunNumber_, Lumisection_t Lumisection_begin_, Lumisection_t Lumisection_end_) : RunNumber(RunNumber_), Lumisection_begin(Lumisection_begin_), Lumisection_end(Lumisection_end_){}
  GoodRunLumisection(GoodRunLumisection const& other) : RunNumber(other.RunNumber), Lumisection_begin(other.Lumisection_begin), Lumisection_end(other.Lumisection_end){}

  bool testEvent(RunNumber_t const& irun, Lumisection_t const& ilumi) const{ return (irun==RunNumber && ilumi>=Lumisection_begin && ilumi<=Lumisection_end); }

  RunNumber_t const& getRunNumber() const{ return RunNumber; }
  Lumisection_t const& getLumisectionBegin() const{ return Lumisection_begin; }
  Lumisection_t const& getLumisectionEnd() const{ return Lumisection_end; }

};

namespace GoodEventFilter{
  std::vector<GoodRunLumisection> getGoodRunLumisectionList();
  bool testEvent(RunNumber_t const& irun, Lumisection_t const& ilumi);
  void printGoodRunLumisectionList();
}



#endif
