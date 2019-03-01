#ifndef DUPLICATEEVENTHANDLER_H
#define DUPLICATEEVENTHANDLER_H

#include <set>
#include "FrameworkVariables.hh"


class EventIdentifier{
protected:
  RunNumber_t run;
  Lumisection_t lumi;
  EventNumber_t evid;

public:
  EventIdentifier(RunNumber_t run_, Lumisection_t lumi_, EventNumber_t evid_);
  EventIdentifier(EventIdentifier const& other);
  ~EventIdentifier(){}

  bool operator<(EventIdentifier const& other) const;
  bool operator==(EventIdentifier const& other) const;
};

class DuplicateEventHandler : public std::set<EventIdentifier>{
public:
  DuplicateEventHandler(){}
  ~DuplicateEventHandler(){}

  bool isUnique(EventIdentifier const& other);
  bool isUnique(EventIdentifier const& other) const;
  bool isUnique(RunNumber_t const& run, Lumisection_t const& lumi, EventNumber_t const& evid);
  bool isUnique(RunNumber_t const& run, Lumisection_t const& lumi, EventNumber_t const& evid) const;
};


#endif
