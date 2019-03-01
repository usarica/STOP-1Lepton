#include "DuplicateEventHandler.h"


using namespace std;


EventIdentifier::EventIdentifier(RunNumber_t run_, Lumisection_t lumi_, EventNumber_t evid_) : run(run_), lumi(lumi_), evid(evid_){}
EventIdentifier::EventIdentifier(EventIdentifier const& other) : run(other.run), lumi(other.lumi), evid(other.evid){}

bool EventIdentifier::operator<(EventIdentifier const& other) const{
  if (run != other.run) return run < other.run;
  else if (lumi != other.lumi) return lumi < other.lumi;
  else if (evid != other.evid) return evid < other.evid;
  else return false;
}
bool EventIdentifier::operator==(EventIdentifier const& other) const{ return (run == other.run && lumi == other.lumi && evid == other.evid); }

bool DuplicateEventHandler::isUnique(EventIdentifier const& other){ return (this->insert(other).second); }
bool DuplicateEventHandler::isUnique(EventIdentifier const& other) const{ return (this->find(other)==this->cend()); }
bool DuplicateEventHandler::isUnique(RunNumber_t const& run, Lumisection_t const& lumi, EventNumber_t const& evid){ return this->isUnique(EventIdentifier(run, lumi, evid)); }
bool DuplicateEventHandler::isUnique(RunNumber_t const& run, Lumisection_t const& lumi, EventNumber_t const& evid) const{ return this->isUnique(EventIdentifier(run, lumi, evid)); }
