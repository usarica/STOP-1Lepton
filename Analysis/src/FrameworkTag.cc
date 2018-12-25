#include <cassert>
#include <stdexcept>
#include "FrameworkTag.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


FrameworkTag::FrameworkTag() :
  rawtag(""),
  parsedtag(0, 0, 0)
{}
FrameworkTag::FrameworkTag(std::string const& tag) :
  rawtag(tag),
  parsedtag(0, 0, 0)
{
  setupParsedTag();
}
FrameworkTag::FrameworkTag(FrameworkTag const& other) :
  rawtag(other.rawtag),
  parsedtag(other.parsedtag)
{}

void FrameworkTag::setupParsedTag(){
  // Framework tags look like 'CMS4_V00-00-01'
  if (rawtag=="") return;
  std::vector<std::string> splittag;
  splitOptionRecursive(rawtag, splittag, '_');
  if (splittag.size()<2){
    MELAerr << "FrameworkTag::setupParsedTag: Raw tag " << rawtag << " cannot be split by \'_\' properly. Please modify FrameworkTag::setupParsedTag to account for this case." << endl;
    assert(0);
    return;
  }
  string splitfurther=splittag.at(1); splittag.clear();
  replaceString<std::string, const char*>(splitfurther, "V", "");
  splitOptionRecursive(splitfurther, splittag, '-', false);
  if (splittag.size()!=3){
    MELAerr << "FrameworkTag::setupParsedTag: Cannot remove \'V\' properly from " << splitfurther << " and split it into 3 components (raw tag: " << rawtag << "). Please modify FrameworkTag::setupParsedTag to account for this case." << endl;
    MELAerr << "- The splittag vector is " << splittag << endl;
    assert(0);
    return;
  }
  vector<unsigned int> tagvars; tagvars.reserve(3);
  for (std::string const& s:splittag){
    unsigned int is;
    try{ is = stoi(s.c_str()); }
    catch (std::invalid_argument& e){
      MELAerr << "FrameworkTag::setupParsedTag: Error in casting " << s << " to an integer!" << endl;
      assert(0);
    }
    tagvars.push_back(is);
  }
  parsedtag = FwkTag_t(tagvars.at(0), tagvars.at(1), tagvars.at(2));
}

bool FrameworkTag::operator==(const FrameworkTag& other) const{ return (parsedtag==other.parsedtag); }
bool FrameworkTag::operator!=(const FrameworkTag& other) const{ return (parsedtag!=other.parsedtag); }
bool FrameworkTag::operator>(const FrameworkTag& other) const{
  if (parsedtag[0]>other.parsedtag[0]) return true;
  else if (parsedtag[0]==other.parsedtag[0] && parsedtag[1]>other.parsedtag[1]) return true;
  else if (parsedtag[0]==other.parsedtag[0] && parsedtag[1]==other.parsedtag[1] && parsedtag[2]>other.parsedtag[2]) return true;
  else return false;
}
bool FrameworkTag::operator<(const FrameworkTag& other) const{
  if (parsedtag[0]<other.parsedtag[0]) return true;
  else if (parsedtag[0]==other.parsedtag[0] && parsedtag[1]<other.parsedtag[1]) return true;
  else if (parsedtag[0]==other.parsedtag[0] && parsedtag[1]==other.parsedtag[1] && parsedtag[2]<other.parsedtag[2]) return true;
  else return false;
}
