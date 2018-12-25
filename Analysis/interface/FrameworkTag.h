#ifndef FRAMEWORKTAG_H
#define FRAMEWORKTAG_H

#include <string>
#include "TNumericUtil.hh"


typedef TNumericUtil::triplet<unsigned int> FwkTag_t;


class FrameworkTag{
protected:
  std::string rawtag;
  FwkTag_t parsedtag;

public:
  FrameworkTag();
  FrameworkTag(std::string const& tag);
  FrameworkTag(FrameworkTag const& other);

  void setupParsedTag();

  const char* getRawTag() const{ return rawtag.c_str(); }
  FwkTag_t const& getParsedTag() const{ return parsedtag; }

  bool operator==(const FrameworkTag& other) const;
  bool operator>(const FrameworkTag& other) const;
  bool operator<(const FrameworkTag& other) const;
  bool operator!=(const FrameworkTag& other) const;
  bool operator>=(const FrameworkTag& other) const{ return (*this==other || *this>other); }
  bool operator<=(const FrameworkTag& other) const{ return (*this==other || *this<other); }

};


#endif
