#ifndef METOBJECT_H
#define METOBJECT_H


class METVariables{
public:
  float met_raw;
  float phi_raw;
  float met;
  float phi;
  float met_JEC;
  float phi_JEC;
  float met_JECup;
  float phi_JECup;
  float met_JECdn;
  float phi_JECdn;

  METVariables();
  METVariables(METVariables const& other);
  METVariables& operator=(const METVariables& other);

  void swap(METVariables& other);

};

class METObject{
public:
  METVariables extras;

  METObject();
  METObject(const METObject& other);
  METObject& operator=(const METObject& other);
  ~METObject();

  void swap(METObject& other);

};

#endif
