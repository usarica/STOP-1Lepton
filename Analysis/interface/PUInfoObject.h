#ifndef PUINFOOBJECT_H
#define PUINFOOBJECT_H


class PUInfoObject{
public:
  int bunchCrossing;
  int nPUVertices;
  float nTrueVertices;

  PUInfoObject();
  PUInfoObject(PUInfoObject const& other);
  ~PUInfoObject(){}

  PUInfoObject& operator=(const PUInfoObject& other);
  void swap(PUInfoObject& other);
};

#endif
