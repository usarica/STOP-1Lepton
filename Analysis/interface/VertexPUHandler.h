#ifndef VERTEXPUHANDLER_H
#define VERTEXPUHANDLER_H

#include <vector>
#include "IvyBase.h"
#include "VertexObject.h"
#include "SecondaryVertexObject.h"
#include "PUInfoObject.h"


class VertexPUHandler : public IvyBase{
protected:
  bool doVertices;
  bool doSecondaryVertices;
  bool doPUInfos;

  std::vector<VertexObject*> vertices;
  std::vector<SecondaryVertexObject*> secondary_vertices;
  std::vector<PUInfoObject*> puinfos;

  void clear(){ for (VertexObject*& vtx:vertices) delete vtx; vertices.clear(); for (SecondaryVertexObject*& vtx:secondary_vertices) delete vtx; secondary_vertices.clear(); for (PUInfoObject*& puinfo:puinfos) delete puinfo; puinfos.clear(); }

  bool constructVertices();
  bool constructSecondaryVertices();
  bool constructPUInfos();

public:
  // Constructors
  VertexPUHandler();

  // Destructors
  ~VertexPUHandler(){ clear(); }

  void setVerticesFlag(bool flag){ doVertices=flag; }
  void setSecondaryVerticesFlag(bool flag){ doSecondaryVertices=flag; }
  void setPUInfosFlag(bool flag){ doPUInfos=flag; }

  bool getVerticesFlag() const{ return doVertices; }
  bool getSecondaryVerticesFlag() const{ return doSecondaryVertices; }
  bool getPUInfosFlag() const{ return doPUInfos; }

  bool constructVertexPUInfos();

  std::vector<VertexObject*> const& getVertices() const{ return vertices; }
  std::vector<SecondaryVertexObject*> const& getSecondaryVertices() const{ return secondary_vertices; }
  std::vector<PUInfoObject*> const& getPUInfos() const{ return puinfos; }

  void bookBranches(BaseTree* tree);

};


#endif
