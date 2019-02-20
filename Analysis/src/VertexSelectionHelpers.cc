#include "VertexSelectionHelpers.h"


bool VertexSelectionHelpers::testGoodVertex(VertexObject const& vtx){
  return (vtx.isValid && !vtx.isFake && vtx.ndof>=vtx_ndof_thr && vtx.rho()<=vtx_rho_thr && vtx.z()<=vtx_z_thr);
}

void VertexSelectionHelpers::setSelectionBits(VertexObject& vtx){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  if (testGoodVertex(vtx)) vtx.setSelectionBit(kGoodVertex);
}
