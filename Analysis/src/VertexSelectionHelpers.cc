#include "VertexSelectionHelpers.h"
#include "HelperFunctions.h"


bool VertexSelectionHelpers::testGoodVertex(VertexObject const& vtx){
  return (vtx.isValid && !vtx.isFake && vtx.ndof>=vtx_ndof_thr && vtx.rho()<=vtx_rho_thr && vtx.z()<=vtx_z_thr);
}

void VertexSelectionHelpers::setSelectionBits(VertexObject& vtx){
  using namespace HelperFunctions;
  if (testGoodVertex(vtx)) set_bit(vtx.selectionBits, kGoodVertex);
}
