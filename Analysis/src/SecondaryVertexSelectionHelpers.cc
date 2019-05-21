#include "SecondaryVertexSelectionHelpers.h"
#include <cmath>


bool SecondaryVertexSelectionHelpers::testSoftB(SecondaryVertexObject const& vtx){
  return (
    vtx.nTracks>=softb_ntrks_thr
    &&
    vtx.IP2D<softb_IP2D_upperthr
    &&
    vtx.SIP3D>softb_SIP3D_thr
    &&
    cos(vtx.anglePV)>softb_cosPVangle_thr
    &&
    vtx.pt()<softb_pt_upperthr
    );
}

void SecondaryVertexSelectionHelpers::setSelectionBits(SecondaryVertexObject& vtx){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  if (testSoftB(vtx)) vtx.setSelectionBit(kSoftB);
}
