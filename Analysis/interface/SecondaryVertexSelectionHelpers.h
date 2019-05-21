#ifndef SECONDARYVERTEXSELECTIONHELPERS_H
#define SECONDARYVERTEXSELECTIONHELPERS_H


#include "SecondaryVertexObject.h"


namespace SecondaryVertexSelectionHelpers{
  constexpr unsigned int softb_ntrks_thr = 4; // >=4
  constexpr float softb_pt_upperthr = 20.; // <20
  constexpr float softb_IP2D_upperthr = 3.; // < 3
  constexpr float softb_SIP3D_thr = 4.; // >4
  constexpr float softb_cosPVangle_thr = 0.98; // >0.98

  enum SelectionBits{
    kSoftB,

    nSelectionBits
  };

  bool testSoftB(SecondaryVertexObject const& vtx);

  void setSelectionBits(SecondaryVertexObject& vtx);

}


#endif
