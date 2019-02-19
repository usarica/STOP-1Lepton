#include <cmath>
#include "IsoTrackSelectionHelpers.h"


bool IsoTrackSelectionHelpers::testVetoId(IsoTrackObject const& part){
  auto const& extras = part.extras;
  return (extras.isPFCand && !extras.hasLepOverlap && extras.charge!=0 && fabs(extras.dz)<0.1);
}
bool IsoTrackSelectionHelpers::testVetoSelection(IsoTrackObject const& part){
  if (!testVetoId(part)) return false;

  float const& pfiso_ch = part.extras.pfIso_ch;
  float const pfiso_ch_thr = std::min(part.pt()*0.1f, 6.f);
  return (pfiso_ch < pfiso_ch_thr);
}
bool IsoTrackSelectionHelpers::testPtEtaSkim(IsoTrackObject const& part){
  return (testVetoSelection(part) && part.pt()>=ptThr_skim_veto && fabs(part.eta())<etaThr_skim_veto);
}

void IsoTrackSelectionHelpers::setSelectionBits(IsoTrackObject& part){
  if (testVetoId(part)) part.setSelectionBit(kVetoID);
  if (testVetoSelection(part)) part.setSelectionBit(kVetoIDIso);
  if (testPtEtaSkim(part)) part.setSelectionBit(kSkimPtEta);
}
