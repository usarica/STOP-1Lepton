#include <cmath>
#include "TauSelectionHelpers.h"


bool TauSelectionHelpers::testVetoId(TauObject const& part){ return part.extras.pfDecayModeFinding; }
bool TauSelectionHelpers::testVetoSelection(TauObject const& part){ return (testVetoId(part) && part.extras.pfIso); }
bool TauSelectionHelpers::testPtEtaSkim(TauObject const& part){ if (testVetoSelection(part)) return (part.pt()>=ptThr_skim_veto && fabs(part.eta())<etaThr_skim_veto); else return false; }

void TauSelectionHelpers::setSelectionBits(TauObject& part){
  static_assert(std::numeric_limits<unsigned long long>::digits >= nSelectionBits);

  if (testVetoId(part)) part.setSelectionBit(kVetoID);
  if (testVetoSelection(part)) part.setSelectionBit(kVetoIDIso);
  if (testPtEtaSkim(part)) part.setSelectionBit(kSkimPtEta);
}
