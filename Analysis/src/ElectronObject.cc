#include "Samples.h"
#include "ElectronObject.h"


ElectronObject::ElectronObject() :
  MELAParticle()
{
  setup_extras();
}
ElectronObject::ElectronObject(int id_) :
  MELAParticle(id_)
{
  setup_extras();
}
ElectronObject::ElectronObject(int id_, TLorentzVector p4_) :
  MELAParticle(id_, p4_)
{
  setup_extras();
}
ElectronObject::ElectronObject(const ElectronObject& other) :
  MELAParticle(other),
  extras(other.extras)
{}
ElectronObject& ElectronObject::operator=(const ElectronObject& other){
  ElectronObject tmp(other);
  swap(tmp);
  return *this;
}
ElectronObject::~ElectronObject(){}

void ElectronObject::setup_extras(){
  extras.setNamedVal<float>("etaSC", 0.f); // Supercluster eta
  extras.setNamedVal<float>("rho", 0.f); // evt_fixgridfastjet_all_rho, actually a constant over the event
  extras.setNamedVal<float>("energySC", 0.f); // Supercluster energy, eSC
  extras.setNamedVal<float>("sigmaIEtaIEta_full5x5", 0.f); // full5x5_sigmaIetaIeta
  extras.setNamedVal<float>("dEtaIn", 0.f); // abs(dEtaIn)
  extras.setNamedVal<float>("dPhiIn", 0.f); // abs(dPhiIn)
  extras.setNamedVal<float>("hOverE", 0.f); // hOverE
  extras.setNamedVal<float>("ecalEnergy", 0.f);
  extras.setNamedVal<float>("eOverPIn", 0.f);
  extras.setNamedVal<float>("dxyPV", 0.f); // abs(d0)
  extras.setNamedVal<float>("dzPV", 0.f); // abs(dz)
  extras.setNamedVal<int>("exp_innerlayers", 0); // expectedMissingInnerHits
  extras.setNamedVal<bool>("conv_vtx_flag", false); // Conversion veto
}

float ElectronObject::EinvOverPinv(){
  float const& ecalEnergy = extras.namedfloats["ecalEnergy"];
  float const& EoverP = extras.namedfloats["eOverPIn"];
  return (1.f-EoverP)/ecalEnergy;
}
