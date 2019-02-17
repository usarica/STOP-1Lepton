#ifndef GENPARTICLEOBJECT_H
#define GENPARTICLEOBJECT_H

#include "ParticleObject.h"


struct GenEventInfo{
  unsigned int processID;
  float qscale;
  float alphaS;
  float xsec;
  float genMET;
  float genMETPhi;
};

class GenParticleVariables{
public:
  bool isPromptFinalState;
  bool isPromptDecayed;
  bool isDirectPromptTauDecayProductFinalState;
  bool isHardProcess;
  bool fromHardProcessFinalState;
  bool fromHardProcessDecayed;
  bool isDirectHardProcessTauDecayProductFinalState;
  bool fromHardProcessBeforeFSR;
  bool isLastCopy;
  bool isLastCopyBeforeFSR;
  int status;

  GenParticleVariables();
  GenParticleVariables(GenParticleVariables const& other);
  GenParticleVariables& operator=(const GenParticleVariables& other);

  void swap(GenParticleVariables& other);

};

class GenParticleObject : public ParticleObject{
public:
  GenParticleVariables extras;

  GenParticleObject();
  GenParticleObject(int id_);
  GenParticleObject(int id_, CMSLorentzVector const& mom_);
  GenParticleObject(const GenParticleObject& other);
  GenParticleObject& operator=(const GenParticleObject& other);
  ~GenParticleObject();

  void swap(GenParticleObject& other);

};

#endif
