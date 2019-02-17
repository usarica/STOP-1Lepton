#ifndef GENERICEVENTANALYZER_H
#define GENERICEVENTANALYZER_H

#include "FrameworkTreeLooperBase.h"


class GenericEventAnalyzer : public FrameworkTreeLooperBase{
protected:
  bool doWeights;
  bool doGenInfo;
  bool doGenParticles;
  bool doEventFilter;
  bool doVertexPUInfos;
  bool doPFCands;
  bool doMuons;
  bool doElectrons;
  bool doPhotons;
  bool doJetMET;
  bool doWriteSelectionVariables;

  bool runEvent(FrameworkTree* tree, float const& externalWgt, SimpleEntry& product);

public:
  GenericEventAnalyzer();
  GenericEventAnalyzer(FrameworkTree* inTree);
  GenericEventAnalyzer(std::vector<FrameworkTree*> const& inTreeList);
  GenericEventAnalyzer(FrameworkSet const* inTreeSet);

  void setWeightsFlag(bool flag){ doWeights = flag; }
  void setGenInfoFlag(bool flag){ doGenInfo = flag; }
  void setGenParticlesFlag(bool flag){ doGenParticles = flag; if (doGenParticles && !doGenInfo) doGenInfo=true; }
  void setEventFilterFlag(bool flag){ doEventFilter = flag; }
  void setVertexPUInfoFlag(bool flag){ doVertexPUInfos = flag; }
  void setPFCandsFlag(bool flag){ doPFCands = flag; }
  void setMuonsFlag(bool flag){ doMuons = flag; }
  void setElectronsFlag(bool flag){ doElectrons = flag; }
  void setPhotonsFlag(bool flag){ doPhotons = flag; }
  void setJetMETFlag(bool flag){ doJetMET = flag; }
  void setWriteSelectionVariables(bool flag){ doWriteSelectionVariables = flag; }

};


#endif
