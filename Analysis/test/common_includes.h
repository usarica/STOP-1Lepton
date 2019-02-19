#include <iostream>
#include <cmath>
#include <cassert>
#include "FrameworkVariables.hh"
#include "Samples.h"
#include "SampleHelpers.h"
#include "HelperFunctions.h"
#include "FrameworkTag.h"
#include "FrameworkOptionParser.h"
#include "FrameworkTree.h"
#include "FrameworkSet.h"
#include "FrameworkTreeLooperBase.h"
#include "WeightsHandler.h"
#include "GenInfoHandler.h"
#include "EventFilterHandler.h"
#include "VertexPUHandler.h"
#include "PFCandHandler.h"
#include "MuonHandler.h"
#include "ElectronHandler.h"
#include "PhotonHandler.h"
#include "JetMETHandler.h"
#include "IsoTrackHandler.h"
#include "TauHandler.h"
#include "MuonScaleFactorHandler.h"
#include "ElectronScaleFactorHandler.h"
#include "PhotonScaleFactorHandler.h"
#include "PUScaleFactorHandler.h"
#include "MuonSelectionHelpers.h"
#include "ElectronSelectionHelpers.h"
#include "PhotonSelectionHelpers.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "IsoTrackSelectionHelpers.h"
#include "TauSelectionHelpers.h"
#include "VertexSelectionHelpers.h"
#include "TopnessCalculator.h"
#include "MELAStreamHelpers.hh"
#include "TSystem.h"
#include "TDirectory.h"
#include "TFile.h"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;
using namespace SampleHelpers;
