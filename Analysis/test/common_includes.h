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
#include "ElectronHandler.h"
#include "ElectronSelectionHelpers.h"
#include "ElectronScaleFactorHandler.h"
#include "MuonHandler.h"
#include "MuonSelectionHelpers.h"
#include "MuonScaleFactorHandler.h"
#include "MELAStreamHelpers.hh"
#include "TSystem.h"
#include "TDirectory.h"
#include "TFile.h"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;
using namespace SampleHelpers;
