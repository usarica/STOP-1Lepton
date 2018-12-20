#ifndef SAMPLES_H
#define SAMPLES_H

#include "HostHelpersCore.h"
#include <string>

// Package directory
#ifndef xstr_lit
#define xstr_lit(s) str_lit(s)
#define str_lit(s) #s
#endif
#ifndef _stop1lpkgpathstr_
#ifndef _stop1lpkgpath_
#define _stop1lpkgpath_ ./
#endif
#define _stop1lpkgpathstr_ xstr_lit(_stop1lpkgpath_)
#endif
const TString STOP1LPKGPATH = _stop1lpkgpathstr_;
const TString STOP1LPKGDATAPATH = STOP1LPKGPATH + "/data/";
const TString CMSTASCOREPKGPATH = STOP1LPKGPATH + "../../cmstas/CORE/";

//constexpr float xsecScale = 1e3;
constexpr float xsecScale = 1;

// LHC sqrts and data period
constexpr unsigned int theSqrts = 13;
extern TString theDataPeriod;

// CMS4 trees
const TString CMS4_EVENTS_TREE_NAME = "Events";
const TString CMS4_RUNS_TREE_NAME = "Runs";

namespace SampleHelpers{
  void setDataPeriod(const TString theDataPeriod_);
}

#endif
