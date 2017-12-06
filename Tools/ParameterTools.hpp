// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Collection of useful routines to modify ParameterLists
///

#include <vector>
#include <memory>
#include <cassert>

#include "Core/ParameterList.hpp"

#ifndef ParameterTools_h
#define ParameterTools_h

#include "TRandom3.h"

///
/// Expands '~' in file path to user directory.
///
inline std::string expand_user(std::string p) {
  std::string path = p;
  if (not path.empty() and path[0] == '~') {
	  assert(path.size() == 1 or path[1] == '/'); // or other error handling
    char const *home = getenv("HOME");
    if (home || home == getenv("USERPROFILE")) {
      path.replace(0, 1, home);
    } else {
      char const *hdrive = getenv("HOMEDRIVE"), *hpath = getenv("HOMEPATH");
      assert(hdrive); // or other error handling
      assert(hpath);
      path.replace(0, 1, std::string(hdrive) + hpath);
    }
  }
  return path;
}

///
/// Loop over ParameterList and set the error for each parameter to \p error.
///
inline void setErrorOnParameterList(ComPWA::ParameterList &list, double error,
                                    bool asym) {
  for (unsigned int i = 0; i < list.GetNDouble(); i++) {
    std::shared_ptr<ComPWA::DoubleParameter> p = list.GetDoubleParameter(i);
    if (p->IsFixed()) {
      p->SetError(0.0);
      continue;
    }
    if (asym)
      list.GetDoubleParameter(i)->SetError(error, error);
    else
      list.GetDoubleParameter(i)->SetError(error);
  }
}

///
/// Loop over ParameterList and set each parameter to a random value. The
/// random value is uniformly choosen within the parameter bounds
///
inline void randomStartValues(ComPWA::ParameterList &fitPar) {
  for (unsigned int i = 0; i < fitPar.GetNDouble(); i++) {
    std::shared_ptr<ComPWA::DoubleParameter> p = fitPar.GetDoubleParameter(i);
    if (p->IsFixed())
      continue;
    double min = -999, max = 999;
    if (p->HasBounds()) {
      min = p->GetMinValue();
      max = p->GetMaxValue();
    }
    p->SetValue(gRandom->Uniform(min, max));
  }
  std::cout << "Randomizing parameter list. New list:" << fitPar << std::endl;
  return;
}


#endif /* ParameterTools_h */
