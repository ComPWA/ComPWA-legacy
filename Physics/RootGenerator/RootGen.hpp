//-------------------------------------------------------------------------------
// Copyright (c) 2013 michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Class to generate 4-Vectors.
/*!
 * @file RootGen.hpp
 * This application uses arbitrary Amplitude implementations and the root
 * phase-space generator to generate a root-file with 4-Vectors distributed
 * according the given amplitude model.
*/

// Standard header files go here
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>
#include <memory>

// Physics Interface header files go here
#include "Physics/Amplitude.hpp"
#include "Core/PWAParameter.hpp"
#include "Core/PWAGenericPar.hpp"

class RootGen{

public:
  RootGen(double inM, std::vector<double> fiM, std::vector<int> pids);

  virtual ~RootGen();

  virtual void generate(unsigned int nEvents, std::shared_pointer<Amplitude> model, std::string filename);

  inline virtual void setBr(const double inBr){initialBr_=inBr;};

private:
  //sys
  double initialM_; // GeV/c²
  double initialBr_; // GeV/c²
  std::vector<double> finalM_ ; // GeV/c²
  std::vector<int> pids_;
};
