//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Test implementation of ControlParameter.
/*! \class PolyFit
 * @file PolyFit.hpp
 * This class derives from ControlParameter, the data-interface of the optimizers. It
 * represents a set of 1-dim data-points, which are created when instantiating
 * this class using a polynomial and smearing the points with a gausian distri-
 * bution. It also provides a draw function to visualize the data-points.
*/

#ifndef _PolyFit_H
#define _PolyFit_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cassert>
#include <memory>

#include "TROOT.h"

#include "Optimizer/ControlParameter.hpp"
#include "Core/ParameterList.hpp"

class TFile;
class TGraph;
class TCanvas;
class TRandom;

class PolyFit : public ComPWA::Optimizer::ControlParameter {

public:

  // create/copy/destroy:
  static std::shared_ptr<ComPWA::Optimizer::ControlParameter> createInstance(double p0, double p1, double p2, double p3, double sigma);

  /** Destructor */
  virtual ~PolyFit();


  double controlParameter(ComPWA::ParameterList& minPar);
  void drawGraph(double a, double b, double c, double d);
  // Getters:
 
protected:


private:
  std::shared_ptr<TFile> _theTFile;
  std::map <unsigned int, TGraph* > _myGraph;

  std::vector< double > _xValue;
  std::vector< double > _yValue;

  double _sigma;

  PolyFit(double p0, double p1, double p2, double p3, double sigma);

};

#endif
