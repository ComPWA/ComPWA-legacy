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
//#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "Optimizer/ControlParameter.hpp"
#include "Core/ParameterList.hpp"

class TFile;
class TGraph;
class TCanvas;
class TRandom;

class PolyFit : public ControlParameter {

public:

  // create/copy/destroy:
  static std::shared_ptr<ControlParameter> createInstance(double p0, double p1, double p2, double p3, double sigma);

  /** Destructor */
  virtual ~PolyFit();


  double controlParameter(ParameterList& minPar);
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
