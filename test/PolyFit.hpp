//! Test implementation of OIFData.hpp.
/*! \class PolyFit
 * @file PolyFit.hpp
 * This class derives from OIFData.hpp, the data-interface of the optimizers. It
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

#include <boost/shared_ptr.hpp>

#include "TROOT.h"
//#include "qft++/topincludes/relativistic-quantum-mechanics.hh"

#include "OIFData.hpp"
#include "PWAParameter.hpp"

using namespace std;

class TFile;
class TGraph;
class TCanvas;
class TRandom;

class PolyFit : public OIFData {

public:

  // create/copy/destroy:

  ///Constructor 
  PolyFit(double p0, double p1, double p2, double p3, double sigma);


  /** Destructor */
  virtual ~PolyFit();


  double controlParameter(vector<PWAParameter<double> >& minPar);
  void drawGraph(double a, double b, double c, double d);
  // Getters:
 
protected:


private:
  shared_ptr<TFile> _theTFile;
  map <unsigned int, TGraph* > _myGraph;

  vector< double > _xValue;
  vector< double > _yValue;

  double _sigma;

};

#endif
