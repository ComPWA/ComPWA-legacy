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


  double controlParameter(const vector<double>& minPar);
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
