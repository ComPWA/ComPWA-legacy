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
//! Estimator for a Fit of Dalitz-Plot Slices.
/*! \class SliceFit
 * @file SliceFit.hpp
 * This class performs a Chi2-Fit on slices along one axis of a dalitz-plot.
 * The Dalitz-Plot is generated directly in the constructor of this Estimator.
 * Data and Model are provided in the constructor using the Amplitude and Data
 * interfaces. The class itself fulfills the Estimator interface.
*/

#ifndef _SLICEFIT_HPP
#define _SLICEFIT_HPP

#include <vector>
#include <memory>
#include <string>

//Root Header
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"

//PWA-Header
#include "Estimator/Estimator.hpp"
#include "Physics/Amplitude.hpp"
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"

class SliceFit : public Estimator {

public:
  /// Default Constructor (0x0)
  SliceFit(std::shared_ptr<Amplitude>, std::shared_ptr<Data>);

  virtual double controlParameter(ParameterList& minPar);

  /** Destructor */
  virtual ~SliceFit();

protected:

private:
  TH2D* dalitzPlot_;

  std::shared_ptr<Amplitude> pPIF_;
  std::shared_ptr<Data> pDIF_;

  /*double M = 3.096916; // GeV/c² (J/psi+)
  double Br = 0.000093; // GeV/c² (width)
  double m1 = 0.; // GeV/c² (gamma)
  double m2 = 0.139570; // GeV/c² (pi)
  double m3 = 0.139570; // GeV/c² (pi)
  double PI = 3.14159; // m/s*/

};

#endif /* _SLICEFIT_HPP */
