// Definition file for some special tensors -*- C++ -*-
/* Copyright 2008 Mike Williams (mwill@jlab.org)
 *
 * This file is part of qft++.
 *
 * qft++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * qft++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with qft++.  If not, see <http://www.gnu.org/licenses/>.
 */
// Author: Mike Williams
#ifndef _SpecialTensors_H
#define _SpecialTensors_H
//_____________________________________________________________________________

#include "Tensor.h"

//_____________________________________________________________________________
/** @file SpecialTensors.h
 *  @brief Definition file for some special tensors derived from Tensor.
 */
//_____________________________________________________________________________
namespace ComPWA {
    namespace QFT {

/** @class MetricTensor 
 *  @author Mike Williams
 *
 *  @brief \f$ g_{\mu\nu} \f$ : The metric tensor for Minkowski space.
 *
 * A derived class of Tensor<double>, a MetricTensor object is a 2nd rank
 * Tensor equal to 
 * \f$ \left(\begin{array}{cccc} 
 * 1&0&0&0\\0&-1&0&0\\0&0&-1&0\\0&0&0&-1 \end{array}\right) \f$
 */
//_____________________________________________________________________________

class MetricTensor: public Tensor<double> {

public:

  /// Constructor
  MetricTensor():Tensor<double>(2){
    this->Element(0,0) = 1.0;
    this->Element(1,1) = -1.0;
    this->Element(2,2) = -1.0;
    this->Element(3,3) = -1.0;
  }

  /** Destructor */
  ~MetricTensor(){}
};
//_____________________________________________________________________________

/** @class LeviCivitaTensor
 *  @author Mike Williams
 *
 *  @brief \f$\epsilon_{\mu\nu\rho\sigma}\f$ : Totally anti-symmetric Levi-Civita tensor
 *
 * This class is derived from Tensor<double>.
 */
//_____________________________________________________________________________

class LeviCivitaTensor: public Tensor<double> {
  
public:

  /// Constructor
  LeviCivitaTensor():Tensor<double>(4){
    this->Element(0,1,2,3) = 1.0;
    this->Element(0,1,3,2) = -1.0;
    this->Element(0,2,1,3) = -1.0;
    this->Element(0,2,3,1) = 1.0;
    this->Element(0,3,1,2) = 1.0;
    this->Element(0,3,2,1) = -1.0;
    this->Element(1,0,2,3) = -1.0;
    this->Element(1,0,3,2) = 1.0;
    this->Element(1,2,0,3) = 1.0;
    this->Element(1,2,3,0) = -1.0;
    this->Element(1,3,0,2) = -1.0;
    this->Element(1,3,2,0) = 1.0;
    this->Element(2,0,1,3) = 1.0;
    this->Element(2,0,3,1) = -1.0;
    this->Element(2,1,0,3) = -1.0;
    this->Element(2,1,3,0) = 1.0;
    this->Element(2,3,0,1) = 1.0;
    this->Element(2,3,1,0) = -1.0;
    this->Element(3,0,1,2) = -1.0;
    this->Element(3,0,2,1) = 1.0;
    this->Element(3,1,0,2) = 1.0;
    this->Element(3,1,2,0) = -1.0;
    this->Element(3,2,0,1) = -1.0;
    this->Element(3,2,1,0) = 1.0;
  }

  /** Destructor */
  ~LeviCivitaTensor(){}

};
//_____________________________________________________________________________

}    
}
#endif /* _SpecialTensors_H */
