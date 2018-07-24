// Tensor_Base class definition file -*- C++ -*-
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
#ifndef _Tensor_Base_H
#define _Tensor_Base_H

#include "qft++/Conversion.h"
#include <cassert>
//_____________________________________________________________________________
/** @file Tensor_Base.h
 *  @brief Tensor_Base class definition file.
 */
//_____________________________________________________________________________
namespace ComPWA {
    namespace QFT {

//_____________________________________________________________________________
/** @class Tensor_Base
 *  @author Mike Williams
 *
 *  @brief Base class for all Tensor template classes.
 *
 * This is the base class for all Tensor instantiations. It is thru inheritance
 * from this class that the Tensor template classes are able to correctly
 * implement multiplication/contraction.
 */
//_____________________________________________________________________________

class Tensor_Base {

protected:
  
  int _rank; ///< Rank of the Tensor

  // protected functions:

  /// Copy @a tbase data members
  void _Copy(const Tensor_Base &__tbase){
    _rank = __tbase._rank;
  }

  /// Assert that @a tbase and @a this have the same rank
  inline void RankAssert(const Tensor_Base &__tbase) {
    assert(this->_rank == __tbase._rank);
  }
  //inline bool RankAssert(const Tensor_Base &__tbase) {
    //assert(this->_rank == __tbase._rank);
  //}

public:

  // create/copy/destroy:

  /** Default Constructor */
  Tensor_Base(){_rank = 0;}

  /** Constructor */
  Tensor_Base(int __rank){_rank = __rank;}

  /**Copy Ctor*/ 
  Tensor_Base(const Tensor_Base &__tbase){this->_Copy(__tbase);}

  /** Destructor */
  virtual ~Tensor_Base(){}

  // Getters:

  /** Returns the rank of the Tensor */ 
  inline int Rank() const {return _rank;}

  /// Is @a tbase the same rank as @a this?
  inline bool RankCheck(const Tensor_Base &__tbase) const {
    return (this->_rank == __tbase._rank);
  }

};
//_____________________________________________________________________________

/// Is T a Tensor? (does it inherit from Tensor_Base?)
#define IsTensor(__T) IsDerived(__T,Tensor_Base)

#include "qft++/Matrix_Base.h"
/// Is T a scalar? (is it NOT a Tensor or Matrix)
#define IsScalar(__T) (!IsMatrix(__T) && !IsTensor(__T))

/// Is T a scalar that should be passed by value?
#define IsScalarValType(__T) (IsScalar(__T) && PassByValue(__T))

/// Is T a scalar that should be passed by reference?
#define IsScalarRefType(__T) (IsScalar(__T) && !PassByValue(__T))

//_____________________________________________________________________________

    }
    
}
#endif /* _Tensor_Base_H */
