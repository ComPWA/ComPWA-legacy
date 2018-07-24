// Matrix_Base class definition file. -*- C++ -*-
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
#ifndef _Matrix_Base_H
#define _Matrix_Base_H

#include <cassert>
#include "qft++/Conversion.h"
//_____________________________________________________________________________
/** @file Matrix_Base.h
 *  @brief Matrix_Base class definition file (also defines IsMatrix macro).
 */
//_____________________________________________________________________________
namespace ComPWA {
    namespace QFT {

//_____________________________________________________________________________
/** @class Matrix_Base
 *  @author Mike Williams
 *
 *  @brief Base class for all Matrix template classes.
 *
 * This is the base class for all Matrix instantiations. It is thru inheritance
 * from this class that the Matrix template classes are able to correctly
 * implement multiplication.
 */
//_____________________________________________________________________________

class Matrix_Base {

protected:

  // attributes:
  int _num_rows; ///< Number of rows in the Matrix
  int _num_cols; ///< Number of columns in the Matrix

  // protected functions:

  /// Copy @a mbase data members
  void _Copy(const Matrix_Base &__mbase){
    _num_rows = __mbase._num_rows;
    _num_cols = __mbase._num_cols;
  }
  
  /// Assert that @a mbase and @a this have the same size
  void _SizeAssert(const Matrix_Base &__mbase) const {
    assert(this->_num_rows == __mbase._num_rows 
	   && this->_num_cols == __mbase._num_cols);
  }

public:

  // create/copy/destroy:

  /// Default Constructor (0x0)
  Matrix_Base(){ 
    _num_rows = 0;
    _num_cols = 0;
  }

  /** Constructor
   * @param rows Number of rows
   *  @param cols Number of columns
   */
  Matrix_Base(int __rows,int __cols){
    _num_rows = __rows;
    _num_cols = __cols;
  }

  /// Copy Constructor
  Matrix_Base(const Matrix_Base &__mbase){
    this->_Copy(__mbase);
  }

  /** Destructor */
  virtual ~Matrix_Base(){}

  // Getters:

  /** Returns the number of rows */
  int NumRows() const {return _num_rows;}

  /** Returns the number of columns */
  int NumCols() const {return _num_cols;}

  /** Returns the number of elements */ 
  int Size() const {
    return _num_rows * _num_cols;
  }

  // functions:

  /// Is @a mbase the same size as this?
  bool SizeCheck(const Matrix_Base &__mbase) const {
    return (__mbase._num_rows == this->_num_rows 
	    && __mbase._num_cols == this->_num_cols);
  }
  
};
//_____________________________________________________________________________

/// Is T a Matrix? (does it inherit from Matrix_Base?)
#define IsMatrix(__T) IsDerived(__T,Matrix_Base)

//_____________________________________________________________________________

    }
}
#endif /* _Matrix_Base_H */
