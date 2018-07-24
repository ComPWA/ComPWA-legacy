// TensorIndex class definition file -*- C++ -*-
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
// Author: Matt Shepherd
#ifndef _TensorIndex_h
#define _TensorIndex_h
//_____________________________________________________________________________

#include <cassert>
#include <iostream>

//_____________________________________________________________________________
/** @file TensorIndex.h
 *  @brief TensorIndex class definition file.
 */
//_____________________________________________________________________________
namespace ComPWA {
    namespace QFT {

/** @class TensorIndex
 *  @author Mike Williams
 *  @author Matt Shepherd
 *
 *  @brief A class that converts \f$(\mu,\nu,...)\f$ to an index for Tensor.
 *
 * This class is a utility class used by Tensor. Its current implementation is
 * limited to tensors of rank 15 or less (well above current memory 
 * capacities).
 * The elements of Tensor<type> T(rank) can be accessed via           
 * TensorIndex index(rank) by simply calling T(index). 
 *
 * <b> Example Usage </b>
 *
 * <!--
 *      TensorIndex index(this->Rank());              
 *      while(index.IsValid()){                       
 *        (*this)(index) "some operations ...";       
 *        ++index;                                                      
 *      }                                                              
 * -->
 * \include TensorIndex.ex
 *
 * Note: Indicies are t = 0, x = 1, y = 2 and z = 3 
 */
//_____________________________________________________________________________

class TensorIndex {

private:

  unsigned int _rank; ///< Tensor rank this deals with
  int _index; ///< Index to Tensor vector of elements

public:

  // create/copy/destroy

  TensorIndex() : _rank(0),_index(0){/** Default Constructor */}

  /// Constructor (connect with rank @a rank tensor)
  TensorIndex(unsigned int __rank) : _index(0){		
    assert( __rank < 16 );
    _rank = __rank;
  }
    
  virtual ~TensorIndex() {/** Destructor */}

  // operators:

  /// Returns the tensor index value for @a level
  inline int operator[](unsigned int __level) const {
    return ((_index >> (__level << 1 ) ) & 0x3);
  }

  /// Returns the current index
  inline int operator()() const {return _index;}

  /** Prefix step up incrementer.
   *
   * If TensorIndex index(rank) is defined, then ++index steps up index
   * to the next valid tensor index.
   *
   * Example: If TensorIndex index(3) = (0,2,1), then ++index = (0,2,2)). 
   *
   * note: Postfix incrementer, index++, is also defined.
   */
  TensorIndex& operator++(){ _index++; return (*this); }

  /// Postfix step up incrementer (see operator++()).
  inline TensorIndex& operator++(int){
    return (++(*this));
  }

  /** Prefix step down incrementer.
   *
   * If TensorIndex index(rank) is defined, then --index steps down index
   * to the next valid tensor index.
   *
   * Example: For TensorIndex index(3) = (0,2,1), --index = (0,2,0). 
   *
   * note: Postfix incrementer, index--, is also defined.
   */
  TensorIndex& operator--(){ _index--; return (*this);}

  /// Postfix step down incrementer (see operator--()).
  inline TensorIndex& operator--(int){
    return (--(*this));
  }

  // permutation functions:

  /** Obtains the next valid permuation of the current stored index.
   *
   * This function is NOT to be used when the TensorIndex object is being
   * used to get elements from a tensor (the normal use). 
   * This function is meant to be used on a TensorIndex object declared for 
   * the sole purpose of performing permutations on the indicies of a Tensor.
   */
  TensorIndex& Permute();

  /** Checks to see if the current index permutation is valid.
   *
   * Valid means all indicies are between 0 and Rank-1. The function doesn't 
   * actually check if any 2 indicies are equal, this is done in 
   * TensorIndex::Permutation.
   */
  bool PermIsValid() const {
    for(int i = 0; i < this->Size(); i++){
      if(((*this)[i] < 0)||((*this)[i] >= this->Size())) return false;
    }
    return true;
  }

  // functions:

  /// Returns the number of entries (coresponds to tensor rank)
  inline int Size() const {return _rank;}

  /// Adjust the rank
  inline void Resize(unsigned int __rank) {     
    assert(__rank < 16 ); 
    _rank = __rank; 
  }

  /// Set the index value
  inline void SetIndex(unsigned int __i,unsigned int __value){
    //assert(__value < 4 && __i < _rank );
    _index = (_index & ~(0x3 << (__i << 1))) + (__value << (__i << 1));
  }

  /// Checks to see if the current tensor index stored is a valid index.
  inline bool IsValid() const {
    if( ( _index >> ( _rank << 1 ) ) != 0 ) return false;
    return true;
  }

  /// Resets all indicies in this object to 0.
  inline void Reset() {_index = 0;}

  /** Print to screen the current set of indicies.
   *  Format is: \f$ (\mu,\nu,\rho,...,\pi) \f$ 
   */
  void Print(std::ostream &__os = std::cout);
};
//_____________________________________________________________________________

    }
    
}
#endif /* _TensorIndex_H */
