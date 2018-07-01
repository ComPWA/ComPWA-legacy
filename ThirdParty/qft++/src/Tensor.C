// Tensor template class source file. -*- C++ -*-
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

//_____________________________________________________________________________
/** @file Tensor.tcc
 *  @brief Tensor template class source file.
 */
//_____________________________________________________________________________
#include <cmath>
#include "qft++/Tensor.h"
#include "qft++/SpecialTensors.h"

namespace ComPWA {
    namespace QFT {

template class Tensor<double>;


//_____________________________________________________________________________

template <typename _Tp>
Tensor<_Tp> Tensor<_Tp>::operator>>(int __shift) const {

  Tensor<_Tp> ret(_rank);
  int i,j;

  if(this->_rank > 1){
    TensorIndex index(_rank);
    TensorIndex ind(_rank);
    while(index.IsValid()){ // loop over elements
      for(i = 0; i < ind.Size(); i++){
	j = i - __shift;
	while(j < 0) j += _rank;
	ind.SetIndex(i,index[j]);
      }
      ret(index) = this->Element(ind);
      index++;
    }
  }
  else ret = (*this);
 
  return ret;
}
//_____________________________________________________________________________

template <typename _Tp>
Tensor<_Tp> Tensor<_Tp>::operator<<(int __shift) const {

  Tensor<_Tp> ret(_rank);
  int i,j;
  if(this->_rank > 1){
    TensorIndex index(_rank);
    TensorIndex ind(_rank);
    while(index.IsValid()){
      for(i = 0; i < ind.Size(); i++){
	j = i + __shift;
	while(j >= _rank) j -= _rank;       
	ind.SetIndex(i,index[j]);
      }
      ret(index) = this->Element(ind);
      index++;
    }
  }
  else ret = (*this);
 
  return ret;
}
//_____________________________________________________________________________

template <typename _Tp>
Tensor<_Tp> Tensor<_Tp>::Permute(int __mu,int __nu) const {

  Tensor<_Tp> ret(_rank);

  if((this->_rank > 1)&&(__mu < this->_rank)&&(__nu < this->_rank)){
    TensorIndex index(_rank);
    TensorIndex ind(_rank);
    while(index.IsValid()){
      ind = index;
      ind.SetIndex(__mu,index[__nu]);
      ind.SetIndex(__nu,index[__mu]);
      ret(index) = this->Element(ind);
      index++;
    }
  }
  else ret = (*this);
  
  return ret;
}
//_____________________________________________________________________________

template <typename _Tp>
Tensor<_Tp> Tensor<_Tp>::Order(const TensorIndex &__order) const {

  if((int)__order.Size() != _rank){
    std::cout << "Error! Attempt to reorder tensor indicies w/ incorrect number of"
	 <<" indicies." << std::endl;
  }
  assert((int)__order.Size() == _rank);
  Tensor<_Tp> ret(_rank);

  if(_rank > 0){
    TensorIndex index(_rank);
    TensorIndex ind(_rank);
    while(index.IsValid()){ // loop over elements
      for(int i = 0; i < _rank; i++) ind.SetIndex(i,index[__order[i]]);
      
      ret(index) = this->Element(ind);
      index++;
    }  
  }
  return ret;
}
//_____________________________________________________________________________

template <typename _Tp>
Tensor<_Tp> Tensor<_Tp>::Symmetric() const {

  int nterms = 0;
  Tensor<_Tp> ret(_rank);

  // if rank < 2 just return the tensor
  if(_rank > 1){    
    TensorIndex order(_rank);
    // get the 1st permutation (0,1,2,...,rank-1)
    order.Permute();
    while(order.PermIsValid()){ // loop over all valid permutations
      //      order.Print(cout);
      if(nterms == 0) ret = this->Order(order);      
      else ret += this->Order(order);      
      nterms++;
      order.Permute();
    }
    ret /= nterms;
  }
  else ret = *this;
  
  return ret;
}
//_____________________________________________________________________________

template <typename _Tp>
Tensor<_Tp> Tensor<_Tp>::AntiSymmetric() const {

  int nterms = 0,ind;
  Tensor<_Tp> ret(_rank);
  double sign;

  // if rank < 2 just return the tensor
  if(_rank > 1){    
    TensorIndex order(_rank);
    sign = 1.0;
    ind = 1;
    // get the 1st permuation (0,1,2,...,rank -1)
    order.Permute();
    while(order.PermIsValid()){ // loop over all valid permuations
      if(nterms == 0) ret = (this->Order(order))*sign;      
      else ret += (this->Order(order))*sign;
      
      ind++;
      // TensorIndex::Permute returns the permuations in such a way that the 
      // sign for the terms go +--++--++...
      if(ind == 2){
	sign *= -1.0;
	ind = 0;
      }
      nterms++;
      order.Permute();
    }
    ret /= nterms;
  }
  else ret = *this;
  
  return ret;
}
//_____________________________________________________________________________

//_____________________________________________________________________________


    } // namespace QFT
}
