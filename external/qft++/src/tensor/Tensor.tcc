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
#ifndef _Tensor_TCC
#define _Tensor_TCC
//_____________________________________________________________________________
/** @file Tensor.tcc
 *  @brief Tensor template class source file.
 */
//_____________________________________________________________________________

template <typename _Tp> template<typename T> 
Tensor<typename MultType<_Tp,T>::Type>
Tensor<_Tp>::Contract(const Tensor<T> &__tensor,int __num_indicies) const {

  int ind1max,ind2st,sumsize,nterm,rank;
  Tensor<typename MultType<_Tp,T>::Type> ret;
  MetricTensor g;
  double gFactors;
  typename MultType<_Tp,T>::Type element;
  ind1max = 0;
  ind2st = 0;
  sumsize = 0;

  if(_rank == 0){ // if this is rank 0, just multiply __tensor by this' value
    ret.SetRank(__tensor._rank);
    ret = _data[0] * __tensor;
    return ret;
  }
  else if(__tensor._rank == 0){//__tensor is rank 0...
    ret.SetRank(_rank);
    ret = (*this) * __tensor._data[0];
    return ret;
  }

  if(__num_indicies > _rank || __num_indicies > __tensor._rank){
    cout << "<Tensor::Contract> Error! Can't contract " << __num_indicies
	 << " between a " << _rank << " rank and a " << __tensor._rank
	 << " rank tensor." << endl;
    abort();
  }
  
  if(__num_indicies > 0) rank = _rank + __tensor._rank - 2*__num_indicies;
  else rank = abs(_rank - __tensor._rank);

  // set ret's rank and create the summed TensorIndex
  ret.SetRank(rank);
  TensorIndex indSummed;

  if(__num_indicies < 0){
    // check to see which tensor has the smaller rank (calculate how many 
    // summed indicies are needed)
    if(_rank <= __tensor._rank) sumsize = _rank;
    else sumsize = __tensor._rank;
  }
  else sumsize = __num_indicies;
  indSummed.Resize(2*sumsize);

  TensorIndex ind1(this->_rank);
  TensorIndex ind2(__tensor._rank);

  int size1 = ind1.Size();
  int size2 = ind2.Size();

  if(rank > 0){ // rank > 0, we need to do all the loops
    TensorIndex index(rank);

    while(index.IsValid()){ // loop over ret's elements

      // check to see if this will have any free indicies
      if((size1 - sumsize) > 0) ind1max = size1 - sumsize;
      else ind1max = 0;
       
      // set this index (except last ??(number of summed indicies) indicies)
      for(int i = 0; i < ind1max; i++) ind1.SetIndex(i,index[i]);

      // set __tensor index (except 1st ??(summed indicies) indicies)
      ind2st = sumsize;
      for(int i = ind2st; i < size2; i++) 
	ind2.SetIndex(i,index[ind1max+(i-ind2st)]);

      nterm = 0;
      while(indSummed.IsValid()){ // loop over summed indicies

	gFactors = g(indSummed[0],indSummed[0 + indSummed.Size()/2]);
	// get the needed amount of metric tensor factors
	for(int i = 1; i < indSummed.Size()/2; i++){
	  gFactors *= g(indSummed[i],indSummed[i + indSummed.Size()/2]);
	}
	if(gFactors != 0.0){
	  nterm++;
	  // set up last ?? this and 1st ?? __tensor indicies
	  for(int i = ind1max; i < size1;i++) 
	    ind1.SetIndex(i,indSummed[i-ind1max]);
	  for(int i = size1 - ind1max; i < indSummed.Size(); i++)
	    ind2.SetIndex(i - (size1 - ind1max),indSummed[i]);

	  // multiply the metric tensor factor by this and __tensor elements
	  element = (this->Element(ind1))*(__tensor(ind2))*gFactors;

	  // add to each element this*__tensor*g*g...*g with correct # of g's
	  if(nterm == 1) ret(index) = element;	    
	  else ret(index) += element;	    
	}
	++indSummed;
      }
      // reset the summed indicies, step up index to next ret element
      indSummed.Reset();
      ++index;
    }
  }
  else{ // both are same rank tensors (R is rank 0)
    nterm = 0;

    // loop over summed indicies (only loop needed in this case)
    while(indSummed.IsValid()){
      gFactors = g(indSummed[0],indSummed[0 + indSummed.Size()/2]);
      
      // get the needed amount of metric tensor factors
      for(int i = 1; i < indSummed.Size()/2; i++)
	gFactors *= g(indSummed[i],indSummed[i + indSummed.Size()/2]);
      
      if(gFactors != 0.0){
	nterm++;
	for(int i = 0; i < indSummed.Size()/2 ;i++){
	  ind1.SetIndex(i,indSummed[i]);
	  ind2.SetIndex(i,indSummed[i + indSummed.Size()/2]);
	}
	element = (this->Element(ind1))*(__tensor(ind2))*gFactors;
	if(nterm == 1) ret() = element;
	else ret() += element;	  
      }
      indSummed++;
    } 
  }      
  return ret;
}
//_____________________________________________________________________________

template <typename _Tp> 
void Tensor<_Tp>::Boost(double __bx,double __by,double __bz){
   
  // check to see if bx,by,bz are all less than 1
  if(abs(__bx) >= 1 || abs(__by) >= 1 || abs(__bz) >= 1)
    cout << "Error! Attempt to boost using invalid boost vector." << endl;
  assert((abs(__bx) < 1)&&(abs(__by)<1)&&(abs(__bz)<1));

  Tensor<double> lt(2); // Lorentz transformation tensor
  double gamma = 1.0/sqrt(1.0 - __bx*__bx - __by*__by - __bz*__bz);
  double gamFact = (gamma*gamma)/(gamma + 1.0);

  // set up the Lorentz transformation tensor
  lt.Zero();

  lt(0,0) = gamma;
  lt(0,1) = gamma*__bx;
  lt(0,2) = gamma*__by;
  lt(0,3) = gamma*__bz;
  
  lt(1,1) = (__bx*__bx*gamFact)+1;
  lt(1,2) = __bx*__by*gamFact;
  lt(1,3) = __bx*__bz*gamFact;
  
  lt(2,2) = (__by*__by*gamFact)+1;
  lt(2,3) = __by*__bz*gamFact;
  
  lt(3,3) = (__bz*__bz*gamFact)+1;
  
  lt(1,0) = lt(0,1);
  lt(2,0) = lt(0,2);
  lt(2,1) = lt(1,2);
  lt(3,0) = lt(0,3);
  lt(3,1) = lt(1,3);
  lt(3,2) = lt(2,3);

  this->Transform(lt);
}
//_____________________________________________________________________________

template <typename _Tp> 
void Tensor<_Tp>::Rotate(double __alpha,double __beta,double __gamma){
    
  double ca = cos(__alpha);
  double sa = sin(__alpha);
  double cb = cos(__beta);
  double sb = sin(__beta);
  double cg = cos(__gamma);
  double sg = sin(__gamma);

  Tensor<double> lt(2); // Lorentz transformation tensor

  lt.Zero();

  lt(0,0) = 1.0;
  
  lt(1,1) = ca*cb*cg - sa*sg;
  lt(1,2) = sa*cb*cg + ca*sg;
  lt(1,3) = -sb*cg;

  lt(2,1) = -ca*cb*sg - sa*cg;
  lt(2,2) = -sa*cb*sg + ca*cg;
  lt(2,3) = sb*sg;

  lt(3,1) = ca*sb;
  lt(3,2) = sa*sb;
  lt(3,3) = cb;

  this->Transform(lt);
}
//_____________________________________________________________________________

template <typename _Tp> 
void Tensor<_Tp>::RotateX(double __alpha){
  double ca = cos(__alpha);
  double sa = sin(__alpha);
  Tensor<double> lt(2); // Lorentz transformation tensor
  lt.Zero();

  lt(0,0) = 1.0;
  lt(1,1) = 1.0;
  lt(2,2) = ca;
  lt(2,3) = -sa;
  lt(3,2) = sa;
  lt(3,3) = ca;

  this->Transform(lt);
}
//_____________________________________________________________________________

template <typename _Tp> 
void Tensor<_Tp>::RotateY(double __alpha){
  double ca = cos(__alpha);
  double sa = sin(__alpha);
  Tensor<double> lt(2); // Lorentz transformation tensor
  lt.Zero();

  lt(0,0) = 1.0;
  lt(1,1) = ca;
  lt(1,3) = sa;
  lt(2,2) = 1.0;
  lt(3,1) = -sa;
  lt(3,3) = ca;

  this->Transform(lt);
}
//_____________________________________________________________________________

template <typename _Tp> 
void Tensor<_Tp>::RotateZ(double __alpha){
  double ca = cos(__alpha);
  double sa = sin(__alpha);
  Tensor<double> lt(2); // Lorentz transformation tensor
  lt.Zero();

  lt(0,0) = 1.0;
  lt(1,1) = ca;
  lt(1,2) = -sa;
  lt(2,1) = sa;
  lt(2,2) = ca;
  lt(3,3) = 1.0;

  this->Transform(lt);
}
//_____________________________________________________________________________
template <typename _Tp>
void Tensor<_Tp>::Print(std::ostream& __os) const {
 
  if(_rank == 0) __os << "{Rank = 0 " << _data[0] << " }";
  else if(_rank == 1){
    __os << "{Rank = 1 ( " ;
    for(int mu = 0; mu < 3; mu++) __os << _data[mu] << ",";
    __os << _data[3] << ") } ";
  }
  else if(_rank == 2){
    int index;
    __os << "{Rank = 2 ";
    for(int mu = 0; mu < 4; mu++){
      __os << "(";
      for(int nu = 0; nu < 4; nu++){
	index = 4*nu + mu;
	__os << _data[index];
	if(nu < 3) __os << ",";
      }
      __os << ")";
      if(mu < 3) __os << ",";
    }
    __os << "}";
  }
  else{
    cout << "<Tensor::Print(ostream&)> Error! Can NOT print a Tensor with "
	 << " Rank > 2." << endl;
  }
}
//_____________________________________________________________________________

template <typename _Tp> template <typename T>
Tensor<typename MultType<_Tp,T>::Type>
Tensor<_Tp>::operator%(const Tensor<T> &__tensor) const {
  
  int rank = this->_rank + __tensor._rank;
  Tensor<typename MultType<_Tp,T>::Type> ret(rank);
  
  // if either tensor is rank 0, just return this*__tensor
  if((_rank == 0) || (__tensor._rank == 0)) return (*this)*__tensor;  
  else{ // we actually have to do some work
    TensorIndex index(rank);
    TensorIndex ind1(this->_rank);
    TensorIndex ind2(__tensor._rank);

    int size1 = ind1.Size();
    //    int size2 = ind2.Size(); 

    while(index.IsValid()){ // loop over ret's elements

      // set up this' indicies
      for(int i = 0; i < size1; i++) ind1.SetIndex(i,index[i]);     
      // set up _tensor's indicies
      for(int i = size1; i < index.Size(); i++) 
	ind2.SetIndex(i-size1,index[i]);
      
      // set element to product of this and __tensor elements
      ret(index) = (this->Element(ind1))*(__tensor(ind2));
      index++;
    }
  }
  return ret;
}
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
    cout << "Error! Attempt to reorder tensor indicies w/ incorrect number of"
	 <<" indicies." << endl;
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

template<typename _Tp> void Tensor<_Tp>::Transform(const Tensor<double> &__lt){

  if(__lt.Rank() != 2)
    cout << "Error! Lorentz transformation tensor NOT rank 2." << endl;
  assert(__lt.Rank() == 2);
  int rank = this->Rank();
  if(rank > 0) { // if rank 0 no transformation needed
    TensorIndex index(rank);
    TensorIndex indSummed(rank);
    int nterm;
    double lamFact;
    // make a copy
    Tensor<_Tp> copy(*this);
   
    while(index.IsValid()){  // loop over elements of this tensor
      nterm = 0;
      while(indSummed.IsValid()){
	// get the appropriate number of Lambda_mu_nu factors
	lamFact = __lt(index[0],indSummed[0]);
	for(int i = 1; i < rank; i++) lamFact *= __lt(index[i],indSummed[i]);
	if(lamFact != 0.0){
	  nterm++;
	  // add to each element this*Lambda*Lambda*...*Lambda 
	  if(nterm == 1) (*this)(index) = lamFact*(copy(indSummed));	
	  else (*this)(index) += lamFact*(copy(indSummed));	  
	}
	++indSummed;
      }
      // reset summed indicies, step up index to next element
      indSummed.Reset();
      ++index;
    }
  }
}
//_____________________________________________________________________________

#endif /* _Tensor_TCC */
