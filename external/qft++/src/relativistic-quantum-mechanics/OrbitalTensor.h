// OrbitalTensor class definition file. -*- C++ -*-
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
#ifndef _OrbitalTensor_H
#define _OrbtialTensor_H
//_____________________________________________________________________________
/** @file OrbitalTensor.h
 *  @brief OrbitalTensor class definition file
 */
//_____________________________________________________________________________

#include "../../include/tensor.h"
#include <iostream>

using namespace std;
//_____________________________________________________________________________
/** @class OrbitalTensor
 *  @author Mike Williams
 *
 *  @brief \f$L^{(\ell)}_{\mu_1\mu_2\ldots\mu_{\ell}}\f$ : Orbital angular momentum tensors.
 *
 * OrbitalTensor inherits from Tensor<double>. The only extension to Tensor 
 * implemented in OrbitalTensor can be found in the member function SetP4. 
 * This function takes 2 4-momenta (in the form of Vector4<double>) and 
 * constructs out of them a state of pure orbital angular momentum (with 
 * \f$ \ell\f$  = the rank of the OrbitalTensor object). 
 *
 * Orbital angular momentum tensors coupling particles originating from the 
 * decay of a composite particle with 4-momentum \f$ P \f$ satify the 
 * Rarita-Schwinger conditions:
 *
 * \f[ P^{\mu_i}L^{(\ell)}_{\mu_1\mu_2\ldots\mu_i\ldots\mu_{\ell}}(p_{ab})=0\f]
 * \f[ L^{(\ell)}_{\mu_1\mu_2\ldots\mu_i\ldots\mu_j\ldots\mu_{\ell}}(p_{ab}) 
 *  = L^{(\ell)}_{\mu_1\mu_2\ldots\mu_j\ldots\mu_i\ldots\mu_{\ell}}(p_{ab}) \f]
 * \f[ g^{\mu_i\mu_j}
 * L^{(\ell)}_{\mu_1\mu_2\ldots\mu_i\ldots\mu_j\ldots\mu_{\ell}}(p_{ab}) = 0\f]
 *
 * for any \f$\mu_i,\mu_j\f$, which insures that they have \f$(2\ell+1)\f$ 
 * independent elements. 
 *
 * If the 4-momenta of the decay products are, \f$ p_a\f$
 * and \f$ p_b \f$, then the orbital tensors for \f$\ell = 0,1,2,3 \f$ are,
 * \f[ L^{(0)}(p_{ab}) =  1 \f]
 * \f[ L^{(1)}_{\mu}(p_{ab}) = \tilde{p}^{ab}_{\mu} \f] 
 * \f[ L^{(2)}_{\mu_1\mu_2}(p_{ab}) = \frac{3}{2} \left(\tilde{p}^{ab}_{\mu_1}
 *      \tilde{p}^{ab}_{\mu_2} - \frac{1}{3}\tilde{p}_{ab}^2
 *      \tilde{g}_{\mu_1\mu_2} \right)\f]
 * \f[L^{(3)}_{\mu_1\mu_2\mu_3}(p_{ab})=\frac{5}{2}\left(\tilde{p}^{ab}_{\mu_1}
 *      \tilde{p}^{ab}_{\mu_2}\tilde{p}^{ab}_{\mu_3} - \frac{1}{5}
 *      \tilde{p}_{ab}^2(\tilde{g}_{\mu_1\mu_2}\tilde{p}^{ab}_{\mu_3} 
 *      + \tilde{g}_{\mu_1\mu_3}\tilde{p}^{ab}_{\mu_2} 
 *      + \tilde{g}_{\mu_2\mu_3}\tilde{p}^{ab}_{\mu_1}) \right)\f]
 * where \f$p_{ab} = \frac{1}{2}(p_a - p_b)\f$, 
 * \f$\tilde{g}_{\mu_1\mu_2} = g_{\mu_1\mu_2}-\frac{P_{\mu_1}P_{\mu_2}}{P^2}\f$
 *  and \f$\tilde{p}^{ab}_{\mu} = \tilde{g}_{\mu\nu}p_{ab}^{\nu} \f$.
 *
 * Note: Normalization used here follows Anisovich, et. al, which is given by
 * \f[ L^{(\ell)}_{\mu_1\mu_2\ldots\mu_{\ell}}(p_{ab}) = \alpha(\ell)(-)^{\ell}
 *     P^{(\ell)}_{\mu_1\mu_2\ldots\mu_{\ell}\nu_1\nu_2\ldots\nu_{\ell}}(P)
 *     p_{ab}^{\nu_1} p_{ab}^{\nu_2} \ldots p_{ab}^{\nu_{\ell}} \f]
 *
 * where \f$ \alpha(\ell) = \frac{(2\ell - 1)!!}{\ell !} \f$
 */
//_____________________________________________________________________________

class OrbitalTensor : public Tensor<double> {

public:

  // create/copy/destroy:

  /// Default Constructor (rank 0)
  OrbitalTensor() : Tensor<double>::Tensor() {}

  /// Constructor 
  /** Construct and OrbitalTensor for \f$\ell\f$ = @a rank.*/
  OrbitalTensor(int __rank) : Tensor<double>::Tensor(__rank) {}

  /// Copy Constructor
  OrbitalTensor(const OrbitalTensor &__orb):Tensor<double>::Tensor(__orb){}

  /// Destructor
  virtual ~OrbitalTensor(){}

  // operators:

  /// Assignment operator
  OrbitalTensor& operator=(const Tensor<double> &__tensor){
    this->Tensor<double>::operator=(__tensor);
    return *this;
  }

  /// Assignment operator
  OrbitalTensor& operator=(double __x){
    this->Tensor<double>::operator=(__x);
    return *this;
  }

  // functions:

  /** Construct the orbital angular momentum tensor for \f$X\rightarrow ab\f$
   * 
   * \param p4a 4-momentum of decay particle a
   * \param p4b 4-momentum of decay particle b
   * 
   * See the class description above for details.
   */  
  void SetP4(const Vector4<double> &__p4a,const Vector4<double> &__p4b){
    int rank = this->Rank();
    if(rank == 0) *this = 1.;
    else{
      Tensor<double> kperp(1);
      Tensor<double> gbar(2);
      MetricTensor g;

      gbar = g - ((__p4a + __p4b)%(__p4a + __p4b))/((__p4a + __p4b).M2());
      kperp = 0.5*(__p4a - __p4b)*gbar;

      double kperp_2 = kperp*kperp;

      if(rank == 1) *this = kperp;
      else if(rank == 2) *this = 1.5*(kperp%kperp - (kperp_2/3.)*gbar);
      else if(rank == 3) { 
	// the factor of 3 is because Symmetric normalizes by number of terms
	*this = (5./2.)*(kperp%kperp%kperp 
			 - (kperp_2/5.)*(3.*(gbar%kperp).Symmetric()));
      }
      else if(rank == 4){
	// factors of 6 and 3 are to cancel Symmetric's normalization
	*this = (35./8.)*(kperp%kperp%kperp%kperp 
			  - (kperp_2/7.)*(6.*(gbar%kperp%kperp).Symmetric())
			  +(kperp_2*kperp_2/35.)*(3.*(gbar%gbar).Symmetric()));
      }
      else{ // rank > 4 
	int rank = 5;
	Tensor<double> z(6);
	OrbitalTensor x(4);
	x.SetP4(__p4a,__p4b);

	while(rank <= this->Rank()){
	  z.SetRank(rank + 1);
	  z = (x%gbar).Permute(0,rank-1);
	  for(int i = 1; i < rank; i++){
	    z += (x%gbar).Permute(i,rank-1);
	    for(int j = 0; j < rank; j++){
	      if(j < i){
		z += (-2./(2.*rank-1.))*(((gbar%x).Permute(0,i)).Permute(1,j));
	      }
	    }
	  }
	  z *= (2.*rank-1.)/(rank*rank);
	  x.SetRank(rank);
	  x = z*kperp;
	  rank++;
	}
	(*this) = x;
      }
    }
  }

};
//_____________________________________________________________________________

#endif /* _OrbtialTensor_H */
