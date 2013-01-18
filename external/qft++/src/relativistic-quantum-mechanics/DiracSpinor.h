// DiracSpinor class definition file -*- C++ -*-
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
#ifndef _DiracSpinor_H
#define _DiracSpinor_H
//_____________________________________________________________________________
/** @file DiracSpinor.h
 *  @brief DiracSpinor class definition file
 */
//_____________________________________________________________________________
#include <iostream>
#include <complex>
#include "PolVector.h"
#include "Utils.h"
#include "Spin.h"
#include "../../include/matrix.h"
#include "../../include/tensor.h"
#include "DiracMatricies.h"

using namespace std;
//_____________________________________________________________________________
/** @class DiracSpinor
 *  @author Mike Williams
 *
 *  @brief \f$u_{\mu_1\mu_2\ldots\mu_{S-1/2}}(p,m)\f$: Spin-S spinors for any half-integer spin.
 *
 * This class handles particles with half-integer spin. The spin of the 
 * particle is set when the constructor is called. The class stores 
 * \f$ (2S+1)\f$ spinors in a std::vector of 
 * Matrix<Tensor<complex<double> > >'s along with the particle's 4-momentum 
 * (Vector4<double>), the spin-S projection operator as a 
 * Matrix<Tensor<complex<double> > > along with the mass and spin (Spin). 
 * The spin-S spinors are \f$4\times 1\f$ column matricies
 * of rank \f$S-1/2\f$ tensors that satisfy the Rarita-Schwinger conditions,
 * \f[(\gamma^{\mu}p_{\mu}-w) u_{\mu_1\mu_2\ldots\ldots\mu_n}(p,m) = 0 \f]
 * \f[u_{\mu_1\mu_2\ldots\mu_i\ldots\mu_j\ldots\mu_n}(p,m) 
 * = u_{\mu_1\mu_2\ldots\mu_j\ldots\mu_i\ldots\mu_n}(p,m) \f]
 * \f[p^{\mu_i}u_{\mu_1\mu_2\ldots\mu_i\ldots\mu_n}(p,m) = 0 \f]
 * \f[\gamma^{\mu_i}u_{\mu_1\mu_2\ldots\mu_i\ldots\mu_n}(p,m) = 0 \f]
 * \f[g^{\mu_i\mu_j}u_{\mu_1\mu_2\ldots\mu_i\ldots\mu_j\ldots\mu_n}(p,m) =0 \f]
 * for any \f$\mu_i,\mu_j\f$. This reduces the number of independent elements
 * from \f$4^{S+1/2}\f$ to \f$(2S+1)\f$.
 * 
 * For a spin-1/2 particle, the spinors are,
 * \f[ u(p,m) = \sqrt{E+w}\left( \begin{array}{c} \chi(m) \\ 
 * \frac{\vec{\sigma}\cdot\vec{p}}{E + w} \chi(m) \end{array}  \right) \f]
 * where \f$p,m \f$ are the 4-momentum,spin projection of the particle and
 * \f[ \chi(+\frac{1}{2}) = \left(\begin{array}{c}1\\0\end{array}\right) 
 * \qquad \chi(-\frac{1}{2}) = \left(\begin{array}{c}0\\1\end{array}\right) \f]
 *
 * Higher spin spinors are built from tensor products of the spin-1/2 spinors
 * and the spin-(n = S-1/2) polarization vectors according to,
 * \f[u_{\mu_1\mu_2\ldots\ldots\mu_n}(p,m) = \sum\limits_{m_nm_{\frac{1}{2}}}
 * (nm_n\frac{1}{2}m_{\frac{1}{2}}|Sm)\epsilon_{\mu_1\mu_2\ldots\ldots\mu_n}
 * (p,m_n)u(p,m_{\frac{1}{2}})\f]
 * where \f$(nm_n\frac{1}{2}m_{\frac{1}{2}}|Sm)\f$ are the Clebsch-Gordon
 * coefficents coupling spin-(n = S-1/2) and spin-1/2 to spin-S
 * (see Clebsch()). These are set in the function SetP4(). Access to the 
 * spinors is provided through the () operator. Since they are stored as 
 * Matrix<Tensor<complex<double> > >'s, they can be used with other Matrix,
 * Tensor objects and inherited objects of either very efficently.
 *
 * <b> Example Usage </b>
 *
 * The amplitude for \f$ \Delta \rightarrow p \pi \f$ can be written as,
 * \f[ \bar{u}(p_p,m_p)p_p^{\mu}u_{\mu}(p_{\Delta},m_{\Delta}) \f]
 * where \f$p_x,m_x\f$ are the 4-momentum,spin projection of particle x. We 
 * can write this in the code as,
 * \include DiracSpinor_Delta.ex
 *
 *
 * The spin-S projection operator for half integer spins is defined as,
 * \f[ P^{(S)}_{\mu_1\mu_2\ldots\mu_{S-1/2}\nu_1\nu_2\ldots\nu_{S-1/2}}(P) 
 * = \sum\limits_{M} u_{\mu_1\mu_2\ldots\mu_{S-1/2}}(P,M) 
 * \bar{u}_{\nu_1\nu_2\ldots\nu_{S-1/2}}(P,M) \f] 
 * where \f$P,M\f$ are the 4-momentum, spin projection and 
 * \f$-S\leq M \leq S\f$
 *
 * For a spin-1/2 particle, this takes on the familar form of,
 * \f[ P^{(\frac{1}{2})}(P) = \frac{1}{2w}(\gamma^{\mu}P_{\mu} + w) \f]
 * where \f$ w \f$ is the mass. For details on how this is calculated for 
 * higher spins, see SetProjector(). This operator is stored internally in the 
 * class for performance purposes. For higher spins, the projection operator 
 * becomes a rather large object that becomes costly to pass by value.
 * DiracSpinor sets _projector when SetP4() is called and returns a constant
 * reference to it via the Projector() member function. DiracSpinor also 
 * provides member functions to return the Feynman and Briet-Wigner propagators
 * (see Propagator() and PropagatorBW()). Since these aren't stored by the 
 * class, they must be returned by value. Thus, for higher spins this will be 
 * slow. If the projector is going to be fully contracted with other Tensor's
 * to form a scalar (ie. calculating an amplitude), the best approach is to 
 * contract Projector() with the other Tensor's and then divide the result by 
 * the appropriate propagator.
 * 
 * <b> Example Usage </b>
 *
 * The amplitude for \f$ p\pi \rightarrow \frac{3}{2}^- \rightarrow p\pi \f$
 * can be written as,
 * \f[ \bar{u}(p_f,m_f)\gamma^5 p_f^{\mu} P^{(3/2)}_{\mu\nu}(P) 
 * \gamma^5 p_i u(p_i,m_i) \f]
 * where \f$p_f(p_i)\f$ are the final(initial) momenta, \f$m_f(m_i)\f$ are the 
 * final(initial) spin projections and \f$P = p_f + p^f_{\pi}\f$ is the 
 * momenta of the \f$3/2^-\f$. We can write this in the code as,
 *
 * \include DiracSpinor_Amp.ex
 *
 */
//_____________________________________________________________________________

class DiracSpinor {

private:
  // attributes:
  
  /// Spinors, _spinors[n] = \f$ u_{\mu_1\mu_2\ldots\mu_{S-1/2}}(P,n-S) \f$
  vector<Matrix<Tensor<complex<double> > > > _spinors; 
  Vector4<double> _p4; ///< The particle's 4-momentum
  double _mass; ///< The particle's mass
  Spin _spin; ///< Spin of the particle
  Matrix<Tensor<complex<double> > > _projector; ///< Spin S Projection Operator

  // private functions:

  /** Initialize for spin */
  void _Init(const Spin &__spin);
   
  /** Copy @a dspin to @a this */
  void _Copy(const DiracSpinor &__dspin);
   
  /** Set the spin-S projection operator.
   *
   * For a spin-1/2 particle, the projector is set to be,
   * \f[ P^{(\frac{1}{2})}(P) = \frac{1}{2w}(\gamma^{\mu}P_{\mu} + w) \f]
   * where \f$P\f$ is the 4-momentum, \f$ w \f$ is the mass,
   * and \f$M = \pm \frac{1}{2} \f$ is the spin projection.
   *
   * For higher spins, the projector is calculated using,
   * \f[ p^{(S)}_{\mu_1\mu_2\ldots\mu_n\nu_1\nu_2\ldots\nu_n}(P)
   * = \frac{1}{2w}\sum\limits_M u_{\mu_1\mu_2\ldots\mu_n}(P,M)
   * \bar{u}_{\nu_1\nu_2\ldots\nu_n}(P,M)\f]
   *
   * where \f$ n = S - 1/2 \f$.
   *
   */
  void _SetProjector();

public:

  // create/copy/destroy:

  /// Default Constructor (spin 1/2)
  DiracSpinor(){
    this->_Init(Spin(1,2));
    _mass = 0.;
  }

  /// Constructor
  DiracSpinor(const Spin &__spin){
    this->_Init(__spin);
    _mass = 0.;
  }

  /// Copy Constructor
  DiracSpinor(const DiracSpinor &__dspin){
    this->_Copy(__dspin);
  }

  /** Destructor */
  virtual ~DiracSpinor(){}

  // Setters:

  /// Set the spin of the spinor
  void SetSpin(const Spin &__spin){
    this->Clear();
    this->_Init(__spin);
    _mass = 0.;
  }

  /** Set the DiracSpinor object for a given 4-momentum.
   *
   * @param p4 The particle's 4-momentum
   * @param mass Pole mass (used for normalization)
   *
   * For spin-1/2 particles, the spinors are defined as,
   * \f[ u(P,M) = \sqrt{E+w}\left( \begin{array}{c} \chi(M) \\ 
   * \frac{\vec{\sigma}\cdot\vec{P}}{E + w} \chi(M) \end{array}  \right) \f]
   * where \f$P,M \f$ are the 4-momentum,spin projection of the particle and
   * \f[ \chi(+\frac{1}{2}) = \left(\begin{array}{c}1\\0\end{array}\right) 
   * \qquad \chi(-\frac{1}{2}) =\left(\begin{array}{c}0\\1\end{array}\right)\f]
   *
   * Higher spin spinors are built from tensor products of the spin-1/2 spinors
   * and the spin-(n = S-1/2) polarization vectors according to,
   * \f[u_{\mu_1\mu_2\ldots\ldots\mu_n}(P,M) = \sum\limits_{m_nm_{\frac{1}{2}}}
   * (nm_n\frac{1}{2}m_{\frac{1}{2}}|SM)\epsilon_{\mu_1\mu_2\ldots\ldots\mu_n}
   * (P,m_n)u(P,m_{\frac{1}{2}})\f]
   * where \f$(nm_n\frac{1}{2}m_{\frac{1}{2}}|SM)\f$ are the Clebsch-Gordon
   * coefficents coupling spin-(n = S-1/2) and spin-1/2 to spin-S
   * (see Clebsch()).
   */
  void SetP4(const Vector4<double> &__p4,double __mass);

  // Getters:
  
  /** Returns \f$u_{\mu_1\mu_2\ldots\mu_{S-1/2}}(p,m_z)\f$ 
   *
   * @param mz Spin projection along the z-axis (\f$ -S \leq m_z \leq S\f$)
   */
  const Matrix<Tensor<complex<double> > >& operator()(const Spin &__mz) const {
    int i = (int)(_spin + __mz);
    if(i < 0 || i > 2*_spin) {
      cout << "Error! Attempt to access non-existent spinor (mz out of range)."
	   << endl;
    }
    assert(i >= 0 && i <= (int)(2*_spin));
    return _spinors[i];
  }

  /** Returns the spin projection operator for spin S
   *
   * The spin projection operator is defined as,
   * \f[ P^{(S)}_{\mu_1\mu_2\ldots\mu_{S-1/2}\nu_1\nu_2\ldots\nu_{S-1/2}}(P) 
   * = \frac{1}{2w} \sum\limits_{M} u_{\mu_1\mu_2\ldots\mu_{S-1/2}}(P,M) 
   * \bar{u}_{\nu_1\nu_2\ldots\nu_{S-1/2}}(P,M) \f] 
   * where \f$P,M\f$ are the 4-momentum, spin projection and 
   * \f$-S\leq M \leq S\f$.
   *
   * <b> Example: </b>
   *
   * \f[ P^{(\frac{3}{2})}_{\mu\nu}(P) = \sum\limits_{M}u_{\mu}(P,M)
   * \bar{u}_{\nu}(P,M) = \frac{1}{2w}(P_{\rho}\gamma^{\rho}+w)\left[-
   * \tilde{g}_{\mu\nu} + \frac{1}{3}\tilde{g}_{\mu\alpha}\gamma^{\alpha}
   * \tilde{g}_{\nu\beta}\gamma^{\beta}\right] \f]
   * where \f$\tilde{g}_{\mu\nu} = g_{\mu\nu} - \frac{P_{\mu}P_{\nu}}{w^2} \f$
   * and \f$w\f$ is the mass.
   *
   * See SetProjector() for details on how this is calculated.
   */
  const Matrix<Tensor<complex<double> > >& Projector() const {
    return _projector;
  }

  /** Returns the Feynman propagator for spin S
   *
   * \returns \f[\frac{2w P^{(S)}_{\mu_1\mu_2\ldots\mu_{S-1/2}
   *                           \nu_1\nu_2\ldots\nu_{S-1/2}}(P)}{P^2 - w^2}\f]
   *
   * with \f$w\rightarrow\f$_mass, \f$P\rightarrow\f$_p4 and 
   * \f$P^{(S)}\rightarrow\f$Projector().
   *
   * Note: This is slow, for speed use Projector() and do the division after
   *       any matrix/tensor operations.
   */
  Matrix<Tensor<complex<double> > > Propagator() const {
    Matrix<Tensor<complex<double> > > prop(4,4);  
    prop = (this->Projector()*2.*_mass)/(_p4.M2() - _mass*_mass);
    return prop;
  }

  /** Returns the Breit-Wigner propagator for spin S
   *
   * @param width Breit-Wigner width.
   * @returns \f[ \frac{2w P^{(S)}_{\mu_1\mu_2\ldots\mu_{S-1/2}
   *        \nu_1\nu_2\ldots\nu_{S-1/2}}(P) w\Gamma}{P^2 - w^2 + iw\Gamma}\f]
   * with \f$\Gamma\rightarrow\f$width, \f$w\rightarrow\f$_mass, 
   * \f$P\rightarrow\f$_p4 and \f$P^{(S)}\rightarrow\f$Projector().
   *
   * Note: This is slow, for speed use Projector() and do the division after
   *       any matrix/tensor operations.
   */
  Matrix<Tensor<complex<double> > > PropagatorBW(double __width) const {
    Matrix<Tensor<complex<double> > > prop(4,4);  
    prop = (this->Projector()*2.*_mass)*BreitWigner(_p4,_mass,__width);
    return prop;
  }
  
  /// Return the particle's 4-momentum
  inline const Vector4<double>& GetP4() const {
    return _p4;
  }

  /// Return the particle's mass
  inline double GetMass() const {
    return _mass;
  }

  // Functions:

  /** Boost the spinor.
   *
   * The boost vector is given by \f$\vec{\beta} = (bx,by,bz)\f$
   *
   * For spin-1/2 particles,
   * boosting is done by constructing the boost operator for spinors,
   * \f[ D(\vec{\beta}) = \left(\begin{array}{cc} cosh(\alpha/2) & 
   * \vec{\sigma}\cdot\hat{P} sinh(\alpha/2) \\ \vec{\sigma}\cdot\hat{P} 
   * sinh(\alpha/2) & cosh(\alpha/2) \end{array} \right) \f]
   * where \f$ tanh(\alpha) = \beta \f$. Then the boosted spinor is given by,
   * \f[ u(P,M) = D(\vec{\beta}) u(P',M) \f]
   * where \f$P(P')\f$ is the momentum in the boosted(original) frame.
   *
   * For higher spins, a spin-1/2 spinor, \f$u(P',M)\f$, and a 
   * spin-(n = S - 1/2) polarization vector, 
   * \f$\epsilon_{\mu_1\mu_2\ldots\mu_n}(P',M) \f$, are constructed in the 
   * current frame. These are then boosted using the   
   * boost vector \f$\vec{\beta}\f$ and the spin-S spinor is reconstructed out
   * of the boosted spinor and polarization vector according to the equation
   * given in SetP4(). 
   */
  void Boost(double __bx,double __by,double __bz);

  /// Boost the spinor to p4's rest frame (see Boost(double,double,double)).
  void Boost(const Vector4<double> &__p4) {
    double p = __p4.P();
    this->Boost(-__p4.X()/p,-__p4.Y()/p,-__p4.Z()/p);
  }

  /// Clear the DiracSpinor object
  void Clear() {
    if(!_spinors.empty()) _spinors.clear();
    _spin = 0;
  }

  /// Zero the spinors only.
  void Zero() {
    int size = (int)_spinors.size();
    for(int i = 0; i < size; i++) _spinors[i].Zero();
  }

  /// Returns \f$u_{\mu_1\mu_2\ldots\mu_{S-1/2}}(P,\lambda) \f$ 
  Matrix<Tensor<complex<double> > > Helicity(const Spin &__lam) const {
    Matrix<Tensor<complex<double> > > hel(4,1);
    SetRank(hel,(int)(_spin - 1/2.));
    hel.Zero();
    for(Spin m = -_spin; m <= _spin; m++){
      hel += Wigner_D(_p4.Phi(),_p4.Theta(),0.,_spin,m,__lam)
	*_spinors[(int)(_spin + m)];
    }
    return hel;
  }

  // static functions:
  
  /** Returns the spin-<em>j</em> projector as a 2 x @a rank tensor.
   *
   * @param j Spin of the projector
   * @param rank Tensor rank of projector will be 2 x rank
   * @param p4 4-momentum
   * @param mass Mass 
   * @param projector Will be set to the projector
   *
   * Builds the spin-<em>j</em> projection operator in @a rank space.
   */
  static void Projector(const Spin &__j,int __rank,const Vector4<double> &__p4,
			double __mass,
			Matrix<Tensor<complex<double> > > &__projector);

};
//_____________________________________________________________________________

#endif /* _DiracSpinor_H  */
