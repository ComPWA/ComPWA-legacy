// DiracAntiSpinor class definition file -*- C++ -*-
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
#ifndef _DiracAntiSpinor_H
#define _DiracAntiSpinor_H
//_____________________________________________________________________________
/** @file DiracAntiSpinor.h
 *  @brief DiracAntiSpinor class definition file
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
/** @class DiracAntiSpinor
 *  @author Mike Williams
 *
 *  @brief \f$ v(p,m) \f$: spin-1/2 anti-particle spinors.
 *
 * This class handles spin-1/2 anti-particles. The spinors are given by,
 * \f[ v(p,m) = \sqrt{E+w} \left( \begin{array}{c} \frac{\vec{p} 
 * \cdot \vec{\sigma}}{E+w} \chi(m) \\ \chi(m) \end{array}\right) \f]
 * where \f$w\f$ is the mass, \f$E\f$ is the energy and,
 * \f[ \chi(+\frac{1}{2}) = \left(\begin{array}{c} 0\\1 \end{array}\right)
 * \qquad \chi(-\frac{1}{2}) = \left(\begin{array}{c} 1\\0 \end{array}\right)
 * \f] 
 *
 * They can be accessed through the use of the () operator. DiracAntiSpinor 
 * also provides member functions to return the anti-particle projection 
 * operator, Projector(), and the Feynman and Breit-Wigner propagators 
 * (Propagator() and PropagatorBW()).
 *
 * <b> Example Usage: </b>
 *
 * The amplitude for \f$ e^+ e^- \rightarrow \gamma \rightarrow \mu^+ \mu^- \f$
 *  is proportional to
 * \f[ \bar{u}(p_{\mu^+},m_{\mu^+}) \gamma_{\rho} v(p_{\mu^-},m_{\mu^-}) 
 * \frac{1}{P^2} \bar{v}(p_{e^+},m_{e^+}) \gamma^{\rho}u(p_{e^-},m_{e^-}) \f]
 * where \f$p\f$'s are the momenta, \f$m\f$'s are the spin projections and 
 * \f$P\f$ is the momenta of the virtual photon. We can write this in the code
 * as follows,
 *
 * \include DiracAntiSpinor.ex
 *
 */
//_____________________________________________________________________________

class DiracAntiSpinor {

private:
  // data members:

  Matrix<complex<double> > _spinors[2]; ///< spin up/down anti-particle spinors
  Vector4<double> _p4; ///< The particle's 4-momentum
  double _mass; ///< The particle's mass
  Matrix<complex<double> > _projector; ///< Spin-1/2 projector

  // private functions:

  /// Initialize matricies
  void _Init(){
    for(int i = 0; i < 2; i++) _spinors[i].Resize(4,1);
    _projector.Resize(4,4);
  }

  /// Copy @a aspin to @a this
  void _Copy(const DiracAntiSpinor &__aspin){
    this->_Init();
    _p4 = __aspin._p4;
    _mass = __aspin._mass;
    _projector = __aspin._projector;
    for(int i = 0; i < 2; i++) _spinors[i] = __aspin._spinors[i];
  }
  
  /** Set _projector for the current _p4.
   *
   * The spin-1/2 anti-particle projection operator is defined as,
   * \f[ P^{(\frac{1}{2})}(P) = \frac{1}{2w}
   * \sum\limits_{M}v(P,M)\bar{v}(P,M) 
   * = \frac{1}{2w}(\gamma^{\mu}P_{\mu} - w) \f]
   * where \f$P\f$ is the 4-momentum, \f$ w \f$ is the mass,
   * and \f$M = \pm \frac{1}{2} \f$ is the 
   * spin projection.
   */
  void _SetProjector(){
    IdentityMatrix<double> I(4);
    DiracGamma gamma;
    _projector = (_p4*gamma - I*_mass)/(2*_mass);
  }  


public:

  // create/copy/destroy:

  /// Default Constructor (spin 1/2)
  DiracAntiSpinor(){
    _mass = 0.;
    this->_Init();
  }

  /// Copy Constructor
  DiracAntiSpinor(const DiracAntiSpinor &__aspin){
    this->_Copy(__aspin);
  }
  
  /** Destructor */
  virtual ~DiracAntiSpinor() {}

  // Setters:

  /** Set the DiracAntiSpinor object for a given 4-momentum.
   *
   * @param p4 4-momentum
   * @param mass mass (used to normalize the spinor)
   *
   * The anti-particle spinor can be written as,
   * \f[ v(p,m) = \sqrt{E+w} \left( \begin{array}{c} \frac{\vec{p} 
   * \cdot \vec{\sigma}}{E+w} \chi(m) \\ \chi(m) \end{array}\right) \f]
   * where \f$w\f$ is the mass, \f$E\f$ is the energy and,
   * \f[ \chi(+\frac{1}{2}) = \left(\begin{array}{c} 0\\1 \end{array}\right)
   * \qquad \chi(-\frac{1}{2}) = \left(\begin{array}{c} 1\\0 \end{array}\right)
   *  \f]
   *
   * This function also sets the projection operator (see SetProjector()).
   */
  void SetP4(const Vector4<double> &__p4,double __mass);

  // Getters:

  /// Returns \f$ v(m_z) \f$.
  const Matrix<complex<double> >& operator()(const Spin &__mz) const{
    int index = (int)(__mz + 1/2.);
    if(index < 0 || index > 1) 
      cout << "Error! Attempt to access non-existent spinor." << endl;
    assert(index >= 0 && index <= 1);
    return _spinors[index];
  }

  /** Returns the spin-1/2 anti-particle projection operator.
   *
   * The spin-1/2 anti-particle projection operator is defined as,
   * \f[ P^{(\frac{1}{2})}(P) = \sum\limits_{M}v(P,M)\bar{v}(P,M) 
   * = \frac{1}{2w}(\gamma^{\mu}P_{\mu} - w) \f]
   * where \f$P\f$ is the 4-momentum, \f$ w \f$ is the mass,
   * and \f$M = \pm \frac{1}{2} \f$ is the 
   * spin projection.
   */
  const Matrix<complex<double> >& Projector() const {
    return _projector;
  }

  /** Returns the Feynman propagator for a spin-1/2 anti-particle.
   *
   * \returns \f[ \frac{2w P^{(\frac{1}{2})}(P)}{P^2 - w^2} \f]
   * where \f$w\rightarrow\f$ _mass, \f$P\rightarrow\f$ _p4 and 
   * \f$P^{(\frac{1}{2})}(P)\f$ is returned by Projector().
   */
  Matrix<complex<double> > Propagator() const {
    Matrix<complex<double> > prop(4,4);    
    prop = (this->_projector*2.*_mass)/(_p4.M2() - _mass*_mass);
    return prop;
  }

  /** Returns the Breit-Wigner propagator for a spin-1/2 anti-particle.
   *
   * @param width Width of the Breit-Wigner.
   *
   * \returns \f[ \frac{2w P^{(\frac{1}{2})}(P) w\Gamma}
   * {P^2 - w^2 + iw\Gamma} \f]
   * where \f$\Gamma\rightarrow\f$ width, \f$w\rightarrow\f$ _mass, 
   * \f$P\rightarrow\f$ _p4 and 
   * \f$P^{(\frac{1}{2})}(P)\f$ is returned by Projector().
   */
  Matrix<complex<double> > PropagatorBW(double __width) const{
    Matrix<complex<double> > prop(4,4);  
    prop = (this->_projector*2.*_mass)*BreitWigner(_p4,_mass,__width);
    return prop;
  }

  // Functions:
  
  /** Boost the spinor.
   * 
   * The boost vector is defined as \f$\vec{\beta} = (bx,by,bz)\f$
   *
   * Boosting is done by constructing the boost operator for spinors,
   * \f[ D(\vec{\beta}) = \left(\begin{array}{cc} cosh(\alpha/2) & 
   * \vec{\sigma}\cdot\hat{P} sinh(\alpha/2) \\ \vec{\sigma}\cdot\hat{P} 
   * sinh(\alpha/2) & cosh(\alpha/2) \end{array} \right) \f]
   * where \f$ tanh(\alpha) = \beta \f$. Then the boosted spinor is given by,
   * \f[ u(P,M) = D(\vec{\beta}) u(P',M) \f]
   * where \f$P(P')\f$ is the momentum in the boosted(original) frame.
  */
  void Boost(double __bx,double __by,double __bz);

  /// Boost the spinor to p4's rest frame (see Boost(double,double,double)).
  void Boost(const Vector4<double> &__p4){
    double p = __p4.P();
    this->Boost(-__p4.X()/p,-__p4.Y()/p,-__p4.Z()/p);
  }

};
//_____________________________________________________________________________

#endif /* _DiracAntiSpinor_H */
