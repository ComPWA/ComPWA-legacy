// PolVector class definition file -*- C++ -*-
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
#ifndef _PolVector_H
#define _PolVector_H
//_____________________________________________________________________________
/** @file PolVector.h
 *  @brief PolVector class definition file
 */
//_____________________________________________________________________________
// Headers:
#include <iostream>
#include <complex>
#include <vector>
#include "../../include/tensor.h"
#include "Utils.h"
#include "Spin.h"

using namespace std;
//_____________________________________________________________________________
/** @class PolVector
 *  @author Mike Williams
 *
 * @brief \f$\epsilon_{\mu_1\mu_2\ldots\mu_S}\f$ : Spin-S polarization vectors for any integer spin.
 *
 * This class handles particles with integer spin. The spin of the particle is
 * set when the constructor is called. The class stores the 2S + 1 polarization
 * vectors in a std::vector of Tensor<complex<double> >'s along with
 * the particle's 4-momentum (in a Vector4<double>), the spin-S projector as a
 * Tensor<complex<double> > and the mass and spin. The spin-S
 * polarization vectors are rank S tensors that satisfy the Rarita-Schwinger
 * conditions,
 * \f[ P^{\mu_i}\epsilon_{\mu_1\mu_2\ldots\mu_i\ldots\mu_S}(P,M) = 0 \f]
 * \f[\epsilon_{\mu_1\mu_2\ldots\mu_i\ldots\mu_j\ldots\mu_S}(P,M) 
 *     = \epsilon_{\mu_1\mu_2\ldots\mu_j\ldots\mu_i\ldots\mu_S}(P,M) \f]
 * \f[ g^{\mu_i\mu_j}\epsilon_{\mu_1\mu_2\ldots\mu_i\ldots\mu_j\ldots\mu_S}
 *     (P,M) = 0 \f]
 * for any \f$\mu_i,\mu_j\f$. This reduces the number of independent elements 
 * from \f$4^S\f$ to \f$(2S+1)\f$.
 *
 * For a particle with spin S = 1 (_spin = 1), there are 2 cases. If 
 * mass > 0, then the rest frame polarization vectors are constructed as,
 * \f[ \epsilon^{\mu}(0,m) = (0,\vec{\epsilon}(m)) \f] where
 * \f[ \vec{\epsilon}(\pm 1) = \mp \frac{1}{\sqrt{2}}(1,\pm i,0),\hspace{0.2in}
 *       \vec{\epsilon}(0) = (0,0,1) \f]
 * These are then boosted to the frame where the particle has 4-momentum
 * p4. For the mass = 0 case (the photon), the \f$ M_z = 0 \f$ polarization
 * vector is set to be 0, while the \f$ M_z = \pm 1 \f$ polarization 
 * vectors are set to be the helicity states for the photon. For example,
 * if \f$ \vec{p} = p \hat{z} \f$, then 
 * \f[ \epsilon^{\mu}(p,\pm 1) = \mp \frac{1}{\sqrt{2}}(0,1,\pm i,0) \f]
 * 	
 * For higher spin particles (_spin > 1), the polarization vectors are 
 * built according to,
 * \f[\epsilon_{\mu_1\mu_2\ldots\mu_S}(P,M) 
 *     = \sum\limits_{m_{S-1},m_1}((S-1)m_{S-1}1m_1|Sm)
 *       \epsilon_{\mu_1\mu_2\ldots\mu_{S-1}}(P,m_{S-1})
 *       \epsilon_{\mu_S}(P,m_1) \f]
 *
 * where \f$ ((S-1)m_{S-1}1m_1|Sm) \f$ is the Clebsch-Gordon coefficents
 * coupling spin-(S-1) and spin-1 to spin-S (see Clebsch()). This is all done
 * in SetP4().
 *
 * Access to the polarization vectors is provided through the () operator. 
 * Since they are stored as Tensor objects, they can be used with other 
 * Tensor's or inherited classes of Tensor very efficently.
 *
 * <b> Example Usage </b>
 *
 * The amplitude for \f$\omega\rightarrow\pi^+\pi^-\pi^0\f$ can be written as,
 * \f$ \epsilon_{\mu\nu\rho\sigma}p_{\pi^+}^{\mu}p_{\pi^-}^{\nu}
 *      p_{\pi^0}^{\rho}\omega^{\sigma}(p_{\omega},m_{\omega})\f$ where 
 * \f$p_{\pi}\f$ are the pion 4-momenta, \f$\epsilon_{\mu\nu\rho\sigma}\f$ is 
 * the Levi-Civita tensor and \f$\omega^{\sigma}(p_{\omega},m_{\omega})\f$ is 
 * the polarization vector of the \f$\omega\f$ meson with 4-momentum 
 * \f$p_{\omega}\f$ and spin projection \f$m_{\omega}\f$.
 *
 * This can be written in the code as,
 * \include PolVector_omega.ex
 *
 * The spin-S projection operator is defined as,
 *\f[P^{(S)}_{\mu_1\mu_2\ldots\mu_S\nu_1\nu_2\ldots\nu_S}(P) 
 *      = \sum\limits_{M}\epsilon_{\mu_1\mu_2\ldots\mu_S}(P,M)
 *      \epsilon^*_{\nu_1\nu_2\ldots\nu_S}(P,M) \f]
 * for S > 1 and simply as,
 * \f$ -g_{\mu\nu} + \frac{P_{\mu}P_{\nu}}{w^2} \f$
 * where \f$w\f$ is the mass, for a massive spin-1 particle.
 * If the mass = 0, then the projector is set to \f$-g_{\mu\nu}\f$. See 
 * SetProjector() for details. This operator is stored internally in the class
 * for performance purposes. For higher spins, the projection operator 
 * becomes a rather large object that becomes costly to pass by value.
 * PolVector sets _projector when SetP4() is called and returns a constant
 * reference to it via the Projector() member function. 
 *
 * PolVector also provides
 * member functions to return the Feynman and Briet-Wigner propagators 
 * (see Propagator() and PropagatorBW()). Since these aren't stored by the 
 * class, they must be returned by value. Thus, for higher spins this will be 
 * slow. If the projector is going to be fully contracted with other Tensor's
 * to form a scalar (ie. calculating an amplitude), the best approach is to 
 * contract Projector() with the other Tensor's and then divide the result by 
 * the appropriate propagator.
 *
 * <b>Example Usage:</b>
 *
 * The amplitude for \f$ 1^- + 0^- \rightarrow 2^+ \rightarrow 0^- 0^- \f$ can
 * be written as,
 * \f[   p_{1\mu}p_{2\nu} P^{(2)\mu\nu\alpha\beta}(P) q_1^{\rho}q_2^{\delta}
 *     \epsilon^{\sigma}(q_1,m_1)\epsilon_{\rho\delta\sigma\alpha}q_{2\beta}\f]
 * where \f$ p(q) \f$'s are the final(initial) momenta, 
 * \f$ \epsilon^{\sigma}(q_1,m_1) \f$ is the polarization vector for the 
 * \f$ 1^- \f$ particle and 
 * \f$ P^{(2)\mu\nu\alpha\beta}(P) \f$ is the spin-2 projection operator. We 
 * can write this in the code as follows,
 *
 * \include PolVector_amp.ex
 *
 */
//_____________________________________________________________________________

class PolVector {

private:
  // data members:

  vector<Tensor<complex<double > > > _pols; ///< polarization tensors
  Vector4<double> _p4; ///< particle's 4-momentum
  double _mass; ///< particle's mass (pole mass, see SetProjector())
  Tensor<complex<double> > _projector; ///< spin S projection operator
  Spin _spin; ///< particle's spin (S)

  // private functions:

  /** Initialize for spin */
  void _Init(const Spin &__spin);

  /// Copy @a pvect
  void _Copy(const PolVector &__pvect){
    _spin = __pvect._spin;
    _mass = __pvect._mass;
    _p4 = __pvect._p4;
    _projector = __pvect._projector;
    int size = (int)__pvect._pols.size();
    _pols.resize(size);
    for(int i = 0; i < size; i++) _pols[i] = __pvect._pols[i];
  }

  /** Set the Spin projection operator.
   *
   * For spin S = _spin = 1, if _mass > 0 the _projector is set to be 
   * \f$ -g_{\mu\nu} + \frac{P_{\mu}P_{\nu}}{w^2} \f$
   * where \f$w \rightarrow\f$ _mass and \f$P \rightarrow\f$ _p4. <br>
   * If _mass = 0, then _projector is set to \f$-g_{\mu\nu}\f$.
   *
   * For S > 1, _projector is set to be,
   * \f[P^{(S)}_{\mu_1\mu_2\ldots\mu_S\nu_1\nu_2\ldots\nu_S}(P) 
   * = \sum\limits_{M}\epsilon_{\mu_1\mu_2\ldots\mu_S}(P,M)
   * \epsilon^*_{\nu_1\nu_2\ldots\nu_S}(P,M) \f]
   * where \f$M\f$ is the spin projection.
   */
  void _SetProjector();

  /// Boost the pol vectors using boost vector \f$ \vec{\beta}=(bx,by,bz) \f$
  inline void _BoostPolVectors(double __bx,double __by,double __bz){
    int size = (int)(2*_spin + 1);
    for(int i = 0; i < size; i++) _pols[i].Boost(__bx,__by,__bz);  
  }

public:
  
  // create/copy/destroy:

  /** Default Constructor (spin 1) */
  PolVector(){this->_Init(1);}

  /** Ctor */ 
  PolVector(const Spin &__spin){this->_Init(__spin);}

  /// Copy Constructor
  PolVector(const PolVector &__pvect){
    this->_Copy(__pvect);
  }

  /** Destructor */
  virtual ~PolVector(){}

  // Setters:

  /** Sets the PolVector object for a given 4-momentum.  
   *
   * @param p4 particle's 4-momentum
   * @param mass particle's mass (pole mass see SetProjector())   
   *
   * For a particle with spin S = 1 (_spin = 1), there are 2 cases. If 
   * mass > 0, then the rest frame polarization vectors are constructed as,
   * \f[ \epsilon^{\mu}(0,m) = (0,\vec{\epsilon}(m)) \f] where
   * \f[ \vec{\epsilon}(\pm 1) = \mp \frac{1}{\sqrt{2}}(1,\pm i,0),
   * \hspace{0.2in} \vec{\epsilon}(0) = (0,0,1) \f]
   * These are then boosted to the frame where the particle has 4-momentum
   * p4. For the mass = 0 case (the photon), the \f$ M_z = 0 \f$ polarization
   * vector is set to be 0, while the \f$ M_z = \pm 1 \f$ polarization 
   * vectors are set to be the helicity states for the photon. For example,
   * if \f$ \vec{p} = p \hat{z} \f$, then 
   * \f[ \epsilon^{\mu}(p,\pm 1) = \mp \frac{1}{\sqrt{2}}(0,1,\pm i,0) \f]
   * 	   
   *
   * For higher spin particles (_spin > 1), the polarization vectors are 
   * built according to,
   * \f[\epsilon_{\mu_1\mu_2\ldots\mu_S}(P,M) 
   * = \sum\limits_{m_{S-1},m_1}((S-1)m_{S-1}1m_1|Sm)
   *\epsilon_{\mu_1\mu_2\ldots\mu_{S-1}}(P,m_{S-1})\epsilon_{\mu_S}(P,m_1)\f]
   *
   * where \f$ ((S-1)m_{S-1}1m_1|Sm) \f$ are the Clebsch-Gordon coefficents
   * coupling spin-(S-1) and spin-1 to spin-S (see Clebsch()).
   */
  void SetP4(const Vector4<double> &__p4,double __mass);

  // Getters:

  /** Returns the spin projection operator for spin S = _spin.
   *
   * The spin projection operator for spin S is defined as,
   * \f[P^{(S)}_{\mu_1\mu_2\ldots\mu_S\nu_1\nu_2\ldots\nu_S}(P) 
   * = \sum\limits_{M}\epsilon_{\mu_1\mu_2\ldots\mu_S}(P,M)
   * \epsilon^*_{\nu_1\nu_2\ldots\nu_S}(P,M) \f]
   * where \f$P\f$ is the 4-momentum and \f$M \epsilon[-S,S] \f$ is the spin 
   * projection.
   *
   * <b>Example:</b>
   * \f[P^{(2)}_{\mu_1\mu_2\nu_1\nu_2}(P) = \sum\limits_{M}
   * \epsilon_{\mu_1\mu_2}(P,M)\epsilon^*_{\nu_1\nu_2}(P,M) = \frac{1}{2}
   * (\tilde{g}_{\mu_1\nu_1}\tilde{g}_{\mu_2\nu_2}+\tilde{g}_{\mu_1\nu_2}
   * \tilde{g}_{\mu_2\nu_1}) - \frac{1}{3}\tilde{g}_{\mu_1\mu_2}
   * \tilde{g}_{\nu_1\nu_2}\f]
   * where \f$\tilde{g} = g_{\mu\nu} - \frac{P_{\mu}P_{\nu}}{P^2} \f$
   */
  inline const Tensor<complex<double> >& Projector() const {
    return _projector;
  }

  /** Returns the Feynman propagator for spin S = _spin.
   *
   * \returns
   * \f[\frac{P^{(S)}_{\mu_1\mu_2\ldots\mu_S\nu_1\nu_2\ldots\nu_S}(P)}
   * {P^2 - w^2} \f]
   *
   * where \f$w \rightarrow \f$ _mass,\f$P \rightarrow \f$ _p4 and 
   * \f$ P^{(S)}_{\mu_1\mu_2\ldots\mu_S\nu_1\nu_2\ldots\nu_S}(P) \f$ 
   * is returned by Projector().
   *
   * Note: This is slow, for speed use Projector() and do the division after
   *       any matrix/tensor operations.
   */
  Tensor<complex<double> > Propagator() const;

  /** Returns the Breit-Wigner propagator for spin S = _spin.
   *
   * @param width Width of the Breit-Wigner.
   *
   * \returns
   * \f[\frac{P^{(S)}_{\mu_1\mu_2\ldots\mu_S\nu_1\nu_2\ldots\nu_S}(P)}
   * {P^2 - w^2 + iw\Gamma} \f]
   *
   * where \f$\Gamma \rightarrow \f$ width, \f$w \rightarrow \f$ _mass,
   * \f$P \rightarrow \f$ _p4 and 
   * \f$ P^{(S)}_{\mu_1\mu_2\ldots\mu_S\nu_1\nu_2\ldots\nu_S}(P) \f$ is 
   * returned by Projector().
   *
   * Note: This is slow, for speed use Projector() and do the division after
   *       any matrix/tensor operations.
   */
  Tensor<complex<double> > PropagatorBW(double __width = 0.) const;

  /** Returns \f$ \epsilon_{\mu_1\mu_2\ldots\mu_S}(m) \f$.
   *
   * \param m Spin projection along the z-axis (\f$ -S \leq m \leq S \f$) 
   */
  const Tensor<complex<double> >& operator()(const Spin  &__m) const {
    int i = (int)(_spin + __m);
    assert(i >= 0 && i <= (int)(2*_spin));
    return _pols[i];
  }

  /// Returns the spin of the particle
  const Spin& GetSpin() const {
    return _spin;
  }

  /// Returns the mass of the particle
  double GetMass() const {
    return _mass;
  }

  /// Returns the 4-momentum of the particle
  const Vector4<double>& GetP4() const {
    return _p4;
  }

  // Functions:

  /// Clear the PolVector object.
  inline void Clear() {
    if(!_pols.empty()) _pols.clear();
    _spin.SetSpin(0);
  } 

  /// Zeros the Polarization vectors only (_pols[i])
  inline void Zero() {
    int size = (int)(2*_spin + 1);
    for(int i = 0; i < size; i++) _pols[i].Zero();
  }

  /// Boost the particle using boost vector \f$ \vec{\beta}=(bx,by,bz) \f$
  inline void Boost(double __bx,double __by,double __bz){
    int size = (int)(2*_spin + 1);
    for(int i = 0; i < size; i++) _pols[i].Boost(__bx,__by,__bz);  
    _p4.Boost(__bx,__by,__bz);
    this->_SetProjector();
  }

  /** Boost the particle to p4's rest frame.
   * @param p4 a valid 4-momentum (must have \f$\beta < 1 \f$ 
   */
  inline void Boost(const Vector4<double> &__p4){
    int size = (int)(2*_spin + 1);
    for(int i = 0; i < size; i++) _pols[i].Boost(__p4);  
    _p4.Boost(__p4);
    this->_SetProjector();
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
			double __mass,Tensor<complex<double> > &__projector);
};
//_____________________________________________________________________________

#endif /* _PolVector_H */
