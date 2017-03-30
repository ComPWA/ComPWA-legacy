 
            
          
        
      
    
  


//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//		Peter Weidenkaff -
//-------------------------------------------------------------------------------

#ifndef KINEMATICS_HPP_
#define KINEMATICS_HPP_

#include <vector>
#include <string>
#include <complex>

#include "Core/Event.hpp"
#include "Core/PhysConst.hpp"
#include "Core/Spin.hpp"

namespace ComPWA {

class dataPoint;

static const char *formFactorTypeString[] = {"noFormFactor", "BlattWeisskopf",
                                             "CrystalBarrel"};

enum formFactorType { noFormFactor = 0, BlattWeisskopf = 1, CrystalBarrel = 2 };

class Kinematics {
public:
  //! singleton pattern
  static Kinematics *Instance();

  //! converts Event to dataPoint
  virtual void EventToDataPoint(const ComPWA::Event &ev,
                                dataPoint &point) const = 0;

  //! Checks of data point is within phase space boundaries
  virtual bool IsWithinPhsp(const dataPoint &point) const = 0;

  //! calculated the PHSP volume of the current decay by MC integration
  virtual double GetPhspVolume();

  //! Get number of variables
  virtual unsigned int GetNVars() const { return _varNames.size(); }

  //! Get final state
  virtual std::vector<int> GetFinalState() { return _finalState; }

  //! Get inital state
  virtual std::vector<int> GetInitialState() { return _initialState; }

  /** Calculate Break-up momentum squared
   *
   * Calculate Break-up momentum at energy @param sqrtS for particles with
   * masses @param ma and @param mb .
   * From PDG2014 Eq.46-20a. Below threshold the function is analytically
   * continued.
   * @param sqrtS center-of-mass energy
   * @param ma mass particle A
   * @param mb mass particle B
   * @return |break-up momentum|
   */
  static double qSqValue(double sqrtS, double ma, double mb);

  /** Calculate Break-up momentum
   *
   * Calculate Break-up momentum at energy @param sqrtS for particles with
   * masses @param ma and @param mb .
   * From PDG2014 Eq.46-20a. Below threshold the function is analytically
   * continued.
   * @param sqrtS center-of-mass energy
   * @param ma mass particle A
   * @param mb mass particle B
   * @return |break-up momentum|
   */
  static std::complex<double> qValue(double sqrtS, double ma, double mb);
  /** Two body phsp factor
   *
   * From PDG2014 Eqn.47-2
   * @param sqrtS invariant mass of particles A and B
   * @param ma Mass of particle A
   * @param mb Mass of particle B
   * @return
   */
  static std::complex<double> phspFactor(double sqrtS, double ma, double mb);

protected:
  std::vector<int> _initialState;
  std::vector<int> _finalState;

  //! Internal names of variabes
  std::vector<std::string> _varNames;
  //! Latex titles for variables
  std::vector<std::string> _varTitles;

  // Singleton stuff
  static Kinematics *_inst;

  //! Constructor
  Kinematics(std::vector<int> initial = std::vector<int>(),
             std::vector<int> finalS = std::vector<int>())
      : _initialState(initial), _finalState(finalS),
        is_PS_area_calculated_(false), PS_area_(0.0){

                                       };

  //! Delete Copy constructor
  Kinematics(const Kinematics &) = delete;

  //! Default destructor
  virtual ~Kinematics(){};

  //! Delete assignment operator
  void operator=(const Kinematics &) = delete;

  virtual double calculatePSArea() = 0;
  bool is_PS_area_calculated_;
  double PS_area_;
};

} /* namespace ComPWA */
#endif /* KINEMATICS_HPP_ */
