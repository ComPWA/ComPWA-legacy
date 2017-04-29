//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//----------------------------------------------------------------------------------

#ifndef PARTIALDECAY_HPP_
#define PARTIALDECAY_HPP_

#include <memory>
#include <boost/property_tree/ptree.hpp>

#include "Physics/Resonance.hpp"
#include "Physics/HelicityFormalism/AbstractDynamicalFunction.hpp"
#include "Physics/HelicityFormalism/AmpWignerD.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

class PartialDecay : public ComPWA::Physics::Resonance {

public:
  PartialDecay(){};

  //! Clone function
  virtual PartialDecay *Clone(std::string newName = "") const {
    auto tmp = new PartialDecay(*this);
    tmp->SetName(newName);
    return tmp;
  }

  //! Implementation of interface for streaming info about the strategy
  virtual std::string to_str() const { return std::string("PartialDecay"); }

  virtual void GetParameters(ParameterList &list);

  /**! Evaluate decay */
  std::complex<double> Evaluate(const dataPoint &point) const {
    std::complex<double> result =
        std::polar(_magnitude->GetValue(), _phase->GetValue());
    result *= _angD->Evaluate(point, _dataPos + 1, _dataPos + 2);
    result *= _dynamic->Evaluate(point);

    return result;
  };

  //! Check of tree is available
  virtual bool HasTree() const { return true; }

  /**! Setup function tree */
  virtual std::shared_ptr<FunctionTree>
  GetTree(ParameterList &sample, ParameterList &toySample, std::string suffix);

  /**
   Factory for PartialDecay

   @param pt Configuration tree
   @return Constructed object
   */
  static std::shared_ptr<ComPWA::Physics::Resonance>
  Factory(const boost::property_tree::ptree &pt);

  static boost::property_tree::ptree Save(std::shared_ptr<Resonance> obj);

  /**
   Get WignerD function

   @return Shared_ptr<AmpWignerD>
   */
  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> GetWignerD() {
    return _angD;
  }

  /**
   Set WignerD function

   @param w WignerD function
   */
  void SetWignerD(
      std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> w) {
    _angD = w;
  }

  /**
   Get dynamical function (e.g. Breit-Wigner parametrization)

   @return Shared_ptr<AbstractDynamicalFunction>
   */
  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AbstractDynamicalFunction>
  GetDynamicalFunction() {
    return _dynamic;
  }

  /**
   Set dynamical function

   @param f Dynamical function
   */
  void SetDynamicalFunction(
      std::shared_ptr<
          ComPWA::Physics::HelicityFormalism::AbstractDynamicalFunction>
          f) {
    _dynamic = f;
  }

  //! Set position of variables within dataPoint
  void SetDataPosition(int pos) { _dataPos = pos; }

  //! Get position of variables within dataPoint
  int GetDataPosition() const { return _dataPos; }

  //! Set position of variables within dataPoint
  void SetSubSystem(SubSystem sys) {
    _subSystem = sys;
    _dataPos = 3*dynamic_cast<HelicityKinematics *>(Kinematics::Instance())
                   ->GetDataID(_subSystem);
    if (_dynamic) {
      auto invMassLimit =
          dynamic_cast<HelicityKinematics *>(Kinematics::Instance())
              ->GetInvMassBounds(_subSystem);
      _dynamic->SetLimits(invMassLimit);
      _dynamic->SetDataPosition(_dataPos);
    } else {
      LOG(error) << "PartialDecay::SetSubSystem() | Dynamic function not set "
                    "yet so we can not set limits and data position.";
    }
  }

  //! Get position of variables within dataPoint
  SubSystem GetSubSystem() const { return _subSystem; }

  /*! Set phase space sample
   * We use the phase space sample to calculate the normalization. The sample
   * should be without efficiency applied.
   */
  virtual void
  SetPhspSample(std::shared_ptr<std::vector<ComPWA::dataPoint>> phspSample) {
    _dynamic->SetPhspSample(phspSample);
  };

protected:
  /**! Position where variables are stored in dataPoint
   * We expect to find the invariant mass of the system at @param _dataPos,
   * cosTheta at @param _dataPos+1 and phi at @param _dataPos+2 */
  int _dataPos;

  SubSystem _subSystem;

  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AmpWignerD> _angD;
  std::shared_ptr<ComPWA::Physics::HelicityFormalism::AbstractDynamicalFunction>
      _dynamic;

  // TODO: Operator* to construct SequentialTwoBodyDecay
};

} /* namespace HelicityFormalism */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PARTIALDECAY_ */
