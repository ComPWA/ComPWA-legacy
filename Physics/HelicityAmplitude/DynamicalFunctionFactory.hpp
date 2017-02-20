//-------------------------------------------------------------------------------
// Copyright (c) 2013 Stefan Pflueger.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//   Stefan Pflueger - initial API and implementation
//-------------------------------------------------------------------------------

#ifndef PHYSICS_DYNAMICALDECAYFUNCTIONS_DYNAMICALFUNCTIONFACTORY_HPP_
#define PHYSICS_DYNAMICALDECAYFUNCTIONS_DYNAMICALFUNCTIONFACTORY_HPP_

#include <memory>

#include <boost/assign.hpp>

#include "Core/PhysConst.hpp"

#include "Physics/DynamicalDecayFunctions/AbstractDynamicalFunction.hpp"
#include "Physics/HelicityAmplitude/ParticleStateDefinitions.hpp"

namespace ComPWA {
namespace Physics {
namespace HelicityFormalism {

enum class DynamicalInfoTypes {
  TOP_NODE, RELATIVE_BREIT_WIGNER
};

struct ParticleStateComperator {
    bool operator()(const ParticleStateInfo& a, const ParticleStateInfo& b) const {
      if (a.pid_information_ < b.pid_information_)
        return true;
      else if (a.pid_information_ > b.pid_information_)
        return false;
      return false;
    }
};

const boost::unordered_map<DynamicalInfoTypes, std::string> DynamicalTypeToString =
    boost::assign::map_list_of(DynamicalInfoTypes::TOP_NODE, "topNode")(
        DynamicalInfoTypes::RELATIVE_BREIT_WIGNER, "relBW");

const boost::unordered_map<std::string, DynamicalInfoTypes> StringToDynamicalType =
    boost::assign::map_list_of("topNode", DynamicalInfoTypes::TOP_NODE)("relBW",
        DynamicalInfoTypes::RELATIVE_BREIT_WIGNER);


class DynamicalFunctionFactory
{
  std::map<ParticleStateInfo,
      std::shared_ptr<AbstractDynamicalFunction>,
	  ParticleStateComperator> dynamical_function_list_;

  std::shared_ptr<AbstractDynamicalFunction> generateRelativisiticBreitWigner(
      const ParticleStateInfo& state_info,
      const ExternalParameters& external_parameters);

public:
  DynamicalFunctionFactory();
  virtual ~DynamicalFunctionFactory();

  std::shared_ptr<AbstractDynamicalFunction> generateDynamicalFunction(
      const ParticleStateInfo& state_info,
      const ExternalParameters& external_parameters);
};

} /* namespace DynamicalFunctions */
} /* namespace Physics */
} /* namespace ComPWA */

#endif /* PHYSICS_DYNAMICALDECAYFUNCTIONS_DYNAMICALFUNCTIONFACTORY_HPP_ */
