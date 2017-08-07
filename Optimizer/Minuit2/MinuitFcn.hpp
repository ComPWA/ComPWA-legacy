//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Minuit2 function to be optimized.
/*! \class MinuitFcn
 * @file MinuitFcn.hpp
 * Based on the Minuit2 FcnBase. This class uses the ControlParameter interface
 * for the
 * optimization.
 */

#ifndef _OIFMinuitFcn_HPP
#define _OIFMinuitFcn_HPP

#include <vector>
#include <memory>
#include <string>
#include <map>

#include "Core/Estimator.hpp"
#include "Core/ParameterList.hpp"

#include "Minuit2/FCNBase.h"

namespace ROOT {
namespace Minuit2 {

class MinuitFcn : public FCNBase {

public:
  MinuitFcn(std::shared_ptr<ComPWA::IEstimator> theData,
            ComPWA::ParameterList &parList);
  virtual ~MinuitFcn();

  double operator()(const std::vector<double> &x) const;

  double Up() const;

  inline void setNameID(const unsigned int id, const std::string &name) {
    auto result =
        _parNames.insert(std::pair<unsigned int, std::string>(id, name));
    if (!result.second) {
      std::stringstream ss;
      ss << "MinuitFcn::setNameID(): Could not create entry in ID-name map for "
            "id="
         << id << " and name=" << name;
      throw std::runtime_error(ss.str());
    }
  };

  inline std::string parName(const unsigned int id) {
    return _parNames.at(id);
  };

private:
  //// pointer to the ControlParameter (e.g. Estimator)
  std::shared_ptr<ComPWA::IEstimator> _myDataPtr;
  
  /// List of Parameters the ControlParameter needs
  ComPWA::ParameterList &_parList;
  
  /// mapping of minuit ids to ComPWA names 
  std::map<unsigned int, std::string> _parNames;
};

} // namespace Minuit2
} // namespace ROOT

#endif
