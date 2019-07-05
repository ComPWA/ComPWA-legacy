// Copyright (c) 2019 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//===----------------------------------------------------------------------===//
///
/// \file
/// This file contains functions which will update parameters' value, range etc.
/// of a property tree according providing values/FitParameters.
///
//===----------------------------------------------------------------------===//
#ifndef COMPWA_TOOLS_UPDATE_PTREE_PARAMETER_HPP_
#define COMPWA_TOOLS_UPDATE_PTREE_PARAMETER_HPP_

#include "Core/FitParameter.hpp"

#include <boost/property_tree/ptree.hpp>

#include <string>
#include <vector>

namespace ComPWA {
namespace Tools {

/// Update range of specified parameters of a ptree.
/// \param Tree Property tree which contains parameters.
/// \param ParameterType Parameters with same type
/// (e.g., Magnitude) will be updated.
/// \param Min Minimum value of the range.
/// \param Max Maximum value of the range.
void updateParameterRangeByType(boost::property_tree::ptree &Tree,
                                const std::string &ParameterType, double Min,
                                double Max);

/// Update range of specified parameters of a ptree.
/// \param Tree Property tree which contains parameters.
/// \param ParameterName Parameters with same name will be updated.
/// \param Min Minimum value of the range.
/// \param Max Maximum value of the range.
void updateParameterRangeByName(boost::property_tree::ptree &Tree,
                                const std::string &ParameterName, double Min,
                                double Max);

/// Update value of specified parameters of a ptree.
/// \param Tree Property tree which contains parameters.
/// \param ParameterName Parameters with same name will be updated.
/// \param Value New value of parameter.
void updateParameterValue(boost::property_tree::ptree &Tree,
                          const std::string ParameterName, double Value);

/// Fix specified parameters of a ptree.
/// \param Tree Property tree which contains parameters.
/// \param ParameterName Parameters with same name will be updated.
/// \param Value New value of parameter. Default means fix to current value.
void fixParameter(boost::property_tree::ptree &Tree,
                  const std::string ParameterName, double Value = -999);

/// Release specified parameters of a ptree.
/// \param Tree Property tree which contains parameters.
/// \param ParameterName Parameters with same name will be updated.
/// \param Value New value of parameter. Default means the current value.
void releaseParameter(boost::property_tree::ptree &Tree,
                      const std::string ParameterName, double Value = -999);

/// Update specified parameters of a ptree.
/// \param Tree Property tree which contains parameters.
/// \param KeyType The type(e.g., "Name" or "Type") of 'key' used to find
/// parameter in the ptree.
/// \param KeyValue The value of 'key' used to find parameter in the ptree.
/// \param Value New value of parameter.
/// \param Fix Parameter will be fixed or not.
/// \param Min Lower range of the parameter
/// \param Max Upper range of the parameter
/// \param UpdateValue If update value of the parameter or not.
/// \param UpdateFix If update fix status of the parameter or not.
/// \param UpdateRange If update range of the parameter or not.
void updateParameter(boost::property_tree::ptree &Tree,
                     const std::string &KeyType, const std::string &KeyValue,
                     double Value, bool Fix, double Min, double Max,
                     bool UpdateValue, bool UpdateFix, bool UpdateRange);

/// Update value, range, fix status of parameters of a ptree, the new values
/// comes from \p fitParameters.
/// \param Tree Property tree which contains parameters.
/// \param FitParameters Target parameters' vector. Parameters appear both in
/// \p Tree and \p FitParameters will be updated according the parameter in
/// \p FitParameters
void updateParameter(boost::property_tree::ptree &Tree,
                     const FitParameterList &FitParameters);

} // end namespace Tools
} // end namespace ComPWA

#endif
