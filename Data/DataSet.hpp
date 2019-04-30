// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef DATA_DATASET_HPP_
#define DATA_DATASET_HPP_

#include <string>
#include <vector>

#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"

namespace ComPWA {
class Kinematics;
class Intensity;
namespace Data {

class DataSet {
public:
  virtual ~DataSet() = default;
  DataSet(const std::vector<Event> &Events);
  DataSet(const std::vector<DataPoint> &DataPoints);

  void reduceToPhaseSpace(std::shared_ptr<ComPWA::Kinematics> Kinematics);
  void addIntensityWeights(std::shared_ptr<ComPWA::Intensity> Intensity,
                           std::shared_ptr<ComPWA::Kinematics> Kinematics);

  void convertEventsToDataPoints(std::shared_ptr<Kinematics> Kinematics);
  void convertEventsToParameterList(std::shared_ptr<Kinematics> Kinematics);

  const std::vector<Event> &getEventList() const;
  const std::vector<DataPoint> &getDataPointList() const;
  const ParameterList &getParameterList() const;

  const std::vector<std::string> &getKinematicVariableNames() const;

private:
  void convertDataPointsToParameterList();
  void convertParameterListToDataPoints();

  std::vector<Event> EventList;
  std::vector<DataPoint> DataPointList;

  /// A 'horizontal' list of the kinematic variables. For each variable
  /// (e.g. m23sq, m13sq ...) a MultiDouble is added to ParameterList.
  /// This ParameterList is used in the FunctionTree.
  ParameterList HorizontalDataList;

  std::vector<std::string> KinematicVariableNames;
};

} // namespace Data
} // namespace ComPWA

#endif
