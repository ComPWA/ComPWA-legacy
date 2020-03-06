// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "DataSet.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Logging.hpp"

namespace ComPWA {
namespace Data {

ComPWA::EventCollection
addIntensityWeights(std::shared_ptr<ComPWA::Intensity> Intensity,
                    const EventCollection &DataSample,
                    const ComPWA::Kinematics &Kinematics) {
  auto DataSet = Kinematics.convert(DataSample);
  auto Weights = Intensity->evaluate(DataSet.Data);
  ComPWA::EventCollection NewEventList{DataSample.Pids};
  std::transform(DataSample.Events.begin(), DataSample.Events.end(),
                 Weights.begin(), NewEventList.Events.begin(),
                 [](Event Event, double Weight) {
                   Event.Weight *= Weight;
                   return Event;
                 });
  return NewEventList;
}

} // namespace Data
} // namespace ComPWA
