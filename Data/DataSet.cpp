// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "DataSet.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Logging.hpp"

namespace ComPWA {
namespace Data {

std::vector<Event>
addIntensityWeights(std::shared_ptr<ComPWA::Intensity> Intensity,
                    const std::vector<Event> &Events,
                    const ComPWA::Kinematics &Kinematics) {
  auto dataset = Kinematics.convert(Events);
  auto weights = Intensity->evaluate(dataset.Data);
  std::vector<Event> NewEvents(Events.size());
  std::transform(Events.begin(), Events.end(), weights.begin(),
                 NewEvents.begin(), [](Event evt, double weight) {
                   evt.Weight *= weight;
                   return evt;
                 });
  return NewEvents;
}

} // namespace Data
} // namespace ComPWA
