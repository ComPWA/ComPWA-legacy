// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "DataSet.hpp"
#include "Core/Kinematics.hpp"
#include "Core/Logging.hpp"

namespace ComPWA {
namespace Data {

std::vector<Event> reduceToPhaseSpace(const std::vector<Event> &Events,
                                      const ComPWA::Kinematics &Kinematics) {
  std::vector<Event> tmp;
  LOG(INFO) << "DataSet::reduceToPhaseSpace(): "
               "Remove all events outside PHSP boundary from data sample.";

  std::copy_if(Events.begin(), Events.end(), std::back_inserter(tmp),
               [&](const Event &evt) {
                 DataPoint point = Kinematics.convert(evt);
                 return Kinematics.isWithinPhaseSpace(point);
               });
  LOG(INFO) << "reduceToPhaseSpace(): Removed " << Events.size() - tmp.size()
            << " from " << Events.size() << " Events ("
            << (1.0 - tmp.size() / Events.size()) * 100 << "%).";
  return tmp;
}

std::vector<Event>
addIntensityWeights(std::shared_ptr<ComPWA::Intensity> Intensity,
                    const std::vector<Event> &Events,
                    const ComPWA::Kinematics &Kinematics) {
  auto dataset = convertEventsToDataSet(Events, Kinematics);
  auto weights = Intensity->evaluate(dataset.Data);
  std::vector<Event> NewEvents;
  NewEvents.reserve(Events.size());
  for (size_t i = 0; i < Events.size(); ++i) {
    Event evt(Events[i]);
    evt.Weight *= weights[i];
  }
  return NewEvents;
}

DataSet convertEventsToDataSet(const std::vector<Event> &Events,
                               const ComPWA::Kinematics &Kinematics) {
  DataSet dataset;
  dataset.VariableNames = Kinematics.getKinematicVariableNames();
  for (auto x : dataset.VariableNames) {
    dataset.Data.push_back(std::vector<double>());
  }
  for (auto const &evt : Events) {
    DataPoint p = Kinematics.convert(evt);
    for (size_t i = 0; i < p.KinematicVariableList.size(); ++i) {
      dataset.Data[i].push_back(p.KinematicVariableList[i]);
      // TODO: just do it with iterators and just move both iterators forward on
      // each step!
    }
    dataset.Weights.push_back(evt.Weight);
  }
  return dataset;
}

} // namespace Data
} // namespace ComPWA
