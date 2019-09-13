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
  std::vector<Event> NewEvents(Events.size());
  std::transform(Events.begin(), Events.end(), weights.begin(),
                 NewEvents.begin(), [](Event evt, double weight) {
                   evt.Weight *= weight;
                   return evt;
                 });
  return NewEvents;
}

DataSet convertEventsToDataSet(std::vector<Event>::const_iterator EventsBegin,
                               std::vector<Event>::const_iterator EventsEnd,
                               const ComPWA::Kinematics &Kinematics) {
  auto VariableNames = Kinematics.getKinematicVariableNames();
  std::vector<std::vector<double>> Data(VariableNames.size());

  std::vector<double> Weights;
  for (auto evt = EventsBegin; evt != EventsEnd; ++evt) {
    DataPoint p = Kinematics.convert(*evt);
    auto data_it = Data.begin();
    for (auto kinvar : p.KinematicVariableList) {
      data_it->push_back(kinvar); // warning: past the end access possible
      ++data_it;
    }
    Weights.push_back(evt->Weight);
  }

  return DataSet{
      .Data = Data, .Weights = Weights, .VariableNames = VariableNames};
}

DataSet convertEventsToDataSet(const std::vector<Event> &Events,
                               const ComPWA::Kinematics &Kinematics) {
  return convertEventsToDataSet(Events.begin(), Events.end(), Kinematics);
}

std::vector<DataPoint> convertEventsToDataPoints(const std::vector<Event> &Events,
                               const ComPWA::Kinematics &Kinematics) {
  std::vector<DataPoint> points(Events.size());
  std::transform(
      Events.begin(), Events.end(), points.begin(),
      [&Kinematics](const Event &ev) { return Kinematics.convert(ev); });

  return points;
}


} // namespace Data
} // namespace ComPWA
