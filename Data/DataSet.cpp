// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "DataSet.hpp"
#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Data {

DataSet::DataSet(const std::vector<Event> &Events) : EventList(Events) {}

DataSet::DataSet(const std::vector<DataPoint> &DataPoints)
    : DataPointList(DataPoints) {
  convertDataPointsToParameterList();
}

void DataSet::reduceToPhaseSpace(
    std::shared_ptr<ComPWA::Kinematics> Kinematics) {
  std::vector<Event> tmp;
  LOG(INFO) << "DataSet::reduceToPhaseSpace(): "
               "Remove all events outside PHSP boundary from data sample.";

  auto newend =
      std::remove_if(EventList.begin(), EventList.end(), [&](const Event &evt) {
        DataPoint point = Kinematics->convert(evt);
        if (Kinematics->isWithinPhaseSpace(point))
          return true;
        return false;
      });
  unsigned int PreviousSize(EventList.size());
  // remove_if only reorders to elements to remove to a the back
  // erase actually removes them
  EventList.erase(newend, EventList.end());
  unsigned int NewSize(EventList.size());
  LOG(INFO) << "DataSet::reduceToPhsp(): " << EventList.size() << " from "
            << EventList.size() << "(" << (1.0 - PreviousSize / NewSize) * 100
            << "%) were removed.";
}

void DataSet::addIntensityWeights(
    std::shared_ptr<ComPWA::Intensity> Intensity,
    std::shared_ptr<ComPWA::Kinematics> Kinematics) {
  convertEventsToParameterList(Kinematics);
  std::vector<std::vector<double>> data;
  for (auto x : HorizontalDataList.mDoubleValues()) {
    data.push_back(x->value());
  }
  auto weights = Intensity->evaluate(data);
  for (size_t i = 0; i < EventList.size(); ++i) {
    EventList[i].Weight *= weights[i];
  }
}

void DataSet::convertEventsToDataPoints(
    std::shared_ptr<ComPWA::Kinematics> Kinematics) {
  auto VarNames = Kinematics->getKinematicVariableNames();
  if (VarNames == KinematicVariableNames &&
      DataPointList.size() == EventList.size()) {
    // nothing has changed, the cached values are fine
    return;
  }
  if (VarNames != KinematicVariableNames) {
    if (0 < KinematicVariableNames.size())
      LOG(INFO) << "DataSet::convertEventsToDataPoints(): the kinematic "
                   "variables have changed! recalculating...";
  } else if (DataPointList.size() != EventList.size()) {
    LOG(INFO) << "DataSet::convertEventsToDataPoints(): the event list "
                 "size has changed! recalculating...";
  }

  KinematicVariableNames = VarNames;
  DataPointList.clear();
  DataPointList.reserve(EventList.size());
  for (auto const &evt : EventList) {
    DataPointList.push_back(Kinematics->convert(evt));
  }
}

void DataSet::convertEventsToParameterList(
    std::shared_ptr<ComPWA::Kinematics> Kinematics) {
  convertEventsToDataPoints(Kinematics);
  convertDataPointsToParameterList();
}

const std::vector<Event> &DataSet::getEventList() const { return EventList; }
const std::vector<DataPoint> &DataSet::getDataPointList() const {
  return DataPointList;
}
const ParameterList &DataSet::getParameterList() const {
  return HorizontalDataList;
}

const std::vector<std::string> &DataSet::getKinematicVariableNames() const {
  return KinematicVariableNames;
}

void DataSet::convertDataPointsToParameterList() {
  if (0 < DataPointList.size()) {
    // reset old parameter list
    HorizontalDataList = ParameterList();
    size_t NumberOfKinematicVariables =
        DataPointList[0].KinematicVariableList.size();
    std::vector<std::vector<double>> Data(NumberOfKinematicVariables,
                                          std::vector<double>());
    std::vector<double> Weights;
    Weights.reserve(DataPointList.size());

    for (auto const &point : DataPointList) {
      Weights.push_back(point.Weight);
      for (unsigned int i = 0; i < NumberOfKinematicVariables; ++i)
        Data[i].push_back(point.KinematicVariableList[i]);
    }

    // Add data vector to ParameterList
    for (auto x : Data)
      HorizontalDataList.addValue(MDouble("", x));
    // Adding weight at the end
    HorizontalDataList.addValue(MDouble("Weight", Weights));
  }
}

void DataSet::convertParameterListToDataPoints() {
  if (0 < HorizontalDataList.mDoubleValues().size() - 1) {
    unsigned int NumberOfKinematicVariables =
        HorizontalDataList.mDoubleValues().size() - 1;
    for (unsigned int evt_index = 0;
         evt_index < HorizontalDataList.mDoubleValues()[0]->values().size();
         ++evt_index) {
      DataPoint dp;
      for (unsigned int j = 0; j < NumberOfKinematicVariables; ++j)
        dp.KinematicVariableList.push_back(
            HorizontalDataList.mDoubleValue(j)->values()[evt_index]);

      dp.Weight = HorizontalDataList.mDoubleValue(NumberOfKinematicVariables)
                      ->values()[evt_index];
      DataPointList.push_back(dp);
    }
  }
}

} // namespace Data
} // namespace ComPWA
