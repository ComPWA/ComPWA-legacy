// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.
#include "DataTransformation.hpp"
#include "Core/Kinematics.hpp"

namespace ComPWA {
namespace Data {

void reduceToPhaseSpace(std::vector<Event> &EventList,
                        std::shared_ptr<ComPWA::Kinematics> Kinematics) {
  std::vector<Event> tmp;
  LOG(INFO) << "DataConverter::reduceToPhaseSpace() | "
               "Remove all events outside PHSP boundary from data sample.";

  auto newend =
      std::remove_if(EventList.begin(), EventList.end(), [&](const Event &evt) {
        try {
          DataPoint point = Kinematics->convert(evt);
        } catch (BeyondPhsp &ex) { // event outside phase, remove
          return true;
        }
        return false;
      });
  unsigned int PreviousSize(EventList.size());
  // remove_if only reorders to elements to remove to a the back
  // erase actually removes them
  EventList.erase(newend, EventList.end());
  unsigned int NewSize(EventList.size());
  LOG(INFO) << "Data::reduceToPhsp() | " << EventList.size() << " from "
            << EventList.size() << "(" << (1.0 - PreviousSize / NewSize) * 100
            << "%) were removed.";
}

std::vector<DataPoint>
convertEventsToDataPoints(const std::vector<Event> &Events,
                          std::shared_ptr<ComPWA::Kinematics> Kinematics) {
  std::vector<DataPoint> DataPoints;
  for (auto const &evt : Events) {
    try {
      DataPoints.push_back(Kinematics->convert(evt));
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }
  }
  return DataPoints;
}

ParameterList
convertEventsToParameterList(const std::vector<Event> &Events,
                             std::shared_ptr<ComPWA::Kinematics> Kinematics) {
  std::vector<DataPoint> DataPoints;
  for (auto const &evt : Events) {
    try {
      DataPoints.push_back(Kinematics->convert(evt));
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }
  }

  return convertDataPointsToParameterList(DataPoints);
}

ComPWA::ParameterList
convertDataPointsToParameterList(const std::vector<DataPoint> &DataPointList) {
  ComPWA::ParameterList DataList;
  if (0 < DataPointList.size()) {
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
      DataList.addValue(MDouble("", x));
    // Adding weight at the end
    DataList.addValue(MDouble("Weight", Weights));
  }
  return DataList;
}

std::vector<DataPoint> convertParameterListToDataPoints(
    const ComPWA::ParameterList &DataParameterList) {
  std::vector<DataPoint> DataPoints;
  if (0 < DataParameterList.mDoubleValues().size() - 2) {
    unsigned int NumberOfKinematicVariables =
        DataParameterList.mDoubleValues().size() - 2;
    for (unsigned int evt_index = 0;
         evt_index < DataParameterList.mDoubleValues()[0]->values().size();
         ++evt_index) {
      DataPoint dp;
      for (unsigned int j = 0; j < NumberOfKinematicVariables; ++j)
        dp.KinematicVariableList.push_back(
            DataParameterList.mDoubleValue(j)->values()[evt_index]);

      dp.Weight =
          DataParameterList.mDoubleValue(NumberOfKinematicVariables)
              ->values()[evt_index];
      DataPoints.push_back(dp);
    }
  }
  return DataPoints;
}

/*
std::shared_ptr<Data>
rndSubSet(std::shared_ptr<Kinematics> kin, unsigned int size,
                         std::shared_ptr<ComPWA::Tools::Generator> gen) {
  std::shared_ptr<Data> out(new Data());
  rndReduceSet(kin, size, gen, this, out.get());
  return out;
}

void rndReduceSet(
                                 unsigned int size,
                                 std::shared_ptr<ComPWA::Generator> gen,
                                 Data *in1, Data *out1, Data *in2, Data *out2) {
  if (!in1)
    throw std::runtime_error("rndSubSet() | No input data set!");
  if (!out1)
    throw std::runtime_error("rndSubSet() | No output data set!");
  if (out1->numEvents())
    throw std::runtime_error("rndSubSet() | First output sample not empty!");
  if (in2) {
    if (in1->numEvents() != in2->numEvents())
      throw std::runtime_error(
          "rndSubSet() | Samples have different event count!");
    if (!out2)
      throw std::runtime_error("rndSubSet() | Second output set is NULL!");
    if (out2->numEvents())
      throw std::runtime_error("rndSubSet() | Second output sample not empty!");
  }

  unsigned int totalSize = in1->numEvents();
  unsigned int newSize = totalSize;

  for (unsigned int i = 0; i < totalSize;
       i++) { // count how many events are not within PHSP
    ComPWA::DataPoint point;
    try {
      kin->convert(in1->event(i), point);
    } catch (ComPWA::BeyondPhsp &ex) { // event outside phase, remove
      newSize--;
      continue;
    }
    //    dataPoint point(in1->getEvent(i));
    //    if(!Kinematics::instance()->isWithinPhsp(point)) newSize--;
  }
  double threshold = (double)size / newSize; // calculate threshold
  for (unsigned int i = 0; i < totalSize; i++) {
    ComPWA::DataPoint point;
    try {
      kin->convert(in1->event(i), point);
    } catch (ComPWA::BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }
    //    dataPoint point(in1->getEvent(i)); //use first sample for
    // hit&miss
    //    if(!Kinematics::instance()->isWithinPhsp(point)) continue;
    if (gen->uniform(0, 1) < threshold) {
      out1->add(in1->event(i));
      // write second sample if event from first sample were accepted
      if (in2)
        out2->add(in2->event(i));
    }
  }
  if (out2)
    assert(out1->numEvents() == out2->numEvents());

  LOG(DEBUG) << "DataReader::rndReduceSet() | sample size reduced to "
             << out1->numEvents();
  return;
}*/

} // namespace Data
} // namespace ComPWA
