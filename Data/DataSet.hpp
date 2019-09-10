// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef DATA_DATASET_HPP_
#define DATA_DATASET_HPP_

#include <memory>
#include <vector>

#include "Core/Event.hpp"
#include "Core/Function.hpp"

namespace ComPWA {
class Kinematics;
namespace Data {

using DataList = std::vector<std::vector<double>>;

struct DataSet {
  DataList Data;
  std::vector<double> Weights;
  std::vector<std::string> VariableNames;
};

inline void resize(DataSet &set, size_t size) {
  set.Weights.resize(size);
  for (auto &i : set.Data) {
    i.resize(size);
  }
}

std::vector<Event> reduceToPhaseSpace(const std::vector<Event> &Events,
                                      const ComPWA::Kinematics &Kinematics);

std::vector<Event>
addIntensityWeights(std::shared_ptr<ComPWA::Intensity> Intensity,
                    const std::vector<Event> &Events,
                    const ComPWA::Kinematics &Kinematics);

DataSet convertEventsToDataSet(const std::vector<Event> &Events,
                               const ComPWA::Kinematics &Kinematics);

} // namespace Data
} // namespace ComPWA

#endif
