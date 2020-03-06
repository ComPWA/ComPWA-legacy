// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef DATA_DATASET_HPP_
#define DATA_DATASET_HPP_

#include <memory>

#include "Core/Event.hpp"
#include "Core/Function.hpp"

namespace ComPWA {
class Kinematics;
namespace Data {

struct DataSet {
  ComPWA::DataMap Data;
  std::vector<double> Weights;
};

inline void resize(DataSet &set, size_t size) {
  set.Weights.resize(size);
  for (auto &i : set.Data) {
    i.second.resize(size);
  }
}

ComPWA::EventCollection
addIntensityWeights(std::shared_ptr<ComPWA::Intensity> Intensity,
                    const ComPWA::EventCollection &DataSample,
                    const ComPWA::Kinematics &Kinematics);

} // namespace Data
} // namespace ComPWA

#endif
