// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_TOOLS_PLOTTING_ROOTPLOTDATA_HPP_
#define COMPWA_TOOLS_PLOTTING_ROOTPLOTDATA_HPP_

#include <map>

#include "Core/FitParameter.hpp"
#include "Core/Function.hpp"

#include "TFile.h"

namespace ComPWA {
namespace Data {
struct DataSet;
}
namespace Physics {
class ParticleStateTransitionKinematicsInfo;
}
namespace Tools {
namespace Plotting {

///
/// \class RootPlotData
/// Allows output of a data sample and an Intensity (and optionally its
/// components) into a ROOT file via TTrees. See the appropriate write
/// functions. The Intensity is evaluated using a phase space sample, which
/// is re-weighted accordingly.
///
class RootPlotData {
private:
  TFile RootFile;

public:
  RootPlotData(
      const Physics::ParticleStateTransitionKinematicsInfo &KinematicsInfo,
      const std::string &filename, const std::string &option = "RECREATE");

  void createDirectory(std::string Name);

  void writeData(const Data::DataSet &DataSample,
                 std::string TreeName = "data");
  void writeIntensityWeightedPhspSample(
      const Data::DataSet &PhspSample, ComPWA::Intensity &Intensity,
      std::string TreeName = "intensity_weighted_phspdata",
      std::map<std::string, std::shared_ptr<ComPWA::Intensity>>
          IntensityComponents = {});
  void writeHitMissSample(const Data::DataSet &HitMissSample,
                          std::string TreeName = "hitmiss_data");
};

} // namespace Plotting
} // namespace Tools
} // namespace ComPWA

#endif
