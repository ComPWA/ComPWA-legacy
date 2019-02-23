// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_TOOLS_PLOTTING_ROOTPLOTDATA_HPP_
#define COMPWA_TOOLS_PLOTTING_ROOTPLOTDATA_HPP_

#include "Core/FitParameter.hpp"

#include "TFile.h"

namespace ComPWA {
class Intensity;
namespace Data {
class DataSet;
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

  void writeData(const Data::DataSet &DataSample);
  void writeIntensityWeightedPhspSample(
      const Data::DataSet &PhspSample,
      std::shared_ptr<ComPWA::Intensity> Intensity,
      std::map<std::string, std::shared_ptr<const ComPWA::Intensity>>
          IntensityComponents = {});
  void writeHitMissSample(const Data::DataSet &HitMissSample);
};

} // namespace Plotting
} // namespace Tools
} // namespace ComPWA

#endif
