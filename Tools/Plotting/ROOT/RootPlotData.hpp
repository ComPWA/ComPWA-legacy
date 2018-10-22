// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_TOOLS_PLOTTING_ROOTPLOTDATA_HPP_
#define COMPWA_TOOLS_PLOTTING_ROOTPLOTDATA_HPP_

#include "Core/AmpIntensity.hpp"
#include "Core/FitParameter.hpp"
#include "Core/Kinematics.hpp"
#include "Data/Data.hpp"

namespace ComPWA {
namespace Tools {
namespace Plotting {

///
/// \class RootPlotData
/// Write data and fit result to TTree's. Simple class to output the result of a
/// fit. TTree's are generated with a branch per entry in DataPoint. The
/// Intensity is specified beforehand and also components of the intensity can
/// be defined. For each event in the phase space sample each component is
/// evaluated and its weight is added to the TTree.
///
class RootPlotData {
  std::shared_ptr<ComPWA::Kinematics> Kinematics;

  std::shared_ptr<ComPWA::AmpIntensity> Intensity;
  std::map<std::string, std::shared_ptr<ComPWA::AmpIntensity>>
      AmplitudeComponents;

  std::shared_ptr<ComPWA::Data::Data> Data;
  std::shared_ptr<ComPWA::Data::Data> WeightedPhspMC;
  std::shared_ptr<ComPWA::Data::Data> HitAndMissMC;

  bool CorrectForEfficiency;

public:
  RootPlotData(std::shared_ptr<ComPWA::Kinematics> kin,
               std::shared_ptr<ComPWA::AmpIntensity> intens);

  virtual ~RootPlotData();

  void useEfficiencyCorrection(bool s) { CorrectForEfficiency = s; }
  void setData(std::shared_ptr<ComPWA::Data::Data> dataSample) {
    Data = dataSample;
  }
  void setPhspMC(std::shared_ptr<ComPWA::Data::Data> phsp) {
    WeightedPhspMC = phsp;
  }
  void setHitMissMC(std::shared_ptr<ComPWA::Data::Data> hitMiss) {
    HitAndMissMC = hitMiss;
  }
  void addComponent(std::string componentName, std::string intensityName,
                    std::string title = "");

  /// Create the TTree's, fill and write them to \p fileName.
  /// \p treePrefix is added in front of each TTree name so that multiple
  /// TTree's can be written to the same file. Usual ROOT::TFile options
  /// can be added.
  void write(std::string treePrefix, std::string fileName,
             std::string option = "RECREATE");
};

} // namespace Plotting
} // namespace Tools
} // namespace ComPWA

#endif
