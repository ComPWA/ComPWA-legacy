// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// RootPlot class
///

#ifndef RootPlot_hpp
#define RootPlot_hpp

#include <stdio.h>
#include <TTree.h>

#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/FitParameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/AmpIntensity.hpp"

namespace ComPWA {
namespace Tools {

///
/// \class RootPlot
/// Write data and fit result to TTree's. Simple class to ouput the result of a
/// fit. TTree's are generated with a branch per entry in dataPoint. The
/// Intensity is specified beforehand and also components of the intensity can
/// be defined. For each event in the phase space sample each component is
/// evaluated and its weight is added to the TTree.
///
class RootPlot {
public:
  RootPlot(std::shared_ptr<ComPWA::Kinematics> kin);

  virtual ~RootPlot() {}

  void setDataSample(std::shared_ptr<ComPWA::DataReader::Data> sample) {
    DataSample = sample->dataPoints(Kin);
  }

  void setPhspSample(std::shared_ptr<ComPWA::DataReader::Data> sample) {
    PhspSample = sample->dataPoints(Kin);
  }

  void setDataSample(std::vector<DataPoint> &points) { DataSample = points; }

  void setPhspSample(std::vector<DataPoint> &points) { PhspSample = points; }

  void setIntensity(std::shared_ptr<ComPWA::AmpIntensity> intens);
  
  /// Add sub component of the Intensity. For each event in the phase space
  /// sample each component is evaluated and its value is added to the TTree.
  void addComponent(std::string componentName, std::string intensityName,
                    std::string title = "") {
    std::string ttt = title;
    if (ttt == "")
      ttt = componentName;
    
    if (!PlotComponents.size())
      throw std::runtime_error("RootPlot::addComponent() | AmpIntensity not "
                               "set! Set the full model first using "
                               "SetFitAmp()!");
    
    std::shared_ptr<ComPWA::AmpIntensity> comp;
    try {
      comp = PlotComponents.at(0)
                 ->component(intensityName)
                 ->component(componentName);
    } catch (std::exception &ex) {
      LOG(ERROR) << "RootPlot::addComponent() | Component " << componentName
                 << " of " << componentName <<" not found in AmpIntensity "
                 << PlotComponents.at(0)->name() << ".";
      return;
    }
    PlotComponents.push_back(comp);
    ComponentNames.push_back(ttt);
  }

  /// Create the TTree's, fill and write them to \p fileName.
  /// \p treePrefix is added in front of each TTree name so that multiple
  /// TTree's can be written to the same file. Usual ROOT::TFile options
  /// can be added.
  void write(std::string treePrefix, std::string fileName,
             std::string option = "RECREATE");

protected:
  std::shared_ptr<ComPWA::Kinematics> Kin;
  
  std::vector<std::shared_ptr<ComPWA::AmpIntensity>> PlotComponents;
  
  std::vector<std::string> ComponentNames;

  std::vector<DataPoint> DataSample;
  
  std::vector<DataPoint> PhspSample;
};
} // ns::Tools
} // ns::ComPWA
#endif
