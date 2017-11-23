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
#include "Core/Parameter.hpp"
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

  void SetData(std::shared_ptr<ComPWA::DataReader::Data> sample) {
    s_data = sample->GetDataPoints(kin_);
  }

  void SetPhspSample(std::shared_ptr<ComPWA::DataReader::Data> sample) {
    s_phsp = sample->GetDataPoints(kin_);
  }

  void SetData(std::vector<dataPoint> &points) { s_data = points; }

  void SetPhspSample(std::vector<dataPoint> &points) { s_phsp = points; }

  void SetFitAmp(std::shared_ptr<ComPWA::AmpIntensity> intens);

  /// Add sub component of the Intensity. For each event in the phase space
  /// sample each component is evaluated and its value is added to the TTree.
  void AddComponent(std::string name, std::string title = "") {
    std::string t = title;
    if (t == "")
      t = name;
    
    if (!_plotComponents.size())
      throw std::runtime_error("PlotData::DrawComponent() | AmpIntensity not "
                               "set! Set the full model first using "
                               "SetFitAmp()!");
    
    std::shared_ptr<ComPWA::AmpIntensity> comp;
    try {
        comp = _plotComponents.at(0)->GetComponent(name);
    } catch (std::exception &ex) {
      LOG(error) << "DalitzPlot::DrawComponent() | Component " << name
                 << " not found in AmpIntensity "
                 << _plotComponents.at(0)->Name() << ".";
      return;
    }
    _plotComponents.push_back(comp);
    _componentNames.push_back(t);
  }

  /// Create the TTree's, fill and write them to \p fileName.
  /// \p treePrefix is added in front of each TTree name so that multiple
  /// TTree's can be written to the same file. Usual ROOT::TFile options
  /// can be added.
  void Write(std::string treePrefix, std::string fileName,
             std::string option = "RECREATE");

protected:
  std::shared_ptr<ComPWA::Kinematics> kin_;
  
  std::vector<std::shared_ptr<ComPWA::AmpIntensity>> _plotComponents;
  
  std::vector<std::string> _componentNames;

  std::vector<dataPoint> s_data;
  
  std::vector<dataPoint> s_phsp;
};
} // ns::Tools
} // ns::ComPWA
#endif
