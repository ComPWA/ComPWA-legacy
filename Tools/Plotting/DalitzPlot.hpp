// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PLOTDATA_HPP_
#define PLOTDATA_HPP_

#include <vector>

#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TTree.h>

#include "Core/FitParameter.hpp"
#include "Core/Function.hpp"
#include "Core/Particle.hpp"
#include "Core/FunctionTree/FunctionTreeIntensity.hpp"

namespace ComPWA {
struct DataPoint;
struct Event;
namespace Data {
struct DataSet;
}
namespace Physics {
namespace HelicityFormalism {
class HelicityKinematics;
}
} // namespace Physics
namespace Tools {
namespace Plotting {

using ComPWA::FunctionTree::FunctionTreeIntensity;
///
/// \class DalitzHisto
///  Simple class to create and fill Dalitz plots
///
class DalitzHisto {
public:
  /// Disable copy constructor since TTree is not copyable
  DalitzHisto(const DalitzHisto &that) = delete;

  /// Default move constructor
  DalitzHisto(DalitzHisto &&other) = default; // C++11 move constructor

  DalitzHisto(ComPWA::Physics::HelicityFormalism::HelicityKinematics &helkin,
              std::string name, std::string title, unsigned int bins,
              Color_t color = kBlack);
  /// Switch on/off stats
  void setStats(bool b);
  /// Fill event
  void fill(const ComPWA::Physics::HelicityFormalism::HelicityKinematics &kin,
            const ComPWA::DataPoint &point, double w = 1);
  /// Scale all distributions
  void scale(double w);
  /// Get 1D histogram
  TH1D *getHistogram(unsigned int num);
  /// Get 2D histogram
  TH2D *getHistogram2D(unsigned int num);
  /// set line color
  void setColor(Color_t color);
  /// Write to TFile
  void write();
  /// GetIntegral
  double integral() { return Integral; }

private:
  std::vector<TH1D> Arr;
  std::vector<TH2D> Arr2D;
  std::string Name, Title;
  unsigned int NumBins;

  std::unique_ptr<TTree> Tree;

  // tree branches
  std::vector<double> BranchPoint;
  double BranchEff, BranchWeight;
  double Integral;
};

class DalitzPlot {
public:
  DalitzPlot(ComPWA::Physics::HelicityFormalism::HelicityKinematics &kin,
             std::string name, int bins = 100);

  virtual ~DalitzPlot() = default;

  void useEfficiencyCorrection(bool s) { _correctForEfficiency = s; }

  void setFitAmp(FunctionTreeIntensity& intens, std::string name,
                 std::string title = "", Color_t color = kBlue);

  void setGlobalScale(double s) { _globalScale = s; }

  void fillData(const std::vector<ComPWA::Event> &data);
  void fillPhaseSpaceData(const std::vector<ComPWA::Event> &data,
                          FunctionTreeIntensity& intens,
                          std::string name, std::string title, Color_t color);
  void fillHitAndMissData(const std::vector<ComPWA::Event> &data);

  void plot();

  void drawComponent(FunctionTreeIntensity& component,
                     std::string componentName, std::string title = "",
                     Color_t color = kBlack) {
    _plotComponents.push_back(component);
    _plotHistograms.push_back(
        DalitzHisto(HelKin, componentName, title, _bins, color));
    _plotHistograms.back().setStats(0);
    _plotLegend.push_back(componentName);
  }

private:
  TString Name;

  Physics::HelicityFormalism::HelicityKinematics &HelKin;

  bool _isFilled;

  unsigned int _bins;

  double _globalScale;

  bool _correctForEfficiency;

  void CreateHist(unsigned int id);
  void CreateHist2(unsigned int id);

  TH1D h_weights;
  TGraph m23m13_contour;
  TGraph m23m12_contour;
  TGraph m12m13_contour;

  std::vector<FunctionTreeIntensity> _plotComponents;
  std::vector<DalitzHisto> _plotHistograms;
  std::vector<std::string> _plotLegend;

  DalitzHisto dataDiagrams;
  DalitzHisto phspDiagrams;
  DalitzHisto fitHitMissDiagrams;
};

} // namespace Plotting
} // namespace Tools
} // namespace ComPWA
#endif
