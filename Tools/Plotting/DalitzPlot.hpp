// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PLOTDATA_HPP_
#define PLOTDATA_HPP_

#include <map>
#include <vector>

#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TTree.h>

#include "Core/FitParameter.hpp"
#include "Core/Function.hpp"
#include "Core/FunctionTree/FunctionTreeIntensity.hpp"
#include "Core/Particle.hpp"

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
  DalitzHisto(DalitzHisto &&other) = default;

  DalitzHisto(ComPWA::Physics::HelicityFormalism::HelicityKinematics &helkin,
              std::string name, std::string title, unsigned int bins,
              Color_t color = kBlack);
  /// Switch on/off stats
  void setStats(bool b);

  void fill(const ComPWA::Data::DataSet &sample,
            std::vector<double> Intensities = {});

  /// Scale all distributions
  void scale(double w);
  /// Get 1D histogram
  TH1D *getHistogram(std::string Name);
  /// Get 2D histogram
  TH2D *getHistogram2D(std::pair<std::string, std::string> Names);
  /// set line color
  void setColor(Color_t color);
  /// Write to TFile
  void write();
  /// GetIntegral
  double integral() { return Integral; }

private:
  std::map<std::string, TH1D> Hists1D;
  std::map<std::pair<std::string, std::string>, TH2D> Hists2D;
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
             const std::string &name, int bins = 100);

  virtual ~DalitzPlot() = default;

  void setGlobalScale(double s) { _globalScale = s; }

  void fill(const std::vector<ComPWA::Event> &data, bool normalize = false,
            const std::string &name = "", const std::string &title = "",
            Color_t color = kBlack);
  void fill(const std::vector<ComPWA::Event> &data, Intensity &intens,
            bool normalize = false, const std::string &name = "",
            const std::string &title = "", Color_t color = kBlack);
  void plot();

private:
  TString Name;

  Physics::HelicityFormalism::HelicityKinematics &HelKin;

  unsigned int _bins;

  double _globalScale;

  void CreateHist(std::string Name);

  std::vector<DalitzHisto> _plotHistograms;
};

} // namespace Plotting
} // namespace Tools
} // namespace ComPWA
#endif
