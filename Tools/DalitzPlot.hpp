// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PLOTDATA_HPP_
#define PLOTDATA_HPP_

#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TGraph.h>
#include <TTree.h>

#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/FitParameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/AmpIntensity.hpp"

namespace ComPWA {
namespace Tools {

///
/// \class DalitzHisto
///  Simple class to create and fill Dalitz plots
///
class DalitzHisto {
public:
  ~DalitzHisto() {
    // Can't call delete here since there is a problem with copy/move
    // constructor
    //    if(tree) delete tree;
  }

  /// Disable copy constructor since TTree is not copyable
  DalitzHisto(const DalitzHisto &that) = delete;

  /// Default move constructor
  DalitzHisto(DalitzHisto &&other) = default; // C++11 move constructor

  DalitzHisto(std::shared_ptr<ComPWA::Kinematics> kin, std::string name,
              std::string title, unsigned int bins, Color_t color = kBlack);
  /// Switch on/off stats
  void setStats(bool b);
  /// Fill event
  void fill(std::shared_ptr<ComPWA::Kinematics> kin, ComPWA::Event &event,
            double w = 1);
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
  DalitzPlot(std::shared_ptr<ComPWA::Kinematics> kin, std::string name,
             int bins = 100);

  virtual ~DalitzPlot();

  void useEfficiencyCorrection(bool s) { _correctForEfficiency = s; }

  void setData(std::shared_ptr<ComPWA::DataReader::Data> dataSample) {
    s_data = dataSample;
  }

  void setPhspData(std::shared_ptr<ComPWA::DataReader::Data> phsp) {
    s_phsp = phsp;
  }

  void setHitMissData(std::shared_ptr<ComPWA::DataReader::Data> hitMiss) {
    s_hitMiss = hitMiss;
  }

  void setFitAmp(std::shared_ptr<ComPWA::AmpIntensity> intens,
                 std::string title = "", Color_t color = kBlack);

  void setGlobalScale(double s) { _globalScale = s; }

  void fill(std::shared_ptr<ComPWA::Kinematics> kin);

  void plot();

  void drawComponent(std::string componentName, std::string intensityName,
                     std::string title = "", Color_t color = kBlack) {
    std::string ttt = title;
    if (ttt == "")
      ttt = componentName;

    if (!_plotComponents.size())
      throw std::runtime_error("DalitzPlot::drawComponent() | AmpIntensity not "
                               "set! Set the full model first using "
                               "SetFitAmp()!");
    std::shared_ptr<ComPWA::AmpIntensity> comp;
    try {
      comp = _plotComponents.at(0)
                 ->component(intensityName)
                 ->component(componentName);
    } catch (std::exception &ex) {
      LOG(ERROR) << "DalitzPlot::drawComponent() | Component " << componentName
                 << " of " << componentName << " not found in AmpIntensity "
                 << _plotComponents.at(0)->name() << ".";
      return;
    }
    _plotComponents.push_back(comp);
    _plotHistograms.push_back(
        DalitzHisto(kin_, componentName, title, _bins, color));
    _plotHistograms.back().setStats(0);
    _plotLegend.push_back(ttt);
  }

protected:
  TString Name;

  std::shared_ptr<ComPWA::Kinematics> kin_;
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

  std::vector<std::shared_ptr<ComPWA::AmpIntensity>> _plotComponents;
  std::vector<DalitzHisto> _plotHistograms;
  std::vector<std::string> _plotLegend;

  DalitzHisto dataDiagrams;
  DalitzHisto phspDiagrams;
  DalitzHisto fitHitMissDiagrams;

  std::shared_ptr<ComPWA::DataReader::Data> s_data;
  std::shared_ptr<ComPWA::DataReader::Data> s_phsp;
  std::shared_ptr<ComPWA::DataReader::Data> s_hitMiss;
};

} // ns::Tools
} // ns::ComPWA

#endif
