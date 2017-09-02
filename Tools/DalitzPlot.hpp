// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PLOTDATA_HPP_
#define PLOTDATA_HPP_

#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TGraph.h>

#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/AmpIntensity.hpp"
#include "DataReader/RootReader/RootReader.hpp"

///
/// \class DalitzHisto
///  Simple class to create and fill Dalitz plots
///
class DalitzHisto {
public:
  ~DalitzHisto() {
    // Can't call delete here since there is a problem with copy/move
    // constructor
    //		if(tree) delete tree;
  }

  //! Disable copy constructor since TTree is not copyable
  DalitzHisto(const DalitzHisto &that) = delete;

  //! Default move constructor
  DalitzHisto(DalitzHisto &&other) = default; // C++11 move constructor

  DalitzHisto(std::shared_ptr<ComPWA::Kinematics> kin, std::string name, std::string title, unsigned int bins,
              Color_t color = kBlack);
  //! Switch on/off stats
  void SetStats(bool b);
  //! Fill event
  void Fill(std::shared_ptr<ComPWA::Kinematics> kin, ComPWA::Event &event, double w = 1);
  //! Scale all distributions
  void Scale(double w);
  //! Get 1D histogram
  TH1D *getHistogram(unsigned int num);
  //! Get 2D histogram
  TH2D *getHistogram2D(unsigned int num);
  //! set line color
  void setColor(Color_t color);
  //! Write to TFile
  void Write();
  //! GetIntegral
  double GetIntegral() { return _integral; }

private:
  std::vector<TH1D> _arr;
  std::vector<TH2D> _arr2D;
  std::string _name, _title;
  unsigned int _nBins;

  std::unique_ptr<TTree> _tree;
  // tree branches
  std::vector<double> t_point;
  double t_eff, t_weight;
  double _integral;
  Color_t _color;
};

class DalitzPlot {
public:
  DalitzPlot(std::shared_ptr<ComPWA::Kinematics> kin, std::string name, int bins = 100);

  virtual ~DalitzPlot();

  void UseEfficiencyCorrection(bool s) { _correctForEfficiency = s; }

  void SetData(std::shared_ptr<ComPWA::DataReader::Data> dataSample) {
    s_data = dataSample;
  }

  void SetPhspData(std::shared_ptr<ComPWA::DataReader::Data> phsp) { s_phsp = phsp; }

  void SetHitMissData(std::shared_ptr<ComPWA::DataReader::Data> hitMiss) {
    s_hitMiss = hitMiss;
  }

  void SetFitAmp(std::shared_ptr<ComPWA::AmpIntensity> intens,
                 std::string title = "", Color_t color = kBlack);

  void SetGlobalScale(double s) { _globalScale = s; }

  void Fill(std::shared_ptr<ComPWA::Kinematics> kin);

  void Plot();

  void DrawComponent(std::string name, std::string title = "",
                     Color_t color = kBlack) {
    if (!_plotComponents.size())
      throw std::runtime_error("PlotData::DrawComponent() | AmpIntensity not "
                               "set! Set the full model first using "
                               "SetFitAmp()!");
    std::shared_ptr<ComPWA::AmpIntensity> comp;
    try{
      comp = _plotComponents.at(0)->GetComponent(name);
    } catch (std::exception& ex) {
      LOG(error) << "DalitzPlot::DrawComponent() | Component " << name
                 << " not found in AmpIntensity "
                 << _plotComponents.at(0)->Name() << ".";
      return;
    }
    _plotComponents.push_back(comp);
    _plotHistograms.push_back(DalitzHisto(kin_, name, title, _bins, color));
    _plotHistograms.back().SetStats(0);
    _plotLegend.push_back(title);
  }

protected:
  TString _name;

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

#endif
