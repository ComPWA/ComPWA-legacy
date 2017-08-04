/*Author: Peter Weidenkaff
 * Date: 2013-08-07
 */
#ifndef PLOTDATA_HPP_
#define PLOTDATA_HPP_

#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TGraph.h>

// Core header files go here
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/AmpIntensity.hpp"
#include "DataReader/RootReader/RootReader.hpp"

using namespace ComPWA;

/* \class dalitzHisto
 *! Simple class to create and fill Dalitz plots
 */
class dalitzHisto {
public:
  ~dalitzHisto() {
    // Can't call delete here since there is a problem with copy/move
    // constructor
    //		if(tree) delete tree;
  }

  //! Disable copy constructor since TTree is not copyable
  dalitzHisto(const dalitzHisto &that) = delete;

  //! Default move constructor
  dalitzHisto(dalitzHisto &&other) = default; // C++11 move constructor

  dalitzHisto(std::shared_ptr<Kinematics> kin, std::string name, std::string title, unsigned int bins,
              Color_t color = kBlack);
  //! Switch on/off stats
  void SetStats(bool b);
  //! Fill event
  void Fill(std::shared_ptr<Kinematics> kin, Event &event, double w = 1);
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

class plotData {
public:
  plotData(std::shared_ptr<Kinematics> kin, std::string name, int bins = 100);

  virtual ~plotData();

  void UseEfficiencyCorrection(bool s) { _correctForEfficiency = s; }

  void SetData(std::shared_ptr<DataReader::Data> dataSample) {
    s_data = dataSample;
  }

  void SetPhspData(std::shared_ptr<DataReader::Data> phsp) { s_phsp = phsp; }

  void SetHitMissData(std::shared_ptr<DataReader::Data> hitMiss) {
    s_hitMiss = hitMiss;
  }

  void SetFitAmp(std::shared_ptr<ComPWA::AmpIntensity> intens,
                 std::string title = "", Color_t color = kBlack);

  void SetGlobalScale(double s) { _globalScale = s; }

  void Fill(std::shared_ptr<Kinematics> kin);

  void Plot();

  void DrawComponent(std::string name, std::string title = "",
                     Color_t color = kBlack) {
    if (!_plotComponents.size())
      throw std::runtime_error("PlotData::DrawComponent() | AmpIntensity not "
                               "set! Set the full model first using "
                               "SetFitAmp()!");
    std::shared_ptr<AmpIntensity> comp;
    try{
      comp = _plotComponents.at(0)->GetComponent(name);
    } catch (std::exception& ex) {
      LOG(error) << "plotData::DrawComponent() | Component " << name
                 << " not found in AmpIntensity "
                 << _plotComponents.at(0)->Name() << ".";
      return;
    }
    _plotComponents.push_back(comp);
    _plotHistograms.push_back(dalitzHisto(kin_, name, title, _bins, color));
    _plotHistograms.back().SetStats(0);
    _plotLegend.push_back(title);
  }

protected:
  TString _name;

  std::shared_ptr<Kinematics> kin_;
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

  std::vector<std::shared_ptr<AmpIntensity>> _plotComponents;
  std::vector<dalitzHisto> _plotHistograms;
  std::vector<std::string> _plotLegend;

  dalitzHisto dataDiagrams;
  dalitzHisto phspDiagrams;
  dalitzHisto fitHitMissDiagrams;

  std::shared_ptr<DataReader::Data> s_data;
  std::shared_ptr<DataReader::Data> s_phsp;
  std::shared_ptr<DataReader::Data> s_hitMiss;
};

#endif
