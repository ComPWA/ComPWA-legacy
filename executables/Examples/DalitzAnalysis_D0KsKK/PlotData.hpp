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
#include "Core/Amplitude.hpp"
#include "Core/Particle.hpp"
#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "DataReader/RootReader/RootReader.hpp"
#include "Physics/AmplitudeSum/AmpSumIntensity.hpp"

using namespace ComPWA;

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

  dalitzHisto() : _integral(0.0){};

  dalitzHisto(std::string n, std::string t, unsigned int bins);
  //! Switch on/off stats
  void SetStats(bool b);
  //! Fill event
  void Fill(Event &event, double w = 1);
    //! Fill dataPoint
    void Fill(dataPoint &point, double w = 1);
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

  static TH2Poly *getTH2PolyPull(TH2Poly *hist1, TH2Poly *hist2,
                                 TString name = "residual_");

private:
  std::vector<TH1D> arr;
  std::vector<TH2D> arr2D;
  std::string name, title;
  unsigned int nBins;

  std::unique_ptr<TTree> tree;
  // tree branches
  std::vector<double> t_point;
  double t_eff, t_weight;
  double _integral;
};

class plotData {
public:
  plotData(std::string name, int bins = 100);

  virtual ~plotData();

  void setCorrectEfficiency(bool s) { _correctForEfficiency = s; }

  void setData(std::shared_ptr<DataReader::Data> dataSample) {
    s_data = dataSample;
  }

  void setPhspData(std::shared_ptr<DataReader::Data> phsp) { s_phsp = phsp; }

  void setFitData(std::shared_ptr<DataReader::Data> fit) { s_fit = fit; }

  void setHitMissData(std::shared_ptr<DataReader::Data> hitMiss) {
    s_hitMiss = hitMiss;
  }

  void setFitAmp(std::vector<std::shared_ptr<Amplitude>> ampVec,
                 std::vector<double> fraction);

  TH2Poly *getAdBinHist(int bins = 30);

  void setGlobalScale(double s) { _globalScale = s; }

  void Fill();

  void Plot();

  void DrawComponent(unsigned int i, Color_t color) {
    if (i >= ampHistos.size())
      return;
    plotComponent.push_back(std::make_pair(i, color));
  }

protected:
  TString _name;
  bool _isFilled;
  unsigned int _bins;

  double _globalScale;
  bool _correctForEfficiency;
  void CreateHist(unsigned int id);
  void CreateHist2(unsigned int id);

  TH2Poly *dataAdaptiveBinningHist;
  TH2Poly *fitAdaptiveBinningHist;
  TH1D h_weights;
  TGraph m23m13_contour;
  TGraph m23m12_contour;
  TGraph m12m13_contour;

  std::vector<dalitzHisto> ampHistos;
  std::vector<dalitzHisto> signalComponents;
  std::vector<std::pair<unsigned int, Color_t>> plotComponent;
  dalitzHisto dataDiagrams;
  dalitzHisto phspDiagrams;
  dalitzHisto fitDiagrams;
  dalitzHisto fitHitMissDiagrams;

  std::shared_ptr<DataReader::Data> s_data;
  std::shared_ptr<DataReader::Data> s_phsp;
  std::shared_ptr<DataReader::Data> s_fit;
  std::shared_ptr<DataReader::Data> s_hitMiss;
  std::vector<std::shared_ptr<Amplitude>> _ampVec;
  std::vector<double> _fraction;
};

#endif
