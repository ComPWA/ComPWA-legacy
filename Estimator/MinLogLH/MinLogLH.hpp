// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//! Negative Log Likelihood-Estimator.
/*! \class MinLogLH
 * @file MinLogLH.hpp
 * This class calculates a simple -log(LH) of a intensity and a dataset.
 * Data and Model are provided in the constructor using the Amplitude and Data
 * interfaces. The class itself fulfills the Estimator interface.
 */

#ifndef _MINLOGLH_HPP
#define _MINLOGLH_HPP

#include <vector>
#include <memory>
#include <string>

// PWA-Header
#include "Core/Estimator.hpp"
#include "Core/AmpIntensity.hpp"
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"

namespace ComPWA {
namespace Estimator {

class MinLogLH : public ComPWA::IEstimator {

public:
  MinLogLH(){};

  /// Create instance of MinLogLH.
  /// An unbinned efficiency correction is applied using accSample.
  /// 
  /// \param amp amplitude
  /// \param data data sample
  /// \param phspSample phsp sample for normalization
  /// \param accSample sample of efficiency applied phsp events for unbinned
  /// efficiency correction
  /// \param startEvent use @param data from that position on
  /// \param nEvents number of events to process
  /// \return std::shared_ptr<Data> of existing instance or newly created
  /// instance
  MinLogLH(std::shared_ptr<Kinematics> kin, std::shared_ptr<AmpIntensity> amp,
           std::shared_ptr<DataReader::Data> data,
           std::shared_ptr<DataReader::Data> phspSample,
           std::shared_ptr<DataReader::Data> accSample, unsigned int startEvent,
           unsigned int nEvents);
  
  virtual ~MinLogLH(){};

  virtual double controlParameter(ParameterList &minPar);

  virtual bool HasTree() { return (_tree) ? 1 : 0; }

  virtual void UseFunctionTree(bool onoff) {
    if (onoff && _tree)
      return;     // Tree already exists
    if (!onoff) { // disable tree
      _tree = std::shared_ptr<FunctionTree>();
      return;
    }
    try {
      IniLHtree();
    } catch (std::exception &ex) {
      throw std::runtime_error(
          "MinLogLH::UseFunctionTree()| FunctionTree can not be "
          "constructed! Error: " +
          std::string(ex.what()));
    }
    return;
  }

  //! Get FunctionTree for LH calculation. Check first if its available using
  //! MinLogLH::hasTree().
  virtual std::shared_ptr<FunctionTree> GetTree() {
    if (!_tree) {
      throw std::runtime_error("MinLogLH::GetTree()| FunctionTree does not "
                               "exists. Enable it first using "
                               "UseFunctionTree(true)!");
    }
    return _tree;
  }

  virtual std::shared_ptr<ComPWA::AmpIntensity> GetIntensity() {
    return _intens;
  }

  virtual int GetNEvents() { return nEvts_; }

protected:

  //! Uses ampTree and creates a tree that calculates the full LH
  virtual void IniLHtree();

  //! Sum up all weights in data set
  void CalcSumOfWeights();

private:
  //! Reset instance
  void Reset();

  std::shared_ptr<Kinematics> kin_;

  //! Intensity
  std::shared_ptr<AmpIntensity> _intens;

  //! FunctionTree for Likelihood calculation
  std::shared_ptr<FunctionTree> _tree;

  //! Number of events in data sample
  unsigned int nEvts_;
  //! Number of event in phsp sample
  unsigned int nPhsp_;

  //! Process data sample from position #nStartEvt_ on
  unsigned int nStartEvt_;
  //! Number of events to process in _dataSample sample
  unsigned int nUseEvt_;

  // Samples
  std::shared_ptr<DataReader::Data> _dataSample; //! Data sample
  ParameterList _dataSampleList; //! _dataSample stored as ParameterList
  double _sumOfWeights;          //! Sum of weights in _dataSample

  std::shared_ptr<DataReader::Data>
      _phspSample;               //! phsp sample for normalization
  ParameterList _phspSampleList; //! _phspSample stored as ParameterList

  std::shared_ptr<DataReader::Data>
      _phspAccSample;               //! Phsp sample with applied efficency
  ParameterList _phspAccSampleList; //! _phspAccSample stored as ParameterList
  double _phspAccSampleEff; //! Total efficiency of phsp with applied efficency
};

} /* namespace Estimator */
} /* namespace ComPWA */

#endif /* _MINLOGLH_HPP */
