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

  /// Constructor for MinLohLH.
  /// An unbinned efficiency correction can be applied using \p accSample.
  ///
  /// \param amp amplitude
  /// \param data data sample
  /// \param phspSample phsp sample for normalization
  /// \param accSample sample of efficiency applied phsp events for unbinned
  /// efficiency correction
  /// \param startEvent use @param data from that position on
  /// \param nEvents number of events to process
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
  virtual void IniLHtree();

  void CalcSumOfWeights();

private:
  void Reset();

  std::shared_ptr<Kinematics> kin_;

  /// Amplitude model
  std::shared_ptr<AmpIntensity> _intens;

  std::shared_ptr<FunctionTree> _tree;

  unsigned int nEvts_;
  unsigned int nPhsp_;

  /// Process data sample from position #nStartEvt_ on
  unsigned int nStartEvt_;
  
  /// Number of events to process in _dataSample sample
  unsigned int nUseEvt_;

  // ================== Samples =======================
  
  /// Data sample
  std::shared_ptr<DataReader::Data> _dataSample;
  
  /// _dataSample stored as ParameterList
  ParameterList _dataSampleList;
  
  /// Sum of weights in _dataSample
  double _sumOfWeights;

  /// phsp sample for normalization
  std::shared_ptr<DataReader::Data> _phspSample;
  
  /// _phspSample stored as ParameterList
  ParameterList _phspSampleList;

  /// Phsp sample with applied efficency
  std::shared_ptr<DataReader::Data> _phspAccSample;
  
  /// _phspAccSample stored as ParameterList
  ParameterList _phspAccSampleList;
  
  /// Total efficiency of phsp with applied efficency
  double _phspAccSampleEff;
};

} /* namespace Estimator */
} /* namespace ComPWA */

#endif /* _MINLOGLH_HPP */
