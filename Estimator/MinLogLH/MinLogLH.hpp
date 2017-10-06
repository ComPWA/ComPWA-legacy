// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
///

#ifndef _MINLOGLH_HPP
#define _MINLOGLH_HPP

#include <vector>
#include <memory>
#include <string>

#include "Core/Estimator.hpp"
#include "Core/AmpIntensity.hpp"
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"

namespace ComPWA {
namespace Estimator {

///
/// \class MinLogLH
/// Negative Log Likelihood-Estimator. This class calculates a simple
/// -log(LH) using an AmpIntensity and datasets for data and phase space.
/// The class fulfills the Estimator interface.
///
/// \par Log likelihood
/// The negative log LH is given by:
/// \f[
///    -log \mathcal{L} = - \frac{N}{\sum_{ev} w_{ev}} \sum_{ev} w_{ev}
///    \log T(ev).
/// \f]
/// With the the value of AmpIntensity \f$T(ev)\f$ for a given event \f$ev\f$.
/// The sum over all weights is necessary to normalize the weights to one.
/// Otherwise the error estimate would be incorrect. AmpIntensity is normalized
/// over the phase space region and therefore the likelihood function is also
/// normalized.
///
/// \par Efficiency correction
/// Efficiency correction is performed within AmpIntensity, more details are
/// given in AmpIntensity. Either each event in phspSample has its
/// efficiency stored inside or another phase space sample accSample must be
/// provided to which the efficiency is applied. That means that accSample
/// has passed reconstruction and selection.
///
class MinLogLH : public ComPWA::IEstimator {

public:
  /// Constructor for MinLogLH.
  MinLogLH(std::shared_ptr<ComPWA::Kinematics> kin,
           std::shared_ptr<ComPWA::AmpIntensity> amp,
           std::shared_ptr<ComPWA::DataReader::Data> data,
           std::shared_ptr<ComPWA::DataReader::Data> phspSample,
           std::shared_ptr<ComPWA::DataReader::Data> accSample,
           unsigned int firstEvent, unsigned int nEvents);

  virtual ~MinLogLH() {};
  
  /// Value of minimum log likelhood function.
  virtual double ControlParameter(ComPWA::ParameterList &par);

  /// Trigger the use of a FunctionTree.
  /// If no tree is provided by the AmpIntensity implementation an exception
  /// is thrown.
  virtual void UseFunctionTree(bool onoff = true);

  /// Get the FunctionTree.
  /// If no FunctionTree is available an std::runtime_error exception is thrown.
  virtual std::shared_ptr<ComPWA::FunctionTree> GetTree();

  /// Number of likelihood evaluations
  virtual int NSteps() { return _nCalls; }

protected:
  /// Initialize FunctionTree.
  /// If the AmpIntensity model does not provide a FunctionTree or if a tree
  /// with error is constructed, an exception is thrown.
  virtual void IniLHtree();

  std::shared_ptr<ComPWA::Kinematics> _kin;

  /// AmpIntensity model
  std::shared_ptr<ComPWA::AmpIntensity> _intens;

  std::shared_ptr<ComPWA::FunctionTree> _tree;

  /// Process data sample from position _firstEvent on
  unsigned int _firstEvent;

  /// Number of events to process in _dataSample sample
  unsigned int _nEvents;

  /// Number of likelihood evaluations
  int _nCalls;

  // ================== Samples =======================

  /// Data sample
  std::shared_ptr<ComPWA::DataReader::Data> _dataSample;

  /// _dataSample stored 'horizontally' as ParameterList
  ParameterList _dataSampleList;

  /// Sum of weights in _dataSample
  double _sumOfWeights;

  /// phsp sample for normalization
  std::shared_ptr<ComPWA::DataReader::Data> _phspSample;

  /// _phspSample stored 'horizontally' as ParameterList
  ParameterList _phspSampleList;

  /// Phsp sample with applied efficency
  std::shared_ptr<ComPWA::DataReader::Data> _phspAccSample;

  /// _phspAccSample 'horizontally' stored as ParameterList
  ParameterList _phspAccSampleList;

  /// Total efficiency of phsp with applied efficency
  double _phspAccSampleEff;
};

} // namespace Estimator
} // namespace ComPWA

#endif // _MINLOGLH_HPP
