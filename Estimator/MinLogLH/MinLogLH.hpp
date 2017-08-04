//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//	   Peter Weidenkaff - Weights and background fractions
//-------------------------------------------------------------------------------

//! Negative Log Likelihood-Estimator.
/*! \class MinLogLH
 * @file MinLogLH.hpp
 * This class calculates a simple -log(LH) of a intensity and a dataset.
 * Data and Model are provided in the constructor using the Amplitude and Data
 * interfaces. The class itself fulfills the Estimator interface.
 */

#ifndef _MINLOGLHBKG_HPP
#define _MINLOGLHBKG_HPP

#include <vector>
#include <memory>
#include <string>

// PWA-Header
#include "Estimator/Estimator.hpp"
#include "Core/AmpIntensity.hpp"
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"

namespace ComPWA {
namespace Estimator {

class MinLogLH : public ComPWA::Estimator::Estimator {

public:
  //! Destructor
  virtual ~MinLogLH(){};

  //! Implementation of ControlParameter::controlParameter
  virtual double controlParameter(ParameterList &minPar);

  /** Create instance of MinLogLH.
   * A binned efficiency correction is used. We expect that phspSample_ has
   * efficiency values for
   * each event.
   *
   * @param amp amplitude
   * @param data data sample
   * @param phspSample phsp sample for normalization. Efficiency values for each
   * point needs
   *  to be set beforehand.
   * @param startEvent use @param data from that position on
   * @param nEvents number of events to process
   * @return std::shared_ptr<Data> of existing instance or newly created
   * instance
   */
  static std::shared_ptr<ComPWA::ControlParameter>
  CreateInstance(std::shared_ptr<Kinematics> kin,
                 std::shared_ptr<AmpIntensity> intens,
                 std::shared_ptr<DataReader::Data> data,
                 std::shared_ptr<DataReader::Data> phspSample,
                 unsigned int startEvent = 0, unsigned int nEvents = 0);

  /** Create instance of MinLogLH.
   * An unbinned efficiency correction is applied using accSample.
   *
   * @param amp amplitude
   * @param data data sample
   * @param phspSample phsp sample for normalization
   * @param accSample sample of efficiency applied phsp events for unbinned
   * efficiency correction
   * @param startEvent use @param data from that position on
   * @param nEvents number of events to process
   * @return std::shared_ptr<Data> of existing instance or newly created
   * instance
   */
  static std::shared_ptr<ComPWA::ControlParameter>
  CreateInstance(std::shared_ptr<Kinematics> kin,
                 std::shared_ptr<AmpIntensity> intens,
                 std::shared_ptr<DataReader::Data> data,
                 std::shared_ptr<DataReader::Data> phspSample,
                 std::shared_ptr<DataReader::Data> accSample,
                 unsigned int startEvent = 0, unsigned int nEvents = 0);

  //! Check if tree for LH calculation is available
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

  //! Get intensity
  virtual std::shared_ptr<ComPWA::AmpIntensity> GetIntensity() {
    return _intens;
  }

  //! Get number of events in data set
  virtual int GetNEvents() { return nEvts_; }

protected:
  //! Default Constructor
  MinLogLH(){};

  //! Constructor for a single amplitude
  MinLogLH(std::shared_ptr<Kinematics> kin, std::shared_ptr<AmpIntensity> amp,
           std::shared_ptr<DataReader::Data> data,
           std::shared_ptr<DataReader::Data> phspSample,
           std::shared_ptr<DataReader::Data> accSample, unsigned int startEvent,
           unsigned int nEvents);

  //! Uses ampTree and creates a tree that calculates the full LH
  virtual void IniLHtree();

  //! Sum up all weights in data set
  void CalcSumOfWeights();

private:
  //! Initialize
  void Init();

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

#endif /* _MINLOGLHBKG_HPP */
