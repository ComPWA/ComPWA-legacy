// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <numeric>

#include "Estimator/MinLogLH/MinLogLH.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Kinematics.hpp"
#include "Core/FitResult.hpp"

using namespace ComPWA;
using namespace ComPWA::Estimator;

MinLogLH::MinLogLH(std::shared_ptr<Kinematics> kin,
                   std::shared_ptr<AmpIntensity> intens,
                   std::shared_ptr<DataReader::Data> data,
                   std::shared_ptr<DataReader::Data> phspSample,
                   std::shared_ptr<DataReader::Data> accSample,
                   unsigned int firstEvent, unsigned int nEvents)
    : _kin(kin), _intens(intens), _firstEvent(firstEvent), _nEvents(nEvents),
      _dataSample(data), _phspSample(phspSample), _phspAccSample(accSample),
      _phspAccSampleEff(1.0) {

  int size = _dataSample->GetNEvents();

  // use the full sample of both are zero
  if (!_nEvents && !_firstEvent) {
    _nEvents = size;
    _firstEvent = 0;
  }

  // Are more events requested than available in data sample?
  if (_firstEvent + _nEvents > size)
    _nEvents = size - _firstEvent;

  // Get data as ParameterList
  _dataSampleList = _dataSample->GetListOfData(_kin);
  _phspSampleList = _phspSample->GetListOfData(_kin);
  if (_phspAccSample)
    _phspAccSampleList = _phspAccSample->GetListOfData(_kin);
  else
    _phspAccSampleList = _phspSample->GetListOfData(_kin);

  // Calculation sum of weights of data sample
  _sumOfWeights = 0;
  for (unsigned int evt = _firstEvent; evt < _nEvents + _firstEvent; evt++) {
    Event ev(_dataSample->GetEvent(evt));
    _sumOfWeights += ev.GetWeight();
  }

  LOG(info) << "MinLogLH::Init() |  Size of data sample = " << _nEvents
            << " ( Sum of weights = " << _sumOfWeights << " ).";

  _nCalls = 0; // member of ControlParameter

  return;
}

double MinLogLH::ControlParameter(ParameterList &minPar) {
  double lh = 0;
  if (!_tree) {
    // Calculate \Sum_{ev} log()
    double sumLog = 0;
    // loop over data sample
    for (unsigned int evt = _firstEvent; evt < _nEvents + _firstEvent; evt++) {
      dataPoint point;
      _kin->EventToDataPoint(_dataSample->GetEvent(evt), point);
      double val = _intens->Intensity(point);
      sumLog += std::log(val) * point.GetWeight();
    }
    lh = (-1) * ((double)_nEvents) / _sumOfWeights * sumLog;
  } else {
    _tree->Recalculate();
    std::shared_ptr<DoubleParameter> logLH =
        std::dynamic_pointer_cast<DoubleParameter>(_tree->Head()->Parameter());
    lh = logLH->GetValue();
  }
  _nCalls++;
  return lh; // return -logLH
}

void MinLogLH::UseFunctionTree(bool onoff) {
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

std::shared_ptr<FunctionTree> MinLogLH::GetTree() {
  if (!_tree) {
    throw std::runtime_error("MinLogLH::GetTree()| FunctionTree does not "
                             "exists. Enable it first using "
                             "UseFunctionTree(true)!");
  }
  return _tree;
}

void MinLogLH::IniLHtree() {
  LOG(debug) << "MinLogLH::IniLHtree() | Constructing FunctionTree!";

  // Ensure that a FunctionTree is provided
  if (!_intens->HasTree())
    throw std::runtime_error("MinLogLH::IniLHtree() |  AmpIntensity does not "
                             "provide a FunctionTree!");

  _tree = std::shared_ptr<FunctionTree>(new FunctionTree());
  int sampleSize = _dataSampleList.GetMultiDouble(0)->GetNValues();

  //-log L = (-1)*N/(\sum_{ev} w_{ev}) \sum_{ev} ...
  _tree->CreateHead("LH",
                    std::shared_ptr<Strategy>(new MultAll(ParType::DOUBLE)));
  _tree->CreateLeaf("minusOne", -1, "LH");
  _tree->CreateLeaf("nEvents", sampleSize, "LH");
  _tree->CreateNode("invSumWeights",
                    std::shared_ptr<Strategy>(new Inverse(ParType::DOUBLE)),
                    "LH");
  _tree->CreateNode("sumEvents",
                    std::shared_ptr<Strategy>(new AddAll(ParType::DOUBLE)),
                    "LH");
  _tree->CreateLeaf("SumOfWeights", _sumOfWeights, "invSumWeights");
  _tree->CreateNode("weightLog",
                    std::shared_ptr<Strategy>(new MultAll(ParType::MDOUBLE)),
                    "sumEvents", sampleSize,
                    false); // w_{ev} * log( I_{ev} )
  _tree->CreateLeaf("Weight",
                    _phspAccSampleList.GetMultiDouble("Weight"), "weightLog");
  _tree->CreateNode("Log",
                    std::shared_ptr<Strategy>(new LogOf(ParType::MDOUBLE)),
                    "weightLog", sampleSize, false);
  _tree->InsertTree(_intens->GetTree(_kin, _dataSampleList, _phspAccSampleList,
                                     _phspSampleList, _kin->GetNVars()),
                    "Log");

  _tree->Recalculate();
  if (!_tree->SanityCheck()) {
    throw std::runtime_error("MinLogLH::IniLHtree() | Tree has structural "
                             "problems. Sanity check not passed!");
  }
  LOG(debug) << "MinLogLH::IniLHtree() | "
                "Construction of LH tree finished!";
  return;
}
