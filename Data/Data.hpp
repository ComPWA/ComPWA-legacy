// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains Data interface class
///

#ifndef COMPWA_DATA_DATA_HPP_
#define COMPWA_DATA_DATA_HPP_

#include <vector>
#include <string>
#include <memory>

#include "DataCorrection.hpp"
#include "Core/Event.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Kinematics.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Generator.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Resolution.hpp"

namespace ComPWA {
namespace Data {

///
/// \class Data
/// This class provides the interface to experimental data. It does not provide
/// any functionality for data read in and write out. Since this is in many cases
/// not needed this class is not pure virtual. This functionality is
/// implemented in derived classes (RootReader and AsciiReader).
/// This class is also supposed to handle binned data but the implemenation of
/// this functionality currently a bit deprecated.
///
class Data {
public:
  Data(bool binning = 0, unsigned int maxBins = 0, double maxW = 0.0);
  
  virtual ~Data() {};

  /// Append \p otherSample to the current one.
  virtual void append(Data &otherSample);

  /// Add \p event to sample.
  virtual void add(const Event &event);

  /// Add \p event to sample.
  virtual void add(const std::vector<Event> &events);

  virtual const std::size_t numEvents() const { return Events.size(); }

  /// Get event by its position \p id.
  virtual Event &event(const int id) { return Events.at(id); }

  /// Get a reference to the internal event storage.
  virtual std::vector<Event> &events() { return Events; }

  /// Get 'horizontal' list of dataPoints. For each variable
  /// (e.g. m23sq, m13sq ...) a MultiDouble is added to ParameterList.
  /// This ParameterList is used to build the FunctionTree.
  virtual const ParameterList &dataList(std::shared_ptr<Kinematics> kin);

  /// Get 'vertical' list of dataPoints.
  std::vector<DataPoint> dataPoints(std::shared_ptr<Kinematics> kin) const;

  /// Apply a correction to the sample.
  /// E.g. correction for efficiency differences between data and Monte-Carlo.
  virtual void applyCorrection(DataCorrection &corr);

  /// Remove all events outside phase space boundaries from this sample.
  virtual void reduceToPhsp(std::shared_ptr<Kinematics> kin);

  /// reduce data set to \p newSize.
  virtual void reduce(unsigned int newSize);

  /// Return a random subset of this sample.
  virtual std::shared_ptr<Data> rndSubSet(std::shared_ptr<Kinematics> kin,
                                          unsigned int size,
                                          std::shared_ptr<Generator> gen);

  void setResolution(std::shared_ptr<Resolution> res);

  /// Set efficiency for each events.
  /// Since the efficiency is usually calculated in terms of phase-space variables
  /// we have to pass a Kinematics object.
  virtual void setEfficiency(std::shared_ptr<Kinematics> kin,
                             std::shared_ptr<Efficiency> eff);

  virtual void resetEfficiency(double eff = 1.);

  /// Get maximum weight of all events.
  virtual double maximumWeight() const;

  virtual void resetWeights(double weight = 1.);

  /// Check if events have weights different from 1.0.
  virtual bool hasWeights();

  /// Delete all events.
  /// An empty Data object remains.
  virtual void clear();

  /// \deprecated
  virtual const std::size_t numBins() const { return fBins.size(); }

  /// \deprecated
  virtual const int bin(const int, double &, double &);

protected:
  /// DataPoints are stored as 'horizontal' structure
  ParameterList DataList;

  /// Vector of events
  std::vector<Event> Events;

  /// Maximum weight of events
  double MaximumWeight;

  /// Binning
  bool fBinned;
  unsigned int fmaxBins;
  std::map<int, std::pair<double, double>> fBins;
};

//===== DATASET OPERATIONS

/// Select random sub sample of data sets.
/// A hit&miss procedure is applied to the first sample to select a random set of
/// events. In case an
/// event of the first sample is added to the output sample, an event from the
/// second sample is also
/// added to the corresponding output sample. This function can be used for two
/// sample that are in sync.
/// E.g. a sample with reconstructed values and a sample with the corresponding
/// true values
/// We expect that @param out1 and @param out2 are pointers to empty samples.
/// 
/// \param size size of sub sample
/// \param gen generator
/// \param in1 input sample 1
/// \param in2 input sample 2
/// \param out1 output sub sample from sample 1
/// \param out2 output sub sample from sample 2
inline void rndReduceSet(std::shared_ptr<Kinematics> kin, unsigned int size,
                         std::shared_ptr<Generator> gen, Data *in1, Data *out1,
                         Data *in2 = NULL, Data *out2 = NULL);
                         
} // ns::Data
} // ns::ComPWA

#endif
