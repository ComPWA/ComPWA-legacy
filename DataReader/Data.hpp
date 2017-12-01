// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Contains Data interface class
///

#ifndef DATA_HPP_
#define DATA_HPP_

#include <vector>
#include <string>
#include <memory>

#include "Core/Event.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Kinematics.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Generator.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Resolution.hpp"
#include "DataReader/DataCorrection.hpp"

namespace ComPWA {
namespace DataReader {

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

  virtual Data *Clone() const { return new Data(*this); };

  virtual Data *EmptyClone() const { return new Data(); };

  /// Append \p otherSample to the current one.
  virtual void Add(Data &otherSample);

  /// Add \p event to sample.
  virtual void PushEvent(const Event &event);

  virtual const std::size_t GetNEvents() const { return fEvents.size(); }

  /// Get event by its position \p id.
  virtual Event &GetEvent(const int id) { return fEvents.at(id); }

  /// Get a reference to the internal event storage.
  virtual std::vector<Event> &GetEvents() { return fEvents; }

  /// Get 'horizontal' list of dataPoints. For each variable
  /// (e.g. m23sq, m13sq ...) a MultiDouble is added to ParameterList.
  /// This ParameterList is used to build the FunctionTree.
  virtual const ParameterList &GetListOfData(std::shared_ptr<Kinematics> kin);

  /// Get 'vertical' list of dataPoints.
  std::vector<DataPoint> GetDataPoints(std::shared_ptr<Kinematics> kin) const;

  /// Apply a correction to the sample.
  /// E.g. correction for efficiency differences between data and Monte-Carlo.
  virtual void ApplyCorrection(DataCorrection &corr);

  /// Remove all events outside phase space boundaries from this sample.
  virtual void ReduceToPhsp(std::shared_ptr<Kinematics> kin);

  /// Reduce data set to \p newSize.
  virtual void Reduce(unsigned int newSize);

  /// Return a random subset of this sample.
  virtual std::shared_ptr<Data> RndSubSet(std::shared_ptr<Kinematics> kin,
                                          unsigned int size,
                                          std::shared_ptr<Generator> gen);

  void SetResolution(std::shared_ptr<Resolution> res);

  /// Set efficiency for each events.
  /// Since the efficiency is usually calculated in terms of phase-space variables
  /// we have to pass a Kinematics object.
  virtual void SetEfficiency(std::shared_ptr<Kinematics> kin,
                             std::shared_ptr<Efficiency> eff);

  virtual void ResetEfficiency(double eff = 1.);

  /// Get maximum weight of all events.
  virtual double GetMaxWeight() const;

  virtual void ResetWeights(double weight = 1.);

  /// Check if events have weights different from 1.0.
  virtual bool HasWeights();

  /// Delete all events.
  /// An empty Data object remains.
  virtual void Clear();

  /// Write sample to file.
  /// Method is supposed to be implemented by derived classes (RootReader,
  /// AsciiReader). Sine Data is not a pure virtual class we implement a function
  /// with an error message.
  virtual void WriteData(std::string file = "", std::string trName = "") {
    LOG(error)
        << "Data::writeData() | Base class does not provide functionality"
           "to write data to file.";
  };

  /// \deprecated
  virtual const std::size_t GetNBins() const { return fBins.size(); }

  /// \deprecated
  virtual const int GetBin(const int, double &, double &);

protected:
  /// DataPoints are stored as 'horizontal' structure
  ParameterList dataList;

  /// Vector of events
  std::vector<Event> fEvents;

  /// Maximum weight of events
  double maxWeight;

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
                         
} /* namespace DataReader */
} /* namespace ComPWA */

#endif
