// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

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

/*! \class Data
 * This class provides the interface to experimental data. It does not provide
 * any funktionality for data read in and write out. This funktionality is
 * in derived classes. See @RootReader, @AsciiReader.
 */
class Data {
public:
  //! Default constructor
  Data(bool binning = 0, unsigned int maxBins = 0, double maxW = 0.0);

  //! Default destructor
  virtual ~Data() { /* nothing */
  }

  //! Create clone
  virtual Data *Clone() const { return new Data(*this); };

  //! Create empty clone
  virtual Data *EmptyClone() const { return new Data(); };

  //! Append data sample
  virtual void Add(Data &otherSample);

  //! Add event to data sample
  virtual void PushEvent(const Event &evt);

  //! Get number of events in data sample
  virtual const std::size_t GetNEvents() const { return fEvents.size(); }

  /**! Get event
   *
   * @param id Event id
   * @return reference to event
   */
  virtual Event &GetEvent(const int id) { return fEvents.at(id); }

  /**! Get events
   *
   * @return Vector of all events
   */
  virtual std::vector<Event> &GetEvents() { return fEvents; }

  /**! Get list of data
   * A 'horizontal' list of dataPoints is obtained. Each variable
   * (e.g. m23sq, m13sq ...) is stored as a MultiDouble in ParameterList.
   * This ParameterList is used to build the FunctionTree
   * @return List of data
   */
  virtual const ParameterList &GetListOfData(std::shared_ptr<Kinematics> kin);

  std::vector<dataPoint> GetDataPoints(std::shared_ptr<Kinematics> kin) const;

  /**! Set correction value for all events.
   *
   * @param corr Correction function
   */
  virtual void ApplyCorrection(DataCorrection &corr);

  //! Remove all events outside phase-space boundaries
  virtual void ReduceToPhsp(std::shared_ptr<Kinematics> kin);

  /**! Reduce data set
   * Select first @param newSize events from full sample.
   * @param newSize
   */
  virtual void Reduce(unsigned int newSize);

  /**! Select random subset of events
   *
   * @param size Size of sub set
   * @param gen Generator
   * @return Sub set
   */
  virtual std::shared_ptr<Data> RndSubSet(std::shared_ptr<Kinematics> kin,
                                          unsigned int size,
                                          std::shared_ptr<Generator> gen);

  /**! Set resolution value for all events.
   *
   * @param res Resolution object
   */
  void SetResolution(std::shared_ptr<Resolution> res);

  /**! Set efficiency value for all events.
   *
   * @param eff Efficiency object
   */
  virtual void SetEfficiency(std::shared_ptr<Kinematics> kin,
                             std::shared_ptr<Efficiency> eff);

  //! Reset effciencies of all events
  virtual void ResetEfficiency(double e = 1.);

  //! Get maximum weight
  virtual double GetMaxWeight() const;

  /**! Reset all weights to a default value
   *
   * @param weight default weight
   */
  virtual void ResetWeights(double weight = 1.);

  //! Check of weights are stored
  virtual bool HasWeights();

  //! Clear Data
  virtual void Clear();

  /**! Write sample to file
   *
   * @param file output file name
   * @param trName name of output tree
   */
  virtual void WriteData(std::string file = "", std::string trName = "") {
    LOG(error)
        << "Data::writeData() | Base class does not provide functionality"
           "to write data to file.";
  };

  //! Obsolete?
  virtual const std::size_t GetNBins() const { return fBins.size(); }

  //! Obsolete?
  virtual const int GetBin(const int, double &, double &);

protected:
  // DataPoints are stored as 'horizontal' structure
  ParameterList dataList;

  // Vector of events
  std::vector<Event> fEvents;

  // Maximum weight of events
  double maxWeight;

  // binning
  bool fBinned;
  unsigned int fmaxBins;
  std::map<int, std::pair<double, double>> fBins;
};

//===== DATASET OPERATIONS

/** Select random sub sample of data sets
 * A hit&miss procedure is applied to the first sample to select a random set of
 * events. In case an
 * event of the first sample is added to the output sample, an event from the
 * second sample is also
 * added to the corresponding output sample. This function can be used for two
 * sample that are in sync.
 * E.g. a sample with reconstructed values and a sample with the corresponding
 * true values
 * We expect that @param out1 and @param out2 are pointers to empty samples.
 *
 * @param size size of sub sample
 * @param gen generator
 * @param in1 input sample 1
 * @param in2 input sample 2
 * @param out1 output sub sample from sample 1
 * @param out2 output sub sample from sample 2
 */
inline void rndReduceSet(std::shared_ptr<Kinematics> kin, unsigned int size,
                         std::shared_ptr<Generator> gen, Data *in1, Data *out1,
                         Data *in2 = NULL, Data *out2 = NULL);
                         
} /* namespace DataReader */
} /* namespace ComPWA */

#endif
