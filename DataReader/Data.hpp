//-------------------------------------------------------------------------------
// Copyright (c) 2013 michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Data Interface Base-Class.
/*! \class Data
 * @file Data.hpp
 * This class provides the interface to experimental data. As it is pure
 * virtual,
 * one needs at least one implementation to provide data for the other modules.
 * If
 * a new reader is derived from and fulfills this base-class, no change in other
 * modules are necessary to work with the new dataset.
 */

#ifndef DATA_HPP_
#define DATA_HPP_

#include <vector>
#include <string>
#include <memory>

#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Generator.hpp"
#include "Core/DataPoint.hpp"
#include "Core/Resolution.hpp"
#include "DataReader/DataCorrection.hpp"

namespace ComPWA {
namespace DataReader {

class Data {
public:
  //! Default constructor
  Data(bool binning = 0, unsigned int maxBins = 0, double maxW = 0.0);

  //! Default destructor
  virtual ~Data() { /* nothing */
  }

  //! Create clone
  virtual Data *Clone() const = 0;

  //! Create empty clone
  virtual Data *EmptyClone() const = 0;

  //! Append data sample
  virtual void Add(Data &otherSample);

  //! Add event to data sample
  virtual void pushEvent(const Event &evt);

  //! Get number of events in data sample
  virtual const unsigned int getNEvents() const { return fEvents.size(); }

  /**! Get event
   *
   * @param id Event id
   * @return reference to event
   */
  virtual Event &getEvent(const int id) { return fEvents.at(id); }

  /**! Get events
   *
   * @return Vector of all events
   */
  virtual std::vector<Event> getEvents() { return fEvents; }

  /**! Get list of data
   * A 'horizontal' list of dataPoints is obtained. Each variable
   * (e.g. m23sq, m13sq ...) is stored as a MultiDouble in ParameterList.
   * This ParameterList is used to build the FunctionTree
   * @return List of data
   */
  virtual const ParameterList &getListOfData();

  std::vector<dataPoint> getDataPoints() const;

  /**! Set correction value for all events.
   *
   * @param corr Correction function
   */
  virtual void applyCorrection(DataCorrection &corr);

  //! Remove all events outside phase-space boundaries
  virtual void reduceToPhsp();

  /**! Reduce data set
   * Select first @param newSize events from full sample.
   * @param newSize
   */
  virtual void reduce(unsigned int newSize);

  /**! Select random subset of events
   *
   * @param size Size of sub set
   * @param gen Generator
   * @return Sub set
   */
  virtual std::shared_ptr<Data> rndSubSet(unsigned int size,
                                          std::shared_ptr<Generator> gen);

  /**! Set resolution value for all events.
   *
   * @param res Resolution object
   */
  void setResolution(std::shared_ptr<Resolution> res);

  /**! Set efficiency value for all events.
   *
   * @param eff Efficiency object
   */
  virtual void setEfficiency(std::shared_ptr<Efficiency> eff);

  //! Reset effciencies of all events
  virtual void resetEfficiency(double e = 1.);

  //! Get maximum weight
  virtual double getMaxWeight() const;

  /**! Reset all weights to a default value
   *
   * @param weight default weight
   */
  virtual void resetWeights(double weight = 1.);

  //! Check of weights are stored
  virtual bool hasWeights();

  //! Clear Data
  virtual void Clear();

  /**! Write sample to file
   *
   * @param file output file name
   * @param trName name of output tree
   */
  virtual void writeData(std::string file = "", std::string trName = "") = 0;

  //! Obsolete?
  virtual const unsigned int getNBins() const { return fBins.size(); }

  //! Obsolete?
  virtual const int getBin(const int, double &, double &);

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
void rndReduceSet(unsigned int size, std::shared_ptr<Generator> gen, Data *in1,
                  Data *out1, Data *in2 = NULL, Data *out2 = NULL);
} /* namespace DataReader */
} /* namespace ComPWA */

#endif
