// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.
#include "Data/Data.hpp"

namespace ComPWA {
namespace Data {

Data::Data(bool binning, unsigned int maxBins, double maxW)
    : MaximumWeight(maxW), fBinned(binning), fmaxBins(maxBins) {}

void rndReduceSet(std::shared_ptr<ComPWA::Kinematics> kin, unsigned int size,
                  std::shared_ptr<ComPWA::Generator> gen, Data *in1, Data *out1,
                  Data *in2, Data *out2) {
  if (!in1)
    throw std::runtime_error("rndSubSet() | No input data set!");
  if (!out1)
    throw std::runtime_error("rndSubSet() | No output data set!");
  if (out1->numEvents())
    throw std::runtime_error("rndSubSet() | First output sample not empty!");
  if (in2) {
    if (in1->numEvents() != in2->numEvents())
      throw std::runtime_error(
          "rndSubSet() | Samples have different event count!");
    if (!out2)
      throw std::runtime_error("rndSubSet() | Second output set is NULL!");
    if (out2->numEvents())
      throw std::runtime_error("rndSubSet() | Second output sample not empty!");
  }

  unsigned int totalSize = in1->numEvents();
  unsigned int newSize = totalSize;
  /* 1th method: new Sample has exact size, but possibly events are added twice.
   * We would have to store all used events in a vector and search the vector at
   * every event -> slow */
  /*unsigned int t=0;
   unsigned int d=0;
   while(t<newSize){
   d = (unsigned int) gen->getUniform()*totalSize;
   newSample->pushEvent(fEvents.at(d));
   t++;
   }*/

  /* 2nd method: events are added once only, but total size of sample varies
   * with sqrt(N) */
  for (unsigned int i = 0; i < totalSize;
       i++) { // count how many events are not within PHSP
    ComPWA::DataPoint point;
    try {
      kin->convert(in1->event(i), point);
    } catch (ComPWA::BeyondPhsp &ex) { // event outside phase, remove
      newSize--;
      continue;
    }
    //    dataPoint point(in1->getEvent(i));
    //    if(!Kinematics::instance()->isWithinPhsp(point)) newSize--;
  }
  double threshold = (double)size / newSize; // calculate threshold
  for (unsigned int i = 0; i < totalSize; i++) {
    ComPWA::DataPoint point;
    try {
      kin->convert(in1->event(i), point);
    } catch (ComPWA::BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }
    //    dataPoint point(in1->getEvent(i)); //use first sample for
    // hit&miss
    //    if(!Kinematics::instance()->isWithinPhsp(point)) continue;
    if (gen->uniform(0, 1) < threshold) {
      out1->add(in1->event(i));
      // write second sample if event from first sample were accepted
      if (in2)
        out2->add(in2->event(i));
    }
  }
  if (out2)
    assert(out1->numEvents() == out2->numEvents());

  LOG(DEBUG) << "DataReader::rndReduceSet() | sample size reduced to "
             << out1->numEvents();
  return;
}

std::shared_ptr<Data> Data::rndSubSet(std::shared_ptr<Kinematics> kin,
                                      unsigned int size,
                                      std::shared_ptr<Generator> gen) {
  std::shared_ptr<Data> out(new Data());
  rndReduceSet(kin, size, gen, this, out.get());
  return out;
}

void Data::resetWeights(double w) {
  for (unsigned int i = 0; i < Events.size(); ++i)
    Events.at(i).setWeight(w);
  MaximumWeight = w;
  return;
}

double Data::maximumWeight() const { return MaximumWeight; }

void Data::reduceToPhsp(std::shared_ptr<Kinematics> kin) {
  std::vector<Event> tmp;
  LOG(INFO) << "Data::reduceToPhsp() | "
               "Remove all events outside PHSP boundary from data sample.";

  for (unsigned int evt = 0; evt < Events.size(); ++evt) {
    DataPoint point;
    try {
      kin->convert(Events.at(evt), point);
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }
    tmp.push_back(Events.at(evt));
  }
  LOG(INFO) << "Data::reduceToPhsp() | " << tmp.size() << " from "
            << Events.size() << "("
            << ((double)tmp.size()) / Events.size() * 100 << "%) were kept.";
  Events = tmp;
  return;
}

void Data::resetEfficiency(double e) {
  for (unsigned int evt = 0; evt < Events.size(); ++evt) {
    Events.at(evt).setEfficiency(e);
  }
}

void Data::reduce(unsigned int newSize) {
  if (newSize >= Events.size()) {
    LOG(ERROR)
        << "RooReader::reduce() requested size too large, cant reduce sample!";
    return;
  }
  Events.resize(newSize);
}

void Data::setEfficiency(std::shared_ptr<Kinematics> kin,
                         std::shared_ptr<Efficiency> eff) {
  for (unsigned int evt = 0; evt < Events.size(); ++evt) {
    DataPoint point;
    try {
      kin->convert(Events.at(evt), point);
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }
    //    dataPoint point(fEvents.at(evt));
    double val = eff->evaluate(point);
    Events.at(evt).setEfficiency(val);
  }
}
void Data::clear() { Events.clear(); }

bool Data::hasWeights() {
  bool has = 0;
  for (unsigned int evt = 0; evt < Events.size(); ++evt) {
    if (Events.at(evt).weight() != 1.) {
      has = 1;
      break;
    }
  }
  return has;
}

const ComPWA::ParameterList &
Data::dataList(std::shared_ptr<ComPWA::Kinematics> kin) {
  // dataList already filled? return filled one
  if (DataList.numParameters() != 0)
    return DataList;

  int size = kin->numVariables();
  std::vector<std::vector<double>> data(size, std::vector<double>());
  std::vector<double> eff;
  eff.reserve(Events.size());
  std::vector<double> weight;
  weight.reserve(Events.size());

  auto itr = Events.begin();
  for (; itr != Events.end(); ++itr) {
    DataPoint point;
    try {
      kin->convert(*itr, point);
    } catch (BeyondPhsp &ex) {
      continue;
    }
    eff.push_back(point.efficiency());
    weight.push_back(point.weight());
    for (int i = 0; i < size; ++i)
      data.at(i).push_back(point.value(i));
  }

  // Add data vector to ParameterList
  for (auto i : data)
    DataList.addValue(MDouble("", i));
  // Adding efficiency at the end
  DataList.addValue(MDouble("Efficiency", eff));
  // Adding weight at the end
  DataList.addValue(MDouble("Weight", weight));

  return DataList;
}

std::vector<ComPWA::DataPoint>
Data::dataPoints(std::shared_ptr<ComPWA::Kinematics> kin) const {
  std::vector<DataPoint> vecPoint;
  for (unsigned int i = 0; i < Events.size(); ++i) {
    DataPoint point;
    try {
      kin->convert(Events.at(i), point);
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }
    vecPoint.push_back(point);
  }
  return vecPoint;
}

void Data::setResolution(std::shared_ptr<Resolution> res) {
  for (unsigned int i = 0; i < Events.size(); ++i)
    res->resolution(Events.at(i));
}

void Data::append(Data &otherSample) {
  std::vector<Event> otherEvents = otherSample.events();
  Events.insert(Events.end(), otherEvents.begin(), otherEvents.end());
  if (otherSample.maximumWeight() > MaximumWeight)
    MaximumWeight = otherSample.maximumWeight();
  return;
}

void Data::applyCorrection(DataCorrection &corr) {
  double sumWeightSq = 0;
  for (auto &Event : Events) {
    double w = corr.correction(Event);
    if (w < 0)
      throw std::runtime_error("Data::applyCorrection() | "
                               "Negative weight!");
    sumWeightSq += w * w;
    double oldW = Event.weight();
    if (w * oldW > MaximumWeight)
      MaximumWeight = w * oldW;
    Event.setWeight(w * oldW);
  }
  LOG(INFO) << "Data::applyCorrection() | "
               "Sample corrected! Sum of weights squared is "
            << sumWeightSq;
  return;
}

const int Data::bin(const int i, double &m12, double &weight) {
  if (!fBinned)
    return -1;

  m12 = fBins[i].first;
  weight = fBins[i].second;

  return 1;
}

void Data::add(const Event &evt) {
  Events.push_back(evt);
  if (evt.weight() > MaximumWeight)
    MaximumWeight = evt.weight();
}

void Data::add(const std::vector<Event> &evts) {
  for (auto const &evt : evts)
    add(evt);
}

} // namespace DataReader
} // namespace ComPWA
