 
                
#include "DataReader/Data.hpp"
namespace ComPWA {
namespace DataReader {

Data::Data(bool binning, unsigned int maxBins, double maxW)
    : maxWeight(maxW), fBinned(binning), fmaxBins(maxBins) {}

void rndReduceSet(unsigned int size, std::shared_ptr<Generator> gen, Data *in1,
                  Data *out1, Data *in2, Data *out2) {
  if (!in1)
    throw std::runtime_error("rndSubSet() | No input data set!");
  if (!out1)
    throw std::runtime_error("rndSubSet() | No output data set!");
  if (out1->GetNEvents())
    throw std::runtime_error("rndSubSet() | First output sample not empty!");
  if (in2) {
    if (in1->GetNEvents() != in2->GetNEvents())
      throw std::runtime_error(
          "rndSubSet() | Samples have different event count!");
    if (!out2)
      throw std::runtime_error("rndSubSet() | Second output set is NULL!");
    if (out2->GetNEvents())
      throw std::runtime_error("rndSubSet() | Second output sample not empty!");
  }

  unsigned int totalSize = in1->GetNEvents();
  unsigned int newSize = totalSize;
  /* 1th method: new Sample has exact size, but possibly events are added twice.
   * We would have to store all used events in a vector and search the vector at
   * every event -> slow */
  /*unsigned int t=0;
   unsigned int d=0;
   while(t<newSize){
   d = (unsigned int) gen->getUniform()*totalSize;
   newSample->pushEvent(fEvents[d]);
   t++;
   }*/

  /* 2nd method: events are added once only, but total size of sample varies
   * with sqrt(N) */
  for (unsigned int i = 0; i < totalSize;
       i++) { // count how many events are not within PHSP
    dataPoint point;
    try {
      point = dataPoint(in1->GetEvent(i));
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      newSize--;
      continue;
    }
    //		dataPoint point(in1->getEvent(i));
    //		if(!Kinematics::instance()->isWithinPhsp(point)) newSize--;
  }
  double threshold = (double)size / newSize; // calculate threshold
  for (unsigned int i = 0; i < totalSize; i++) {
    dataPoint point;
    try {
      point = dataPoint(in1->GetEvent(i));
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }
    //		dataPoint point(in1->getEvent(i)); //use first sample for
    //hit&miss
    //		if(!Kinematics::instance()->isWithinPhsp(point)) continue;
    if (gen->GetUniform(0,1) < threshold) {
      out1->PushEvent(in1->GetEvent(i));
      // write second sample if event from first sample were accepted
      if (in2)
        out2->PushEvent(in2->GetEvent(i));
    }
  }
  if (out2)
    assert(out1->GetNEvents() == out2->GetNEvents());

  LOG(debug) << "DataReader::rndReduceSet() | sample size reduced to "
             << out1->GetNEvents();
  return;
}

std::shared_ptr<Data> Data::RndSubSet(unsigned int size,
                                      std::shared_ptr<Generator> gen) {
  std::shared_ptr<Data> out(this->EmptyClone());
  rndReduceSet(size, gen, this, out.get());
  return out;
}

void Data::ResetWeights(double w) {
  for (unsigned int i = 0; i < fEvents.size(); i++)
    fEvents.at(i).SetWeight(w);
  maxWeight = w;
  return;
}

double Data::GetMaxWeight() const { return maxWeight; }

void Data::ReduceToPhsp() {
  std::vector<Event> tmp;
  LOG(info) << "Data::reduceToPhsp() | "
               "Remove all events outside PHSP boundary from data sample.";

  for (unsigned int evt = 0; evt < fEvents.size(); evt++) {
    dataPoint point;
    try {
      point = dataPoint(fEvents.at(evt));
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }
    //		if(Kinematics::instance()->isWithinPhsp(p))
    tmp.push_back(fEvents.at(evt));
  }
  LOG(info) << "Data::reduceToPhsp() | " << tmp.size() << " from "
            << fEvents.size() << "("
            << ((double)tmp.size()) / fEvents.size() * 100 << "%) were kept.";
  fEvents = tmp;
  return;
}
void Data::ResetEfficiency(double e) {
  for (unsigned int evt = 0; evt < fEvents.size(); evt++) {
    fEvents.at(evt).SetEfficiency(e);
  }
}
void Data::Reduce(unsigned int newSize) {
  if (newSize >= fEvents.size()) {
    LOG(error)
        << "RooReader::reduce() requested size too large, cant reduce sample!";
    return;
  }
  fEvents.resize(newSize);
}

void Data::SetEfficiency(std::shared_ptr<Efficiency> eff) {
  for (unsigned int evt = 0; evt < fEvents.size(); evt++) {
    dataPoint point;
    try {
      point = dataPoint(fEvents.at(evt));
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }
    //		dataPoint point(fEvents.at(evt));
    double val = eff->Evaluate(point);
    fEvents.at(evt).SetEfficiency(val);
  }
}
void Data::Clear() { fEvents.clear(); }

bool Data::HasWeights() {
  bool has = 0;
  for (unsigned int evt = 0; evt < fEvents.size(); evt++) {
    if (fEvents.at(evt).GetWeight() != 1.) {
      has = 1;
      break;
    }
  }
  return has;
}

const ParameterList &Data::GetListOfData() {
  // dataList already filled? return filled one
  if (dataList.GetNParameter() != 0)
    return dataList;

  int size = Kinematics::Instance()->GetNVars();
  std::vector<std::vector<double>> data(size, std::vector<double>());
  std::vector<double> eff;
  eff.reserve(fEvents.size());
  std::vector<double> weight;
  weight.reserve(fEvents.size());

  auto itr = fEvents.begin();
  for (; itr != fEvents.end(); ++itr) {
    dataPoint point;
    try {
      point = dataPoint(*itr);
    } catch (BeyondPhsp &ex) {
      continue;
    }
    eff.push_back(point.GetEfficiency());
    weight.push_back(point.getWeight());
    for (int i = 0; i < size; ++i)
      data.at(i).push_back(point.GetValue(i));
  }

  // Add data vector to ParameterList
  for (int i = 0; i < size; ++i) {
    //std::shared_ptr<MultiDouble> tmp(new MultiDouble(
    //  Kinematics::instance()->GetVarNames().at(i), data.at(i)));
    std::shared_ptr<MultiDouble> tmp(new MultiDouble(
        "", data.at(i)));
    dataList.AddParameter(tmp);
  }

  // Adding efficiency at the end
  dataList.AddParameter(
      std::shared_ptr<MultiDouble>(new MultiDouble("Efficiency", eff)));

  // Adding weight at the end
  dataList.AddParameter(
      std::shared_ptr<MultiDouble>(new MultiDouble("Weight", weight)));

  return dataList;
}

std::vector<dataPoint> Data::GetDataPoints() const {
  std::vector<dataPoint> vecPoint;
  for (int i = 0; i < fEvents.size(); i++) {
    dataPoint point;
    try {
      point = dataPoint(fEvents.at(i));
    } catch (BeyondPhsp &ex) { // event outside phase, remove
      continue;
    }
    vecPoint.push_back(point);
  }
  return vecPoint;
}

void Data::SetResolution(std::shared_ptr<Resolution> res) {
  for (int i = 0; i < fEvents.size(); i++)
    res->resolution(fEvents.at(i));
}

void Data::Add(Data &otherSample) {
  std::vector<Event> otherEvents = otherSample.GetEvents();
  fEvents.insert(fEvents.end(), otherEvents.begin(), otherEvents.end());
  if (otherSample.GetMaxWeight() > maxWeight)
    maxWeight = otherSample.GetMaxWeight();
  return;
}

void Data::ApplyCorrection(DataCorrection &corr) {
  double sumWeightSq = 0;
  for (int i = 0; i < fEvents.size(); i++) {
    double w = corr.getCorrection(fEvents.at(i));
    if (w < 0)
      throw std::runtime_error("Data::applyCorrection() | "
                               "Negative weight!");
    sumWeightSq += w * w;
    double oldW = fEvents.at(i).GetWeight();
    if (w * oldW > maxWeight)
      maxWeight = w * oldW;
    fEvents.at(i).SetWeight(w * oldW);
  }
  LOG(info) << "Data::applyCorrection() | "
               "Sample corrected! Sum of weights squared is "
            << sumWeightSq;
  return;
}

const int Data::GetBin(const int i, double &m12, double &weight) {
  if (!fBinned)
    return -1;

  m12 = fBins[i].first;
  weight = fBins[i].second;

  return 1;
}

//! Add event to data sample
void Data::PushEvent(const Event &evt) {
  fEvents.push_back(evt);
  if (evt.GetWeight() > maxWeight)
    maxWeight = evt.GetWeight();
}
    
} /* namespace DataReader */
} /* namespace ComPWA */
