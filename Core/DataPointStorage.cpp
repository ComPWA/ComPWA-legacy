/*
 * DataPointStorage.cpp
 *
 *  Created on: Aug 9, 2016
 *      Author: steve
 */

#include "Core/DataPointStorage.hpp"

#include "Core/Event.hpp"
#include "Core/DataPoint.hpp"

namespace ComPWA {

DataPointStorage::DataPointStorage() :
    number_of_events_(0) {
}

DataPointStorage::~DataPointStorage() {
}

unsigned int DataPointStorage::getNumberOfEvents() const {
  return number_of_events_;
}

void DataPointStorage::layoutDataStorageStructure(
    unsigned int expected_number_of_events, const Event& evt) {
  // reserve storage to have fast filling later on
  dataPoint dp(evt);
  expected_number_of_events_ = expected_number_of_events;
  for (unsigned int i = 0; i < dp.size(); ++i) {
    data_storage_[i].reserve(expected_number_of_events_);
  }
}

void DataPointStorage::clearStorage() {
  number_of_events_ = 0;
  for(auto& data_list : data_storage_) {
    data_list.second.clear();
  }
}

void DataPointStorage::addEvent(const Event& evt) {
  dataPoint dp(evt);

  addDataPoint(dp);
}

void DataPointStorage::addDataPoint(const dataPoint& dp) {
  for (unsigned int i = 0; i < dp.size(); ++i) {
    data_storage_[i].push_back(dp.getVal(i));
  }
  ++number_of_events_;
}

const DataList& DataPointStorage::getDataList(
    unsigned int storage_index) const {
  return data_storage_.at(storage_index);
}

} /* namespace ComPWA */
