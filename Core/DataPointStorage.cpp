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

DataPointStorage::DataPointStorage() {
}

DataPointStorage::~DataPointStorage() {
}

unsigned int DataPointStorage::getNumberOfEvents(
    unsigned int storage_index) const {
  return number_of_events_.at(storage_index);
}

void DataPointStorage::layoutDataStorageStructure(unsigned int storage_index,
    unsigned int expected_number_of_events, const Event& evt) {
  // reserve storage to have fast filling later on
  dataPoint dp(evt);
  expected_number_of_events_[storage_index] = expected_number_of_events;
  number_of_events_[storage_index] = 0;
  for (unsigned int i = 0; i < dp.size(); ++i) {
    data_storage_[storage_index][i].reserve(
        expected_number_of_events_.at(storage_index));
  }
}

void DataPointStorage::layoutDataStorageStructure(unsigned int storage_index,
    unsigned int expected_number_of_events, const dataPoint& dp) {
  // reserve storage to have fast filling later on
  expected_number_of_events_[storage_index] = expected_number_of_events;
  number_of_events_[storage_index] = 0;
  for (unsigned int i = 0; i < dp.size(); ++i) {
    data_storage_[storage_index][i].reserve(
        expected_number_of_events_.at(storage_index));
  }
}

void DataPointStorage::clearStorage() {
  for (auto& entry : number_of_events_) {
    entry.second = 0;
  }
  for (auto& single_storage : data_storage_) {
    for (auto& data_list : single_storage.second) {
      data_list.second.clear();
    }
  }
}

void DataPointStorage::addEvent(unsigned storage_index, const Event& evt) {
  dataPoint dp(evt);

  addDataPoint(storage_index, dp);
}

void DataPointStorage::addDataPoint(unsigned storage_index,
    const dataPoint& dp) {
  auto& chosen_storage = data_storage_.at(storage_index);
  for (unsigned int i = 0; i < dp.size(); ++i) {
    chosen_storage[i].push_back(dp.getVal(i));
  }
  ++(number_of_events_.at(storage_index));
}

const DataList& DataPointStorage::getDataList(unsigned int storage_index,
    unsigned int variable_index) const {
  return data_storage_.at(storage_index).at(variable_index);
}

} /* namespace ComPWA */
