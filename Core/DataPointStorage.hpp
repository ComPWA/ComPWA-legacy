/*
 * DataPointStorage.hpp
 *
 *  Created on: Aug 9, 2016
 *      Author: steve
 */

#ifndef _DATAPOINTSTORAGE_HPP_
#define _DATAPOINTSTORAGE_HPP_

#include <vector>
#include <map>

namespace ComPWA {

class Event;
class dataPoint;

typedef std::vector<double> DataList;

class DataPointStorage {
  std::map<unsigned int, unsigned int> number_of_events_;
  std::map<unsigned int, unsigned int> expected_number_of_events_;
  std::map<unsigned int, unsigned int> number_of_variables_;
  // the first map key is a the storage index (phsp or data).
  // the second map key is the variable index. then we have one data list for each variable
  std::map<unsigned int, std::map<unsigned int, DataList> > data_storage_;

  DataPointStorage();

public:
  static DataPointStorage& Instance() {
    static DataPointStorage instance;
    return instance;
  }
  virtual ~DataPointStorage();

  DataPointStorage(DataPointStorage const&) = delete;
  void operator=(DataPointStorage const&) = delete;

  unsigned int getNumberOfEvents(unsigned int storage_index) const;

  void layoutDataStorageStructure(unsigned int storage_index,
      unsigned int expected_number_of_events_, const Event& evt);
  void layoutDataStorageStructure(unsigned int storage_index,
      unsigned int expected_number_of_events_, const dataPoint& dp);

  void clearStorage();

  void addEvent(unsigned int storage_index, const Event& evt);
  void addDataPoint(unsigned int storage_index, const dataPoint& dp);

  const DataList& getDataList(unsigned int storage_index,
      unsigned int variable_index) const;
};

} /* namespace ComPWA */

#endif /* _DATAPOINTSTORAGE_HPP_ */
