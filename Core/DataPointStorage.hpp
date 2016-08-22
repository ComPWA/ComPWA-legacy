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
  unsigned int number_of_events_;
  unsigned int expected_number_of_events_;
  unsigned int number_of_variables_;
  // the map key is a the variable or storage index. then we have one data list for each variable
  std::map<unsigned int, DataList> data_storage_;

  DataPointStorage();

public:
  static DataPointStorage& Instance() {
    static DataPointStorage instance;
    return instance;
  }
  virtual ~DataPointStorage();

  DataPointStorage(DataPointStorage const&) = delete;
  void operator=(DataPointStorage const&) = delete;

  unsigned int getNumberOfEvents() const;

  void layoutDataStorageStructure(unsigned int expected_number_of_events_, const Event& evt);

  void clearStorage();

  void addEvent(const Event& evt);
  void addDataPoint(const dataPoint& dp);

  const DataList& getDataList(
      unsigned int storage_index) const;
};

} /* namespace ComPWA */

#endif /* _DATAPOINTSTORAGE_HPP_ */
