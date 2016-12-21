/*
 * IteratorTestApp.cpp
 *
 *  Created on: Mar 03, 2016
 *      Author: michel
 */
#include "Core/Parameter.hpp"
#include <sstream>
#include <vector>
#include <boost/iterator/filter_iterator.hpp>

class my_Resonance {
public:
  my_Resonance(bool in_enabled, unsigned int in_id):enabled(in_enabled),id(in_id){

  };
  ~my_Resonance(){};
  bool enabled;
  unsigned int id;
};

struct is_enabled {
  bool operator()(my_Resonance x) { return x.enabled; }
};

int main(int argc, char** argv){

  std::vector<my_Resonance> resos;
  resos.push_back(my_Resonance(true,1));
  resos.push_back(my_Resonance(false,2));
  resos.push_back(my_Resonance(true,3));
  resos.push_back(my_Resonance(true,4));
  resos.push_back(my_Resonance(false,5));

  // Example using filter_iterator
  typedef boost::filter_iterator<is_enabled, std::vector<my_Resonance>::iterator>
    FilterIter;

  is_enabled predicate;
  FilterIter filter_iter_first(predicate, resos.begin(), resos.end());
  FilterIter filter_iter_last(predicate, resos.end(), resos.end());

  for(FilterIter i = filter_iter_first; i!=filter_iter_last; ++i){
    std::cout << (*i).id << std::endl;
  }

  return 0;
}

