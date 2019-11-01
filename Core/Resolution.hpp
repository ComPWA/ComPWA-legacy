// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef CORE_RESOLUTION_HPP_
#define CORE_RESOLUTION_HPP_

namespace ComPWA {

class Resolution {
public:
  Resolution(){};
  virtual ~Resolution(){};
  virtual void resolution(Event &ev) = 0;
  virtual std::vector<double> resolution(std::vector<double> v) = 0;
};

class ZeroResolution : public Resolution {
public:
  ZeroResolution(){};
  virtual ~ZeroResolution(){};
  virtual void resolution(Event &ev) { return; }
  virtual std::vector<double> resolution(std::vector<double> v) {
    std::vector<double> offset(v.size(), 0);
    for (unsigned int i = 0; i < v.size(); i++)
      v.at(i) += offset.at(i);
    return v;
  }
};
} // namespace ComPWA

#endif
