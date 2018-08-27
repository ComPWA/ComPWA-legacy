// Copyright (c) 2014 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PROGRESSBAR_HPP_
#define PROGRESSBAR_HPP_

#include <ctime>
#include <iostream>

namespace ComPWA {

class ProgressBar {
public:
  ProgressBar() : hasStarted(0){};

  ~ProgressBar() { std::cout << std::endl; };

  ProgressBar(std::size_t size, int update = 1);
  
  /// indicate the next step in process
  void next();

protected:
  double timeRemaining();
  double timePassed();
  time_t endTime();
  void update();

  std::size_t numEvents;
  int updateInterval;
  bool hasStarted;
  time_t startTime;

  double currentPercent;
  unsigned int currentEvent;
  double lastUpdate;
};

} // ns::ComPWA

#endif
