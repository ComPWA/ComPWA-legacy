// Copyright (c) 2014 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef PROGRESSBAR_HPP_
#define PROGRESSBAR_HPP_

#include <ctime>
#include <iostream>

namespace ComPWA {

class progressBar {
public:
  progressBar() : hasStarted(0){};
  ~progressBar() { std::cout << std::endl; };
  progressBar(std::size_t size, int update = 1);
  void nextEvent();

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
  int currentEvent;
  double lastUpdate;
};

} /* namespace ComPWA */

#endif /* PROGRESSBAR_HPP_ */
