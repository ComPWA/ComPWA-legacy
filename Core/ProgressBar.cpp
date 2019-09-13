// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/ProgressBar.hpp"

#include <chrono>
#include <iostream>

namespace ComPWA {

ProgressBar::ProgressBar() : ProgressBar(0, 0) {}

ProgressBar::ProgressBar(std::size_t size, int update)
    : numEvents(size), updateInterval(update), hasStarted(0), startTime(0),
      currentPercent(0), currentEvent(0), lastUpdate(0) {
  if (update == 0)
    updateInterval = 1;
}

ProgressBar::~ProgressBar() { std::cout << std::endl; }

void ProgressBar::next(size_t increment) {
  if (!hasStarted) {
    lastUpdate = 0;
    currentEvent = 0;
    time(&startTime);
    hasStarted = 1;
    update();
    fflush(stdout);
  }
  currentEvent += increment;
  if (currentEvent > numEvents)
    currentEvent = numEvents;
  if ((int)((timePassed() - lastUpdate)) > updateInterval)
    update();
  if (currentEvent == numEvents) {
    update();
    std::cout << std::endl;
  }
}
void ProgressBar::update() {
  //	std::cout<<"update()"<<std::endl;
  //	std::cout<<timePassed()<<" "<<lastUpdate<<"
  //"<<(int)((timePassed()-lastUpdate))<<std::endl;
  currentPercent = (double)currentEvent / numEvents * 100;
  int nStars = ((int)(currentPercent / 10 + 0.5));
  char buf[10];
  int i = 0;
  for (; i < nStars; i++)
    buf[i] = '*';
  for (; i < 10; i++)
    buf[i] = '-';
  //	for( ; i<15; i++) buf[i]=' ';
  time_t estEndTime = endTime();
  char timebuf[10];
  strftime(timebuf, 10, "%I:%M%p", localtime(&estEndTime));
  // the escape sequence \e[?25l switches off the courser and \e[?25h switch it
  // on again
  printf("\e[?25l \r %4.2fmin [ %.*s %4.2f%% ] %4.2fmin , end time: %7s        "
         "       \e[?25h",
         timePassed() / 60, 10, buf, currentPercent, timeRemaining() / 60,
         timebuf);
  fflush(stdout);

  lastUpdate = timePassed();
  return;
}

double ProgressBar::timeRemaining() {
  if (currentPercent == 0)
    return 0.0;
  return timePassed() * (1 / currentPercent * 100 - 1);
}
double ProgressBar::timePassed() {
  time_t currentTime;
  time(&currentTime);
  return difftime(currentTime, startTime);
}

time_t ProgressBar::endTime() {
  time_t currentTime;
  time(&currentTime);
  time_t estEndTime = currentTime + timeRemaining();
  return estEndTime;
}

} // namespace ComPWA
