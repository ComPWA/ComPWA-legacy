/*
 * progressBar.hpp
 *
 *  Created on: Aug 6, 2014
 *      Author: weidenka
 */

#ifndef PROGRESSBAR_HPP_
#define PROGRESSBAR_HPP_

#include <ctime>
#include <iostream>

namespace ComPWA {

class progressBar
{
public:
	progressBar() : hasStarted(0){ };
	~progressBar() { std::cout<<std::endl; };
	progressBar(int size, int update=1);
	void nextEvent();
protected:
	double timeRemaining();
	double timePassed();
	time_t endTime();
	void update();

	int numEvents;
	int updateInterval;
	bool hasStarted;
	time_t startTime;

	double currentPercent;
	int currentEvent;
	double lastUpdate;

};

} /* namespace ComPWA */

#endif /* PROGRESSBAR_HPP_ */
