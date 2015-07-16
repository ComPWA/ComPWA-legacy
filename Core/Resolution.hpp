/*
 * Resolution.hpp
 *
 *  Created on: Jul 14, 2015
 *      Author: weidenka
 */

#ifndef CORE_RESOLUTION_HPP_
#define CORE_RESOLUTION_HPP_


class Resolution
{
public:
	Resolution() {};
	virtual ~Resolution() {};
	virtual void resolution(Event& ev) = 0;
	virtual std::vector<double> resolution(std::vector<double> v) = 0;

};

class ZeroResolution : public Resolution
{
public:
	ZeroResolution() {};
	virtual ~ZeroResolution() {};
	virtual void resolution(Event& ev) { return; }
	virtual std::vector<double> resolution(std::vector<double> v) {
		std::vector<double> offset(v.size(),0);
		for(int i=0; i<v.size(); i++) v.at(i)+=offset.at(i);
	}

};

#endif /* CORE_RESOLUTION_HPP_ */
