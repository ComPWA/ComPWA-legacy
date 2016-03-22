/*
 * TablePrinter.cxx
 *
 *  Created on: Jan 10, 2014
 *      Author: weidenka
 */

#ifndef TABLEFORMATER_CXX_
#define TABLEFORMATER_CXX_
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <string>

#include "Core/ParameterList.hpp"

using namespace std;

class TableFormater
{
public:
	TableFormater(std::ostream* output) :
		out(output),curRow(0),curCol(0),totalWidth(0) {
		sep = "|";
		pm = "+-";
		firstSep = sep;
		lastSep = sep;
	};
	virtual ~TableFormater() { }
	virtual void Reset();
	virtual void delim();
	virtual void footer();
	virtual void header();
	virtual void addColumn(std::string title, unsigned int fixlength=999);
	void trimString(std::string& src);
	TableFormater& operator<<(DoubleParameter in);
	template<typename T> TableFormater& operator<<(T in){
		if(curCol==0) *out << firstSep+" ";
		else *out << " "<<sep<<" ";
		*out << std::setw(columnWidth[curCol]) << in;
		curCol++;
		if(curCol==columnWidth.size()) {
			*out<<" "<<lastSep<<std::endl;
			curRow++; curCol=0;
			//			delim();
		}
		return *this;
	};

protected:
	ostream* out;
	std::vector<unsigned int> columnWidth;
	std::vector<std::string> columnTitle;

	unsigned int curRow;
	unsigned int curCol;
	unsigned int totalWidth;
	std::string sep;
	std::string firstSep;
	std::string lastSep;
	std::string pm;
};

class TexTableFormater : public TableFormater
{
public:
	TexTableFormater(std::ostream* output);
	virtual ~TexTableFormater() { }
	virtual void footer();
	virtual void header();
	virtual void delim();
};
#endif /* TABLEFORMATER_CXX_ */
