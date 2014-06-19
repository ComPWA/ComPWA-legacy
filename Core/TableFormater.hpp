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

using namespace std;


class TableFormater {

public:
	TableFormater(std::ostream* output): curRow(0),curCol(0),totalWidth(0){
		out = output;
	};
	void delim(){
		*out<<"|";
		for(unsigned int i=0;i<totalWidth-1; i++) *out << "-" ;
		*out<<"|"<<endl;
	}
	void footer(){ delim(); }
	void header(){
		delim();
		for(unsigned int i=0;i<columnTitle.size();i++) *out << "| "<<std::setw(columnWidth[i])<<columnTitle[i]<<" ";
		*out<<"|"<<std::endl;
		delim();
	};
	void addColumn(std::string title, unsigned int fixlength=999){
		unsigned int length;
		if(fixlength!=999) length = fixlength;
		else length=title.length();
		columnWidth.push_back(length);
		totalWidth+=length+3;
		columnTitle.push_back(title);
	};

	void trimString(std::string& src){
		//remove  all final zeros
		char chr = '0';
//		std::string::size_type pos =  src.find_first_not_of(chr,0);
//		if(pos-1 > 0)
//			src.erase(0,pos-1);
		std::string::size_type pos2 =  src.find_last_not_of(chr,src.length());
		if(pos2+1 < src.length() )
			src.erase(pos2+1,src.length());
	}
	TableFormater& operator<<(DoubleParameter in){
		std::stringstream errStr;
		if(curCol==0) *out << "| ";
		else *out << " | ";
		if(in.HasError()){
			std::string tmp;
			if(in.GetErrorType()==ErrorType::SYM){
				unsigned int halfWidth = (unsigned int)(columnWidth[curCol])/2;//divide column width
				*out << std::setw(halfWidth) << in.GetValue();
				tmp = " +-"+std::to_string((long double) in.GetError()->GetError()); trimString(tmp);
				*out << std::setw(halfWidth) << tmp;
			}
			if(in.GetErrorType()==ErrorType::ASYM){
				unsigned int w = (unsigned int)(columnWidth[curCol])/3;//take 1/3 of column width
				std::shared_ptr<ParError<double>> err = in.GetError();
				tmp = std::to_string((long double) in.GetValue()); trimString(tmp);
				*out << std::setw(w) << tmp;
				tmp = "+"+std::to_string((long double) err->GetErrorHigh());trimString(tmp);
				*out << std::setw(w) << tmp;
				tmp = "-"+std::to_string((long double) err->GetErrorLow());trimString(tmp);
				*out << std::setw(w) << tmp;
			}
		} else {
			*out << std::setw(columnWidth[curCol]) << in.GetValue();
		}
		curCol++;
		if(curCol==columnWidth.size()) {
			*out<<" |"<<std::endl;
			curRow++; curCol=0;
		}
		return *this;
	};
	template<typename T> TableFormater& operator<<(T in){
		if(curCol==0) *out << "| ";
		else *out << " | ";
		*out << std::setw(columnWidth[curCol]) << in;
		curCol++;
		if(curCol==columnWidth.size()) {
			*out<<" |"<<std::endl;
			curRow++; curCol=0;
			//			delim();
		}
		return *this;
	};
private:
	ostream* out;
	std::vector<unsigned int> columnWidth;
	std::vector<std::string> columnTitle;

	unsigned int curRow;
	unsigned int curCol;
	unsigned int totalWidth;
protected:
};


#endif /* TABLEFORMATER_CXX_ */
