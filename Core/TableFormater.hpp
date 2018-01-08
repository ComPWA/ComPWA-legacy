// Copyright (c) 2015 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef TABLEFORMATER_HPP_
#define TABLEFORMATER_HPP_

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <string>

#include "Core/ParameterList.hpp"

namespace ComPWA {

class TableFormater {
public:
  TableFormater(std::ostream *output)
      : OutputStream(output), CurRow(0), CurCol(0), TotalWidth(0) {
    sep = "|";
    pm = "+-";
    firstSep = sep;
    lastSep = sep;
  };
  
  virtual ~TableFormater() {}
  
  virtual void reset();
  
  virtual void delim();
  
  virtual void footer();
  
  virtual void header();
  
  virtual void addColumn(std::string title, unsigned int fixlength = 999);
  
  void trimString(std::string &src);
  
  TableFormater &operator<<(FitParameter in);
  
  template <typename T> TableFormater &operator<<(T in) {
    if (CurCol == 0)
      *OutputStream << firstSep + " ";
    else
      *OutputStream << " " << sep << " ";
    *OutputStream << std::setw(ColumnWidth.at(CurCol)) << in;
    CurCol++;
    if (CurCol == ColumnWidth.size()) {
      *OutputStream << " " << lastSep << std::endl;
      CurRow++;
      CurCol = 0;
      //			delim();
    }
    return *this;
  };

protected:
  std::ostream *OutputStream;
  std::vector<unsigned int> ColumnWidth;
  std::vector<std::string> ColumnTitle;

  unsigned int CurRow;
  unsigned int CurCol;
  unsigned int TotalWidth;
  std::string sep;
  std::string firstSep;
  std::string lastSep;
  std::string pm;
};

class TexTableFormater : public TableFormater {
public:
  TexTableFormater(std::ostream *output);
  virtual ~TexTableFormater() {}
  virtual void footer();
  virtual void header();
  virtual void delim();
};

} // ns::ComPWA

#endif
