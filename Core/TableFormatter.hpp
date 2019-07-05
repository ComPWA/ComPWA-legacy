// Copyright (c) 2015 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef CORE_TABLEFORMATER_HPP_
#define CORE_TABLEFORMATER_HPP_

#include <iomanip>
#include <ostream>
#include <string>
#include <vector>

#include "Core/FitParameter.hpp"

namespace ComPWA {

class TableFormatter {
public:
  TableFormatter(std::ostream *output)
      : OutputStream(output), CurRow(0), CurCol(0), TotalWidth(0) {
    sep = "|";
    pm = "+-";
    firstSep = sep;
    lastSep = sep;
  };

  virtual ~TableFormatter() {}

  virtual void reset();

  virtual void delim();

  virtual void footer();

  virtual void header();

  virtual void addColumn(std::string title, unsigned int fixlength = 999);

  void trimString(std::string &src);

  template <typename T> TableFormatter &operator<<(T in) {
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

} // namespace ComPWA

#endif
