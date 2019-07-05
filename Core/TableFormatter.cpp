// Copyright (c) 2015 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "TableFormatter.hpp"

namespace ComPWA {

void TableFormatter::delim() {
  *OutputStream << sep;
  for (unsigned int i = 0; i < TotalWidth - 1; i++)
    *OutputStream << "-";
  *OutputStream << sep << std::endl;
}

void TableFormatter::footer() { delim(); }

void TableFormatter::header() {
  delim();
  for (unsigned int i = 0; i < ColumnTitle.size(); i++)
    *OutputStream << "| " << std::setw(ColumnWidth.at(i)) << ColumnTitle.at(i)
                  << " ";
  *OutputStream << "|" << std::endl;
  delim();
}

void TableFormatter::addColumn(std::string title, unsigned int fixlength) {
  unsigned int length;
  if (fixlength != 999)
    length = fixlength;
  else
    length = title.length();
  ColumnWidth.push_back(length);
  TotalWidth += length + 3;
  ColumnTitle.push_back(title);
}

void TableFormatter::reset() {
  CurCol = 0;
  CurRow = 0;
  TotalWidth = 0;
  ColumnWidth.clear();
  ColumnTitle.clear();
}

void TableFormatter::trimString(std::string &src) {
  // remove  all final zeros
  char chr = '0';
  std::string::size_type pos2 = src.find_last_not_of(chr, src.length());
  if (pos2 + 1 < src.length())
    src.erase(pos2 + 1, src.length());
}

} // namespace ComPWA
