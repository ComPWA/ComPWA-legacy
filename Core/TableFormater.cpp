// Copyright (c) 2015 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/TableFormater.hpp"

using namespace ComPWA;

void TableFormater::delim() {
  *OutputStream << sep;
  for (unsigned int i = 0; i < TotalWidth - 1; i++)
    *OutputStream << "-";
  *OutputStream << sep << std::endl;
}

void TableFormater::footer() { delim(); }

void TableFormater::header() {
  delim();
  for (unsigned int i = 0; i < ColumnTitle.size(); i++)
    *OutputStream << "| " << std::setw(ColumnWidth.at(i)) << ColumnTitle.at(i)
                  << " ";
  *OutputStream << "|" << std::endl;
  delim();
}

void TableFormater::addColumn(std::string title, unsigned int fixlength) {
  unsigned int length;
  if (fixlength != 999)
    length = fixlength;
  else
    length = title.length();
  ColumnWidth.push_back(length);
  TotalWidth += length + 3;
  ColumnTitle.push_back(title);
}

void TableFormater::reset() {
  CurCol = 0;
  CurRow = 0;
  TotalWidth = 0;
  ColumnWidth.clear();
  ColumnTitle.clear();
}

void TableFormater::trimString(std::string &src) {
  // remove  all final zeros
  char chr = '0';
  std::string::size_type pos2 = src.find_last_not_of(chr, src.length());
  if (pos2 + 1 < src.length())
    src.erase(pos2 + 1, src.length());
}

TableFormater &TableFormater::operator<<(FitParameter in) {
  if (CurCol == 0)
    *OutputStream << firstSep + " ";
  else
    *OutputStream << " " << sep << " ";
  if (in.hasError()) {
    std::string tmp;
    if (in.errorType() == ErrorType::SYM && in.avgError() != 0) {
      unsigned int halfWidth =
          (unsigned int)(ColumnWidth.at(CurCol)) / 2; // divide column width
      *OutputStream << std::setw(halfWidth) << in.value();
      tmp = pm + std::to_string((long double)in.avgError());
      trimString(tmp);
      *OutputStream << std::setw(halfWidth) << tmp;
    } else if (in.errorType() == ErrorType::ASYM) {
      unsigned int w =
          (unsigned int)(ColumnWidth.at(CurCol)) / 3; // take 1/3 of column width
      tmp = std::to_string((long double)in.value());
      trimString(tmp);
      *OutputStream << std::setw(w) << tmp;
      tmp = "-" + std::to_string((long double)in.error().first);
      trimString(tmp);
      *OutputStream << std::setw(w) << tmp;
      tmp = "+" + std::to_string((long double)in.error().second);
      trimString(tmp);
      *OutputStream << std::setw(w) << tmp;
    } else
      *OutputStream << std::setw(ColumnWidth.at(CurCol)) << in.value();
  } else {
    *OutputStream << std::setw(ColumnWidth.at(CurCol)) << in.value();
  }
  CurCol++;
  if (CurCol == ColumnWidth.size()) {
    *OutputStream << " " << lastSep << std::endl;
    CurRow++;
    CurCol = 0;
  }
  return *this;
}

//====== TexTableFormater =========
TexTableFormater::TexTableFormater(std::ostream *output)
    : TableFormater(output) {
  sep = "&";
  pm = "$\\pm$";
  firstSep = "";
  lastSep = "\\\\";
}

void TexTableFormater::footer() {
  delim();
  delim();
  *OutputStream << "\\end{tabular}" << std::endl;
  //  *out <<"\label{...}"<<endl;
  //  *out <<"\caption{...}"<<endl;
};

void TexTableFormater::header() {
  //  *out <<"\begin{table}"<<endl;
  //  *out <<"\centering"<<endl;
  *OutputStream << "\\begin{tabular}{|";
  for (unsigned int i = 0; i < ColumnTitle.size(); ++i)
    *OutputStream << "c|";
  *OutputStream << "}" << std::endl;
  delim();
  *OutputStream << firstSep;
  for (unsigned int i = 0; i < ColumnTitle.size(); ++i) {
    *OutputStream << std::setw(ColumnWidth.at(i)) << ColumnTitle.at(i);
    if (i == ColumnTitle.size() - 1)
      *OutputStream << lastSep;
    else
      *OutputStream << sep;
  }
  *OutputStream << std::endl;
  delim();
  delim();
};

void TexTableFormater::delim() { *OutputStream << "\\hline" << std::endl; }

