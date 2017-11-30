// Copyright (c) 2015 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/TableFormater.hpp"

namespace ComPWA {
void TableFormater::delim() {
  *out << sep;
  for (unsigned int i = 0; i < totalWidth - 1; i++)
    *out << "-";
  *out << sep << std::endl;
}

void TableFormater::footer() { delim(); }

void TableFormater::header() {
  delim();
  for (unsigned int i = 0; i < columnTitle.size(); i++)
    *out << "| " << std::setw(columnWidth[i]) << columnTitle[i] << " ";
  *out << "|" << std::endl;
  delim();
}

void TableFormater::addColumn(std::string title, unsigned int fixlength) {
  unsigned int length;
  if (fixlength != 999)
    length = fixlength;
  else
    length = title.length();
  columnWidth.push_back(length);
  totalWidth += length + 3;
  columnTitle.push_back(title);
}

void TableFormater::Reset() {
  curCol = 0;
  curRow = 0;
  totalWidth = 0;
  columnWidth.clear();
  columnTitle.clear();
}

void TableFormater::trimString(std::string &src) {
  // remove  all final zeros
  char chr = '0';
  std::string::size_type pos2 = src.find_last_not_of(chr, src.length());
  if (pos2 + 1 < src.length())
    src.erase(pos2 + 1, src.length());
}

TableFormater &TableFormater::operator<<(DoubleParameter in) {
  if (curCol == 0)
    *out << firstSep + " ";
  else
    *out << " " << sep << " ";
  if (in.hasError()) {
    std::string tmp;
    if (in.errorType() == ErrorType::SYM && in.error().first != 0) {
      unsigned int halfWidth =
          (unsigned int)(columnWidth[curCol]) / 2; // divide column width
      *out << std::setw(halfWidth) << in.value();
      tmp = pm + std::to_string((long double)in.error().first);
      trimString(tmp);
      *out << std::setw(halfWidth) << tmp;
    } else if (in.errorType() == ErrorType::ASYM) {
      unsigned int w =
          (unsigned int)(columnWidth[curCol]) / 3; // take 1/3 of column width
      tmp = std::to_string((long double)in.value());
      trimString(tmp);
      *out << std::setw(w) << tmp;
      tmp = "-" + std::to_string((long double)in.error().first);
      trimString(tmp);
      *out << std::setw(w) << tmp;
      tmp = "+" + std::to_string((long double)in.error().second);
      trimString(tmp);
      *out << std::setw(w) << tmp;
    } else
      *out << std::setw(columnWidth[curCol]) << in.value();
  } else {
    *out << std::setw(columnWidth[curCol]) << in.value();
  }
  curCol++;
  if (curCol == columnWidth.size()) {
    *out << " " << lastSep << std::endl;
    curRow++;
    curCol = 0;
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
  *out << "\\end{tabular}" << std::endl;
  //	*out <<"\label{...}"<<endl;
  //	*out <<"\caption{...}"<<endl;
};
void TexTableFormater::header() {
  //	*out <<"\begin{table}"<<endl;
  //	*out <<"\centering"<<endl;
  *out << "\\begin{tabular}{|";
  for (int i = 0; i < columnTitle.size(); i++)
    *out << "c|";
  *out << "}" << std::endl;
  delim();
  *out << firstSep;
  for (unsigned int i = 0; i < columnTitle.size(); i++) {
    *out << std::setw(columnWidth[i]) << columnTitle[i];
    if (i == columnTitle.size() - 1)
      *out << lastSep;
    else
      *out << sep;
  }
  *out << std::endl;
  delim();
  delim();
};
void TexTableFormater::delim() { *out << "\\hline" << std::endl; }

} /* namespace ComPWA */
