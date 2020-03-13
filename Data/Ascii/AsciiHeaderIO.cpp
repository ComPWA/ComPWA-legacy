// This file is part of the ComPWA framework, check
// Copyright (c) 2013 The ComPWA Team.
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <algorithm>
#include <fstream>
#include <regex>
#include <sstream>

#include "Core/Exceptions.hpp"
#include "Data/Ascii/AsciiHeaderIO.hpp"

namespace ComPWA {
namespace Data {
namespace Ascii {

/// @cond INTERNAL

namespace std_fix {
std::string tolower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}
} // namespace std_fix

/// Check if a line starts with a float/int value
bool isValueLine(const std::string &line) {
  std::stringstream ss(line);
  float value;
  return (bool)(ss >> value);
}

/// @endcond

std::string AsciiHeader::extractHeaderContent(std::istream &InputStream) {
  /// -# Define regular expressions
  static const std::regex HeaderBegin("<header>", std::regex_constants::icase);
  static const std::regex HeaderEnd("</header>", std::regex_constants::icase);
  /// -# Define buffers
  std::string Line;
  std::stringstream HeaderContent;
  std::smatch RegexMatch;
  /// -# Find header beginning
  auto pos = InputStream.tellg();
  while (std::getline(InputStream, Line)) {
    if (isValueLine(Line)) {
      InputStream.seekg(pos);
      return ""; // means there is no header
    }
    if (std::regex_search(Line, RegexMatch, HeaderBegin)) {
      HeaderContent << RegexMatch.suffix() << std::endl;
      break;
    }
  }
  /// -# Escape if header content contains closing tag already
  Line = HeaderContent.str();
  if (std::regex_search(Line, RegexMatch, HeaderEnd)) {
    Line = RegexMatch.prefix();
    return Line;
  }
  /// -# Abort if end of file is reached already
  if (InputStream.eof()) {
    InputStream.clear();
    InputStream.seekg(pos);
    throw ComPWA::CorruptFile("No opening <header> tag");
  }
  /// -# Import header content
  while (std::getline(InputStream, Line)) {
    if (isValueLine(Line))
      break;
    if (std::regex_search(Line, RegexMatch, HeaderEnd)) {
      HeaderContent << RegexMatch.prefix() << std::endl;
      return HeaderContent.str();
    }
    HeaderContent << Line << std::endl;
  }
  InputStream.clear();
  InputStream.seekg(pos);
  throw ComPWA::CorruptFile("No closing tag after <header>");
}

/// Set data members by reading a YAML-like string (including newline
/// characters).
void AsciiHeader::importYAML(const std::string &HeaderContent) {
  std::stringstream StringStream(HeaderContent);
  std::string Line;
  static const std::regex RegexKeyValue(
      R"(^\s*([^\s]+.*?)\s*:\s+([^\s]+.*?)\s*$)");
  static const std::regex RegexPID(R"(-?\d+)");
  static const std::regex RegexEfirst("e.*[mpz]", std::regex_constants::icase);
  static const std::regex RegexPfirst("[mp].*e", std::regex_constants::icase);
  while (getline(StringStream, Line)) {
    std::smatch RegexMatch;
    std::regex_search(Line, RegexMatch, RegexKeyValue);
    if (RegexMatch.size() < 3)
      continue;
    const std::string Key{std_fix::tolower(RegexMatch[1])};
    const std::string Value{RegexMatch[2]};
    if (Key == "pids") {
      std::sregex_iterator iter(Value.begin(), Value.end(), RegexPID), end;
      std::vector<int>().swap(PIDs); // completely clear vector
      for (; iter != end; ++iter)
        PIDs.push_back(std::stoi(iter->str()));
      continue;
    } else if (Key == "unit") {
      Unit = Value;
    } else if (Key == "order") {
      if (std::regex_match(Value, RegexPfirst)) {
        EnergyFirst = false;
      } else if (std::regex_match(Value, RegexEfirst)) {
        EnergyFirst = true;
      }
    }
  }
}

void AsciiHeader::dumpToYAML(std::ostream &os) const {
  os << "<header>" << std::endl;
  /// -# Write PIDs if available
  if (PIDs.size()) {
    os << "\tPids: [" << PIDs[0];
    for (size_t i = 1; i < PIDs.size(); ++i)
      os << ", " << PIDs[i];
    os << "]" << std::endl;
  }
  /// -# Write order
  os << "\tOrder: ";
  if (EnergyFirst)
    os << "E px py pz" << std::endl;
  else
    os << "px py pz E" << std::endl;
  /// -# Write unit
  os << "\tUnit: " << Unit << std::endl;
  os << "</header>" << std::endl;
}

} // namespace Ascii
} // namespace Data
} // namespace ComPWA