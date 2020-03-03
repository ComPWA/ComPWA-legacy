// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_DATA_ASCIIHEADER_HPP_
#define COMPWA_DATA_ASCIIHEADER_HPP_

#include <string>
#include <vector>

namespace ComPWA {
namespace Data {
namespace Ascii {

/// Representation of data contained in an ASCII header.
/// This object takes care of reading and writing data concerning a set of
/// events from and to a data file. The concept of a header was introduced to
/// allow the user to document data files and facilitate comparison between
/// [Pawian](https://panda-wiki.gsi.de/foswiki/bin/view/PWA/PawianPwaSoftware)
/// and ComPWA.
class AsciiHeader {
public:
  AsciiHeader(std::vector<int> PIDs = {}, std::string Unit = "GeV",
              bool EnergyFirst = false)
      : PIDs(PIDs), Unit(Unit), EnergyFirst(EnergyFirst) {}

  /// Extract the part that is between the XML/HTML tags `<header>...</header>`
  /// including newlines.
  static std::string extractHeaderContent(std::istream &InputStream);

  /// Set data members by reading a YAML-like string (including newlines).
  /// Example:
  /// ```
  /// Pids: [211, 421, -411]
  /// Unit: GeV
  /// Order: px py pz E
  /// ```
  /// Note that even though the syntax within the header is YAML-like, there is
  /// no full YAML support. In addition, key words are **case-insensitive**.
  void importYAML(const std::string &HeaderContent);
  void importYAML(std::istream &InputStream) {
    importYAML(extractHeaderContent(InputStream));
  }

  /// Serialise data members to YAML format, embedded in XML header tags.
  /// @see extractHeaderContent
  /// @see importYAML
  void dumpToYAML(std::ostream &os) const;

  const std::vector<int> getFinalStatePIDs() const { return PIDs; }
  bool isEnergyFirst() const { return EnergyFirst; }
  const std::string &getUnit() const { return Unit; }

private:
  std::vector<int> PIDs;
  std::string Unit;
  bool EnergyFirst;
};

} // namespace Ascii
} // namespace Data
} // namespace ComPWA

#endif
