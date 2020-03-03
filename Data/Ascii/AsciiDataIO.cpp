// This file is part of the ComPWA framework, check
// Copyright (c) 2013 The ComPWA Team.
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <fstream>
#include <sstream>
#include <utility>

#include "Core/Exceptions.hpp"
#include "Core/Logging.hpp"
#include "Data/Ascii/AsciiDataIO.hpp"
#include "Data/Ascii/AsciiHeaderIO.hpp"

namespace ComPWA {
namespace Data {
namespace Ascii {

EventCollection readData(const std::string &InputFilePath,
                         long long NumberEventsToRead) {
  /// -# Open file
  std::ifstream InputStream(InputFilePath);
  if (!InputStream.good())
    throw ComPWA::CorruptFile("Cannot open file " + InputFilePath);

  /// -# Extract header
  AsciiHeader Header;
  Header.importYAML(InputStream);
  auto Position = InputStream.tellg();

  /// -# Determine whether has header
  if (!Header.getFinalStatePIDs().size())
    throw ComPWA::BadConfig("Input data file " + InputFilePath +
                            " misses a header with PIDs");

  /// -# Determine whether weights or not
  bool HasWeights = false;
  double weight, px, py, pz, e;
  std::string Line;
  while (std::getline(InputStream, Line)) {
    if (std::all_of(Line.begin(), Line.end(),
                    [](char c) { return std::isspace(c); }))
      continue;
    std::stringstream StringStream(Line);
    if ((StringStream >> weight) && !(StringStream >> py))
      HasWeights = true;
    break;
  }
  InputStream.seekg(Position);

  /// -# Check number of particles (if it has weights)
  size_t NumberOfParticles = 0;
  if (HasWeights) {
    double weight, px, py, pz, e;
    InputStream >> weight;
    std::string Line;
    while (std::getline(InputStream, Line)) {
      if (std::all_of(Line.begin(), Line.end(),
                      [](char c) { return std::isspace(c); }))
        continue;
      std::stringstream StringStream(Line);
      if (!(StringStream >> e >> px >> py >> pz))
        break;
      ++NumberOfParticles;
    }
    InputStream.seekg(Position);
    if (Header.getFinalStatePIDs().size() != NumberOfParticles)
      throw ComPWA::BadConfig("Number of particles in header is not same as "
                              "number of tuple lines in file \"" +
                              InputFilePath + "\"");
  } else {
    NumberOfParticles = Header.getFinalStatePIDs().size();
  }

  /// -# Check E,p or p,E
  bool IsEnergyFirst = false;
  if (HasWeights)
    InputStream >> weight;
  InputStream >> e >> px >> py >> pz;
  if (e * e - px * px - py * py - pz * pz > 0.)
    IsEnergyFirst = true;
  InputStream.seekg(Position);
  if (Header.isEnergyFirst() != IsEnergyFirst)
    throw ComPWA::BadConfig(
        "Energy-momentum order is not same as stated in header for file \"" +
        InputFilePath + "\"");

  /// -# Import events
  EventCollection ImportedData{Header.getFinalStatePIDs()};

  weight = 1.;
  while (InputStream.good()) {
    if (HasWeights)
      InputStream >> weight;
    std::vector<FourMomentum> FourVectors;
    FourVectors.reserve(NumberOfParticles);
    for (size_t i = 0; i < NumberOfParticles; ++i) {
      if (IsEnergyFirst)
        InputStream >> e >> px >> py >> pz;
      else
        InputStream >> px >> py >> pz >> e;
      FourVectors.push_back(ComPWA::FourMomentum(px, py, pz, e));
    }
    if (!InputStream.fail())
      ImportedData.Events.push_back(Event{FourVectors, weight});
    if (NumberEventsToRead > 0 &&
        ImportedData.Events.size() == (std::size_t)NumberEventsToRead)
      break;
  }
  InputStream.close();

  return ImportedData;
}

void writeData(const EventCollection &DataSample,
               const std::string &OutputFilePath, bool OverwriteFile) {
  if (!DataSample.Events.size())
    throw ComPWA::BadParameter("Cannot write empty event vector");
  ComPWA::Logging log("warning");

  /// -# Determine stream mode
  auto OpenMode = std::ofstream::out;
  if (!OverwriteFile) {
    std::ifstream InputStream(OutputFilePath);
    auto HeaderContent = AsciiHeader::extractHeaderContent(InputStream);
    if (HeaderContent != "") { // if there is a header, check consistency
      AsciiHeader HeaderIn;
      HeaderIn.importYAML(HeaderContent);
      AsciiHeader HeaderOut(DataSample.Pids);
      if (HeaderOut.getFinalStatePIDs() != HeaderIn.getFinalStatePIDs())
        throw ComPWA::BadConfig("Cannot append to file \"" + OutputFilePath +
                                "\", because its PIDs are not the same as the "
                                "event list you try to write");
      if (HeaderOut.getUnit() != HeaderIn.getUnit())
        throw ComPWA::BadConfig("Cannot append to file \"" + OutputFilePath +
                                "\", because there is a mismatch in units");
      if (HeaderOut.isEnergyFirst() != HeaderIn.isEnergyFirst())
        throw ComPWA::BadConfig(
            "Cannot append to file \"" + OutputFilePath +
            "\", because the momentum-energy order is not the same");
    }
    OpenMode |= std::ofstream::app;
  }

  /// -# Open file
  std::ofstream OutputStream;
  OutputStream.open(OutputFilePath, OpenMode);

  /// -# Write header
  if (OverwriteFile) {
    AsciiHeader HeaderOut(DataSample.Pids);
    HeaderOut.dumpToYAML(OutputStream);
    OutputStream << std::endl;
  }

  /// -# Write events
  for (const auto &Event : DataSample.Events) {
    OutputStream << Event.Weight << std::endl;
    for (const auto &FourMom : Event.FourMomenta) {
      OutputStream << FourMom.px() << "\t";
      OutputStream << FourMom.py() << "\t";
      OutputStream << FourMom.pz() << "\t";
      OutputStream << FourMom.e() << std::endl;
    }
  }
}

} // namespace Ascii
} // namespace Data
} // namespace ComPWA