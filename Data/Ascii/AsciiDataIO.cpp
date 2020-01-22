// This file is part of the ComPWA framework, check
// Copyright (c) 2013 The ComPWA Team.
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <fstream>
#include <sstream>
#include <utility>

#include "Data/Ascii/AsciiDataIO.hpp"
#include "Core/Exceptions.hpp"
#include "Core/Logging.hpp"

namespace ComPWA {
namespace Data {
namespace Ascii {

std::string tolower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return s;
}

bool isHeaderLine(std::string line) {
  return (tolower(line).find("header") != std::string::npos);
}

bool isEmptyString(const std::string &line) {
  return std::all_of(line.begin(), line.end(),
                     [](char c) { return std::isspace(c); });
}

std::ifstream openStream(const std::string &Filename) {
  std::ifstream InputStream(Filename);
  if (!InputStream.good())
    throw ComPWA::CorruptFile("Cannot open file " + Filename);
  return InputStream;
}

std::vector<int> extractHeader(std::ifstream &InputStream) {
  std::string line;
  // Find first header keyword
  while (std::getline(InputStream, line))
    if (isHeaderLine(line))
      break;
  // Import PIDs
  std::vector<int> PIDs;
  while (std::getline(InputStream, line)) {
    if (isHeaderLine(line))
      break;
    std::stringstream StringStream(line);
    std::string Word;
    int PID;
    if (StringStream >> Word >> PID) {
      Word = tolower(Word);
      if (Word.find("pid") != std::string::npos)
        PIDs.push_back(PID);
    }
  }
  return PIDs;
}

std::vector<int> extractHeader(const std::string &Filename) {
  auto InputStream = openStream(Filename);
  return extractHeader(InputStream);
}

std::vector<ComPWA::Event> readData(const std::string &InputFilePath,
                                    long long NumberEventsToRead) {
  /// -# Open file
  auto InputStream = openStream(InputFilePath);

  /// -# Extract header
  auto PIDs = extractHeader(InputStream);
  auto NumberOfParticles = PIDs.size();
  if (!PIDs.size())
    throw ComPWA::BadConfig("Input data file " + InputFilePath +
                            " misses a header with PIDs");

  /// -# Determine whether weights or not
  auto Position = InputStream.tellg();
  bool HasWeights = false;
  double weight, px, py, pz, e;
  std::string line;
  while (std::getline(InputStream, line)) {
    if (isEmptyString(line))
      continue;
    std::stringstream StringStream(line);
    if ((StringStream >> weight) && !(StringStream >> py))
      HasWeights = true;
    break;
  }
  InputStream.seekg(Position);

  /// -# Determine E,p or p,E
  bool IsOrderEnergyMomentum = false;
  if (HasWeights)
    InputStream >> weight;
  InputStream >> e >> px >> py >> pz;
  if (e * e - px * px - py * py - pz * pz > 0.)
    IsOrderEnergyMomentum = true;
  InputStream.seekg(Position);

  /// -# Import events
  std::vector<ComPWA::Event> Events;
  weight = 1.;
  while (InputStream.good()) {
    if (HasWeights)
      InputStream >> weight;
    std::vector<Particle> ParticleList;
    ParticleList.reserve(NumberOfParticles);
    for (size_t i = 0; i < NumberOfParticles; ++i) {
      if (IsOrderEnergyMomentum)
        InputStream >> e >> px >> py >> pz;
      else
        InputStream >> px >> py >> pz >> e;
      ParticleList.push_back(ComPWA::Particle(px, py, pz, e, PIDs[i]));
    }
    if (!InputStream.fail())
      Events.push_back({ParticleList, weight});
    if (NumberEventsToRead > 0 &&
        Events.size() == (std::size_t)NumberEventsToRead)
      break;
  }
  InputStream.close();
  return Events;
}

void writeData(const std::vector<ComPWA::Event> &Events,
               const std::string &OutputFilePath, bool AppendToFile) {
  if (!Events.size())
    throw ComPWA::BadParameter("Cannot write empty event vector");
  ComPWA::Logging log("warning");

  /// -# Determine stream mode
  auto OpenMode = std::ofstream::out;
  if (AppendToFile) {
    auto PIDs = extractHeader(OutputFilePath);
    if (!PIDs.size()) {
      LOG(WARNING) << "Will overwrite corrupt output file " << OutputFilePath
                   << std::endl;
      AppendToFile = false;
    } else
      OpenMode |= std::ofstream::app;
  }

  /// -# Open file
  std::ofstream OutputStream;
  OutputStream.open(OutputFilePath, OpenMode);
  if (!OutputStream)
    throw ComPWA::BadConfig("Cannot open " + OutputFilePath);

  /// -# Write header
  if (!AppendToFile) {
    OutputStream << "Header" << std::endl;
    for (auto Particle : Events[0].ParticleList)
      OutputStream << "\tPid: " << Particle.pid() << std::endl;
    OutputStream << "Header" << std::endl << std::endl;
  }

  /// -# Write events
  for (const auto &Event : Events) {
    OutputStream << Event.Weight << std::endl;
    for (const auto &Particle : Event.ParticleList) {
      OutputStream << Particle.fourMomentum().px() << "\t";
      OutputStream << Particle.fourMomentum().py() << "\t";
      OutputStream << Particle.fourMomentum().pz() << "\t";
      OutputStream << Particle.fourMomentum().e() << std::endl;
    }
  }
} // namespace Ascii

} // namespace Ascii
} // namespace Data
} // namespace ComPWA