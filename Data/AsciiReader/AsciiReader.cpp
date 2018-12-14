// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <fstream>
#include <sstream>
#include <utility>

#include "AsciiReader.hpp"
#include "Core/Exceptions.hpp"

namespace ComPWA {
namespace Data {

// Constructors and destructors
AsciiReader::AsciiReader(unsigned int NumberOfParticles_)
    : NumberOfParticles(NumberOfParticles_) {}

AsciiReader::~AsciiReader() {}

std::shared_ptr<std::vector<ComPWA::Event>>
AsciiReader::readData(const std::string &InputFilePath) const {

  std::ifstream currentStream;
  currentStream.open(InputFilePath);

  if (!currentStream)
    throw ComPWA::BadConfig("Can not open " + InputFilePath);

  std::vector<ComPWA::Event> Events;

  while (!currentStream.eof()) {
    double e, px, py, pz;
    ComPWA::Event newEvent;

    for (unsigned int ipart = 0; ipart < NumberOfParticles; ++ipart) {
      currentStream >> px >> py >> pz >> e;
      newEvent.ParticleList.push_back(ComPWA::Particle(px, py, pz, e));
    }

    if (!currentStream.fail()) {
      Events.push_back(newEvent);
      // for ( parts = 0; parts < linesToSkip; parts++ )
      //   currentStream >> px >> py >> pz >> e;
    }
  }
  currentStream.close();
  return std::make_shared<std::vector<ComPWA::Event>>(Events);
}

} // namespace Data
} // namespace ComPWA
