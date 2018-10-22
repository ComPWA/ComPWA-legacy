// Copyright (c) 2013 Florian Feldbauer.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//_____ I N C L U D E S _______________________________________________________

// ANSI C headers
#include <cstdlib>
#include <fstream>
#include <iostream>
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

// AsciiReader *AsciiReader::clone() const {
// TODO: implement virtual functions and uncomment the following
//	return new AsciiReader(*this);
//  return new AsciiReader();
//}

// AsciiReader *AsciiReader::emptyClone() const { return new AsciiReader(); }

void AsciiReader::writeData(std::shared_ptr<ComPWA::Data::Data> Data,
                            const std::string &OutputFilePath) const {
  LOG(ERROR) << "AsciiReader::writeData() is not implemented!";
}

std::shared_ptr<ComPWA::Data::Data>
AsciiReader::readData(const std::string &InputFilePath) const {

  std::ifstream currentStream;
  currentStream.open(InputFilePath.c_str());

  if (!currentStream)
    throw ComPWA::BadConfig("Can not open " + InputFilePath);

  std::shared_ptr<ComPWA::Data::Data> data(new ComPWA::Data::Data);
  std::vector<ComPWA::Event> Events;

  while (!currentStream.eof()) {
    double e, px, py, pz;
    ComPWA::Event newEvent;

    for (unsigned int ipart = 0; ipart < NumberOfParticles; ++ipart) {
      currentStream >> px >> py >> pz >> e;
      newEvent.addParticle(ComPWA::Particle(px, py, pz, e));
    }

    if (!currentStream.fail()) {
      Events.push_back(newEvent);
      // for ( parts = 0; parts < linesToSkip; parts++ )
      //   currentStream >> px >> py >> pz >> e;
    }
  }
  currentStream.close();
  return data;
}
} // namespace Data
} // namespace ComPWA
