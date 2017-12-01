// Copyright (c) 2013 Florian Feldbauer.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//_____ I N C L U D E S _______________________________________________________

// ANSI C headers
#include "AsciiReader.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

#include "Core/Exceptions.hpp"

using namespace ComPWA::DataReader::AsciiReader;

// Constructors and destructors
AsciiReader::AsciiReader(const std::string inConfigFile, const int particles) {

  std::ifstream currentStream;
  currentStream.open(inConfigFile.c_str());

  if (!currentStream)
    throw BadConfig("Can not open " + inConfigFile);

  while (!currentStream.eof()) {
    double e, px, py, pz;
    Event newEvent;

    for (int parts = 0; parts < particles; parts++) {
      currentStream >> px >> py >> pz >> e;
      newEvent.addParticle(Particle(px, py, pz, e));
    }

    if (!currentStream.fail()) {
      fEvents.push_back(newEvent);
      // for ( parts = 0; parts < linesToSkip; parts++ )
      //   currentStream >> px >> py >> pz >> e;
    }
  }
  currentStream.close();
}

AsciiReader::~AsciiReader() { fEvents.clear(); }

AsciiReader *AsciiReader::Clone() const {
  // TODO: implement virtual functions and uncomment the following
  //	return new AsciiReader(*this);
  return new AsciiReader();
}

AsciiReader *AsciiReader::EmptyClone() const { return new AsciiReader(); }

