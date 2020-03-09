// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_DATA_ROOTDATAIO_HPP_
#define COMPWA_DATA_ROOTDATAIO_HPP_

#include <memory>
#include <string>

#include "Core/Event.hpp"

class TTree;

namespace ComPWA {
namespace Data {
/// Namespace with read/write functions for [ROOT
/// files](https://root.cern.ch/input-and-output).
namespace Root {

/// Create a vector of `Event`s from a [ROOT
/// file](https://root.cern.ch/input-and-output). The input file should have at
/// least one event-based
/// [`TTree`](https://root.cern.ch/doc/master/classTTree.html) with:
/// * A branch named `"Particles"` containing
/// [`TClonesArray`](https://root.cern.ch/doc/master/classTClonesArray.html)s.
/// These arrays should contain
/// [`TParticle`](https://root.cern.ch/doc/master/classTParticle.html) objects
/// with a defined 4-momentum.
/// * A branch of `double`s called `"Weight"`.
/// \param InputFileName Input ROOT file(s); can take wildcards, see
/// [`TChain::Add`](https://root.cern.ch/doc/master/classTChain.html).
/// \param TreeName Name of the event-based tree
/// \param NumberEventsToRead Limit the resulting vector to this number of
/// events (optional).
ComPWA::EventCollection readData(const std::string &InputFileName,
                                 const std::string &TreeName,
                                 long long NumberEventsToRead = -1);

/// Write a vector of `Event`s to a [ROOT
/// file](https://root.cern.ch/input-and-output). See `readData` for the
/// structure of the output file.
/// \param OutputSample List of `Event`s including a info header
/// \param OutputFilePath Path to the output ROOT file
/// \param TreeName Name of the event-based output `TTree` in the file
/// \param OverwriteFile Set to `true` if you do *not* want to append
void writeData(const EventCollection &OutputSample,
               const std::string &OutputFilePath, const std::string &TreeName,
               bool OverwriteFile = true);

} // namespace Root
} // namespace Data
} // namespace ComPWA

#endif
