// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_DATA_ASCIIDATAIO_HPP_
#define COMPWA_DATA_ASCIIDATAIO_HPP_

#include <memory>
#include <vector>

#include "Core/Event.hpp"

namespace ComPWA {
namespace Data {
namespace Ascii {

/// Read momentum tuples from an ASCII file.
/// The file should start with a header that defines the final state, like so:
/// ```
/// Header
///   Pid: 211
///   Pid: -211
///   Pid: 22
/// Header
/// ```
/// This header is followed by rows of momentum tuples, grouped per event. In
/// this case, you would have a row for the \f$\pi^+\f$, then for the
/// \f$\pi^-\f$, then one for the \f$\gamma\f$, and finally back to \f$\pi^+\f$.
/// You may choose to start each event group with a weight value, but you don't
/// need to.
EventList readData(const std::string &InputFilePath,
                   long long NumberEventsToRead = -1);

void writeData(const EventList &EvtList, const std::string &OutputFilePath,
               bool AppendToFile = false);

} // namespace Ascii
} // namespace Data
} // namespace ComPWA

#endif
