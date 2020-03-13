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
/// <header>
///   Pids: [211, -211, 22]
///   Order: E px Py pz
///   Unit: GeV
/// </header>
/// ```
/// Note that the section within the `header` tags is YAML syntax.
///
////This header is followed by rows of momentum tuples, grouped per event. In
/// this case, you would have a row for the \f$\pi^+\f$, then for the
/// \f$\pi^-\f$, then one for the \f$\gamma\f$, and finally back to \f$\pi^+\f$.
/// You may choose to start each event group with a weight value, but you don't
/// need to.
EventCollection readData(const std::string &InputFilePath,
                         long long NumberEventsToRead = -1);

void writeData(const EventCollection &DataSample,
               const std::string &OutputFilePath, bool OverwriteFile = true);

} // namespace Ascii
} // namespace Data
} // namespace ComPWA

#endif
