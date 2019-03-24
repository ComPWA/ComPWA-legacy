// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_PHYSICS_AMPLITUDE_HPP_
#define COMPWA_PHYSICS_AMPLITUDE_HPP_

#include "Core/FunctionTree.hpp"
#include "Core/ParameterList.hpp"

namespace ComPWA {

struct DataPoint;

namespace Physics {

///
/// \class Amplitude
/// Amplitude interface describes a complex function which resembles
/// the transition of a initial state to a final state via a certain process.
/// The initial and/or final state can be intermediate states, hence an
/// amplitude can describe part of a complete transition. Note that it is not a
/// transition probability and therefore also not normalized. Quantum
/// mechanically speaking the Amplitude is \f$ A = <FS| T |IS> \f$, with the
/// transition operator \f$ T \f$.
///
class Amplitude : public Optimizable, public FunctionTreeInterface {
public:
  virtual ~Amplitude() = default;

  /// calculates the value of the amplitude at the phase space \p point
  virtual std::complex<double> evaluate(const DataPoint &point) const = 0;
};

///
/// \class NamedAmplitude
/// Abstract base class, which adds a Name to an amplitude. The name is
/// currently required to correctly build function trees (because the tree nodes
/// are uniquely identified by the names).
class NamedAmplitude : public Amplitude {
public:
  NamedAmplitude(const std::string &name) : Name(name){};
  virtual ~NamedAmplitude() = default;

  const std::string &getName() const { return Name; }

private:
  std::string Name;
};

} // namespace Physics
} // namespace ComPWA
#endif
