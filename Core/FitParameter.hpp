#ifndef CORE_FITPARAMETER_HPP_
#define CORE_FITPARAMETER_HPP_

#include "Core/Function.hpp"
#include "Core/Logging.hpp"

#include <ostream>
#include <string>
#include <vector>

#include <boost/serialization/nvp.hpp>

namespace ComPWA {

template <typename T> struct FitParameter {
  FitParameter() = default;
  FitParameter(std::string name, T val, bool isfixed = true)
      : HasBounds(false), IsFixed(isfixed), Value(val), Name(name),
        Bounds(0.0, 0.0) {}
  FitParameter(std::string name, T val, T min, T max, bool isfixed = true)
      : HasBounds(true), IsFixed(isfixed), Value(val), Name(name),
        Bounds(min, max) {}
  bool HasBounds = false;
  bool IsFixed = true;
  T Value;
  std::pair<T, T> Error;
  std::string Name = "";
  std::pair<T, T> Bounds;

  friend std::ostream &operator<<(std::ostream &os,
                                  const FitParameter<double> &x) {
    os << x.Name << ": " << x.Value;
    if (x.Error.first != 0.0 || x.Error.second != 0.0) {
      if (x.Error.first == x.Error.second) {
        os << " +- " << x.Error.second;
      } else {
        os << " + " << x.Error.second << " - " << x.Error.first;
      }
      if (x.HasBounds) {
        os << " Bounds: [" << x.Bounds.first << ", " << x.Bounds.second << "]";
      }
    }
    if (x.IsFixed)
      os << " (fixed)";
    return os;
  }
};

using FitParameterList = std::vector<FitParameter<double>>;

inline bool isValid(const FitParameterList &FitParameters,
                    const std::vector<ComPWA::Parameter> &EstimatorParameters) {
  // validate FitParameterList
  auto ActualParametersIt = EstimatorParameters.begin();
  for (auto const &Par : FitParameters) {
    if (ActualParametersIt == EstimatorParameters.end()) {
      LOG(ERROR) << "Less fit parameters given then actual parameters in the "
                    "function!";
      return false;
    }

    if (Par.Name != ActualParametersIt->Name) {
      LOG(ERROR) << "Wrong ordering of fit parameters!";
      return false;
    }
    ++ActualParametersIt;
  }
  return true;
}

} // namespace ComPWA

namespace boost {
namespace serialization {

template <class Archive>
void serialize(Archive &ar, ComPWA::FitParameter<double> &Par,
               const unsigned int version) {
  ar &BOOST_SERIALIZATION_NVP(Par.Name);
  ar &BOOST_SERIALIZATION_NVP(Par.Value);
  ar &BOOST_SERIALIZATION_NVP(Par.Error);
  ar &BOOST_SERIALIZATION_NVP(Par.IsFixed);
  ar &BOOST_SERIALIZATION_NVP(Par.HasBounds);
  ar &BOOST_SERIALIZATION_NVP(Par.Bounds);
}

} // namespace serialization
} // namespace boost

#endif
