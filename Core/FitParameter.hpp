#ifndef CORE_FITPARAMETER_HPP_
#define CORE_FITPARAMETER_HPP_

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
  bool IsFixed = false;
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
