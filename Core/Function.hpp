#ifndef CORE_FUNCTION_HPP_
#define CORE_FUNCTION_HPP_

#include <string>
#include <unordered_map>
#include <vector>

namespace ComPWA {

struct Parameter {
  std::string Name;
  double Value;
};

using DataMap = std::unordered_map<std::string, std::vector<double>>;

/// Interface template for a general Function of the form
/// OutputType Function(InputTypes)
/// The concept closely follows the mathematical definition of a
/// function/mapping. The parameters are stated by the Function and can be
/// retrieved via getParameters(). The only difference to a mathematical
/// function is that the evaluation and the setting of the parameters are
/// separated. Parameter have to be altered with updateParametersFrom().
template <typename OutputType, typename... InputTypes> class Function {
public:
  virtual ~Function() = default;

  virtual OutputType evaluate(const InputTypes &... args) noexcept = 0;

  /// It is important to input the vector in the same length and order as
  /// defined in the getParameters() method. So in other words, call
  /// getParameters() first, then use this ordering and to input new values in
  /// this method.
  virtual void updateParametersFrom(const std::vector<double> &) = 0;
  virtual std::vector<Parameter> getParameters() const = 0;
};

/// An Intensity is just a Function that takes a list of data vectors and
/// returns a list of intensities (double)
using Intensity = Function<std::vector<double>, DataMap>;

} // namespace ComPWA

#endif
