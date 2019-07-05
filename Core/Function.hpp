#ifndef CORE_FUNCTION_HPP_
#define CORE_FUNCTION_HPP_

#include <vector>

namespace ComPWA {

/// Interface template for a general Function of the form
/// OutputType (InputTypes)
/// There is no logic for the parameters. They are simply double values that are
/// defined within the Function.
template <typename OutputType, typename... InputTypes> class Function {
public:
  virtual ~Function() = default;
  ///
  virtual OutputType evaluate(const InputTypes &... args) = 0;

  /// It is important to input the vector in the same length and order as
  /// defined in the getParameters() method. So in other words, call
  /// getParameters() first, then modify the contents and finally input them in
  /// this method.
  virtual void updateParametersFrom(const std::vector<double> &) = 0;
  virtual std::vector<double> getParameters() const = 0;
};

/// An Intensity is just a Function that takes a list of data vectors and
/// returns a list of intensities (double)
using Intensity =
    Function<std::vector<double>, std::vector<std::vector<double>>>;

} // namespace ComPWA

#endif
