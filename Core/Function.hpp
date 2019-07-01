#ifndef CORE_FUNCTION_HPP_
#define CORE_FUNCTION_HPP_

#include <vector>

namespace ComPWA {

template <typename OutputType, typename... InputTypes> class Function {
public:
  virtual ~Function() = default;
  // we dont make the evaluate function const anymore, because
  // it will allow internal modification like caching
  virtual OutputType evaluate(const InputTypes &... args) = 0;
  // changes parameters to the given values in the list
  // The order of the parameters in the list is important.
  // It has to be the same as returned by getParameters()
  virtual void updateParametersFrom(const std::vector<double> &) = 0;
  // gets a list of parameters defined by this function
  virtual std::vector<double> getParameters() const = 0;
};

// and intensity is just a function which takes a list of data vectors and
// returns a list of intensities (double)
using Intensity = Function<double, std::vector<double>>;

} // namespace ComPWA

#endif
