#ifndef CORE_UTILS_HPP_
#define CORE_UTILS_HPP_

#include "Core/Logging.hpp"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/optional.hpp>
#include <boost/property_tree/ptree.hpp>

#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

namespace ComPWA {
namespace Utils {

/// Check of numbers \p x and \p are equal within \p nEpsion times the numerical
/// limit.
inline bool equal(double x, double y, int nEpsilon) {
  return std::abs(x - y) < std::numeric_limits<double>::epsilon() *
                               std::abs(x + y) * nEpsilon ||
         std::abs(x - y) < std::numeric_limits<double>::min();
}

inline double shiftAngle(double value) {
  double originalVal = value;
  while (value > M_PI)
    value -= 2 * M_PI;
  while (value < -1.0 * M_PI)
    value += 2 * M_PI;
  if (value != originalVal)
    LOG(DEBUG) << "shiftAngle() | Shifting parameter from " << originalVal
               << " to " << value << "!";
  return value;
}

/// split the string into pieces, which are separated by the separator
/// character (default separator: space)
inline std::vector<std::string> splitString(const std::string &str,
                                            char separator = ' ') {
  std::vector<std::string> result;
  std::istringstream ss(str);
  std::string token;

  while (std::getline(ss, token, separator)) {
    result.push_back(token);
  }
  return result;
}

} // namespace Utils
} // namespace ComPWA

struct BoolTranslator {
  typedef std::string internal_type;
  typedef bool external_type;

  // Converts a string to bool
  boost::optional<external_type> get_value(const internal_type &str) {
    if (!str.empty()) {
      // first remove leading and trailing whitespace
      auto strcopy(str);
      using boost::algorithm::iequals;
      using boost::algorithm::trim;
      trim(strcopy);
      if (iequals(strcopy, "true") || iequals(strcopy, "yes") || strcopy == "1")
        return boost::optional<external_type>(true);
      else
        return boost::optional<external_type>(false);
    } else
      return boost::optional<external_type>(boost::none);
  }

  // Converts a bool to string
  boost::optional<internal_type> put_value(const external_type &b) {
    return boost::optional<internal_type>(b ? "true" : "false");
  }
};

namespace boost {
namespace property_tree {

template <typename Ch, typename Traits, typename Alloc>
struct translator_between<std::basic_string<Ch, Traits, Alloc>, bool> {
  typedef BoolTranslator type;
};

} // namespace property_tree
} // namespace boost

#endif
