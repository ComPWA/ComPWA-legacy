#include <iostream>
#include <sstream>
#include <vector>

namespace ComPWA {
namespace Utils {

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
