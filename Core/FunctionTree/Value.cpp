// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Value.hpp"

namespace ComPWA {
namespace FunctionTree {

template class Value<std::complex<double>>;
template class Value<double>;
template class Value<int>;
template class Value<std::vector<std::complex<double>>>;
template class Value<std::vector<double>>;
template class Value<std::vector<int>>;

} // namespace FunctionTree
} // namespace ComPWA
