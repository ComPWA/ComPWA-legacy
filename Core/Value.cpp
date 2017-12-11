// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "Core/Value.hpp"

using namespace ComPWA;

template class ComPWA::Value<std::complex<double>>;
template class ComPWA::Value<double>;
template class ComPWA::Value<int>;
template class ComPWA::Value<std::vector<std::complex<double>>>;
template class ComPWA::Value<std::vector<double>>;
template class ComPWA::Value<std::vector<int>>;
