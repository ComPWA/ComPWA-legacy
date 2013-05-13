#include <memory>

#include "Optimizer/ControlParameter.hpp"

std::shared_ptr<ControlParameter> ControlParameter::Instance() {
    return ControlParameter::instance_;
}

std::shared_ptr<ControlParameter> ControlParameter::instance_ = 0;
