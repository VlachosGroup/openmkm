//! @file hctexceptions.cpp

// This file is part of HeteroCt. See License.txt in the top-level directory. 

#include <sstream>

#include "cantera/base/ctexceptions.h"
//#include "application.h"
#include "cantera/base/global.h"

#include "hctexceptions.h"

namespace HeteroCt
{

// *** Exceptions ***
std::string YAMLParserError::getMessage() const
{
    return fmt::format("Failed to parse YAML node ({}). {}.",
                       node_, failure_mode_);
}

} // namespace HeteroCt
