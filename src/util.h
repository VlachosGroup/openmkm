//! @file util.h

// This file is part of OpenMKM. See License.txt in the top-level directory 
// for license and copyright information.

#ifndef OMKM_UTIL_H
#define OMKM_UTIL_H

#include <vector>

namespace OpenMKM
{

std::vector<double> get_log10_intervals(double end_val, 
                                        double initial_step, 
                                        int intervals=20);

std::vector<double> get_reg_intervals(double start_val, 
                                      double end_val, 
                                      double step);

}

#endif
