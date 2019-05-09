#include <vector>

#ifndef HTRCT_UTIL_H
#define HTRCT_UTIL_H

namespace HeteroCt
{

std::vector<double> get_log10_intervals(double end_val, double initial_step);

std::vector<double> get_reg_intervals(double start_val, double end_val, double step);

}

#endif
