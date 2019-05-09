
#include "util.h"
# include <iostream>

using namespace std;

namespace HeteroCt
{

std::vector<double> get_log10_intervals(double end_val, double initial_step)
{
    std::vector<double> intervals;
    auto scale = initial_step;
    double r_val = 0;
    while (r_val < end_val) {
        for (size_t i=1; i < 10; i++){
            r_val = i*scale;
            if (r_val > end_val)
                break;
            intervals.push_back(r_val);
        }
        if (r_val > end_val){
            intervals.push_back(end_val);
            break;
        }
        scale *= 10;
    }
    return intervals;
}

std::vector<double> get_reg_intervals(double start_val, double end_val, double step)
{
    std::vector<double> intervals;
    double r_val = start_val;
    while (r_val < end_val) {
        r_val += step;
        if (r_val > end_val){
            intervals.push_back(end_val);
            break;
        }
        intervals.push_back(r_val);
    }
    return intervals;
}

}
