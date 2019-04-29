
#include "util.h"
# include <iostream>

using namespace std;

namespace HeteroCt
{

std::vector<double> get_log_intervals(double end_val)
{
    std::vector<double> intervals;
    auto scale = 1e-6;
    double r_val = 0;
    while (r_val < end_val) {
        for (size_t i=1; i < 10; i++){
            r_val = i*scale;
            if (r_val > end_val)
                break;
            intervals.push_back(r_val);
        }
        if (r_val > end_val)
            break;
        scale *= 10;
    }
    return intervals;
}

}
