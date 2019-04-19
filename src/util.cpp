
#include "util.h"

namespace HeteroCt
{

std::vector<double> get_times(double end_time)
{
    std::vector<double> times;
    auto scale = 1e-6;
    auto r_time=0;
    while (r_time < end_time) {
        for (size_t i=1; i < 10; i++){
            r_time = i*scale;
            if (r_time > end_time)
                break;
            times.push_back(r_time);
        }
        if (r_time > end_time)
            break;
        scale *= 10;
    }
    return times;
}

}
