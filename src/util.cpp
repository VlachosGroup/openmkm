
#include "util.h"
# include <iostream>

using namespace std;

namespace OpenMKM
{

std::vector<double> get_log10_intervals(double end_val, double initial_step, int interval_no)
{
    std::vector<double> steps;
    double step_sz = initial_step * 9/interval_no;
    double r_val = 0;
    //steps.push_back(initial_step);
    while (r_val < end_val) {
        for (size_t i=0; i < interval_no; i++){
            r_val = initial_step + i*step_sz;
            if (r_val > end_val)
                break;
            steps.push_back(r_val);
        }
        if (r_val > end_val){
            steps.push_back(end_val);
            break;
        }
        initial_step *= 10;
        step_sz = initial_step * 9/interval_no;
    }
    return steps;
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
