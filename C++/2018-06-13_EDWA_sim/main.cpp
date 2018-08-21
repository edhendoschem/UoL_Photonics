#include "headers.h"
#include "utility_functions.h"
//Pending recent
//Add Runge Kutta to Pp_f and Ps
//Refactor wl_map to use int instead of double

//Pending long term
//Add function to calculate overlap, scattering and change of the upconversion index
//Current sim does not accept 1480 pump
    
int main(int argc, char **argv)
{
    Simulation::Init_params a;

    
    Simulation::Result test;
    
    
    test.simulate_debug(1.0e20, 9.0e20, 1000);
    //test.report_step(0);
    //test.report_step(test.data.size()/4);
    //test.report_step(test.data.size()*2/4);
    //test.report_step(test.data.size()*3/4);
    //test.report_step(test.data.size()-1);
    
    /*
    test.advance_step(0, 1.0e20, 9.0e20, 1000);
    test.report_step(0);
    test.regress_step(1, 1.0e20, 9.0e20, 1000);
    test.report_step(1);
    test.report_step(0);
    */
    
    /*
    auto f1 = [] (std::vector<double> X)
    {
        return 3.0*X[0] - cos(X[1]*X[2]) - 0.5;
    };    
    
    auto f2 = [] (std::vector<double> X)
    {
        return X[0] * X[0] - 81.0 * (X[1] + 0.1) * (X[1] + 0.1) + sin(X[2]) + 1.06;
    };    
    
    auto f3 = [] (std::vector<double> X)
    {
        return exp(-X[0]*X[1]) + 20.0 * X[2] + (10.0 * PI - 3) / 3.0;
    };
    
    std::vector<double (*)(std::vector<double>)> fs;
    fs.push_back(f1);
    fs.push_back(f2);
    fs.push_back(f3);
    std::vector<double> x0 {0.1, 0.1, -0.1};
    */
    
    return 0;
}
