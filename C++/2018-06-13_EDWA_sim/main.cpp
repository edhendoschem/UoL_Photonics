#include "headers.h"
#include "utility_functions.h"
//Pending recent


//Pending long term
//Add function to calculate overlap scattering
//Current sim does not accept 1480 pump, incorporate that posibility
    
int main(int argc, char **argv)
{
    Simulation::Init_params a{};
    
    a.l = 0.05;
    /*
    a.NEr /= 2.0;
    a.NYb /= 2.0;
    a.recalculate_constants();
    */
    
    a.NEr = 4.45e26;
    a.NYb = 4.45e26;
    a.recalculate_constants();
    Simulation::Result r{a};
    r.simulate(false);
    r.save_data("output.txt", true);
    r.plot_data("output.txt","plot_script.txt");
    
    //std::vector<std::array<double, 4>> results = find_ratio(a);
    
    
 
    return 0;
}
