#include "headers.h"
#include "utility_functions.h"
//Pending recent


//Pending long term
//Add function to calculate overlap, scattering and change of the upconversion index
//Current sim does not accept 1480 pump
    
int main(int argc, char **argv)
{
    Simulation::Init_params a;

    
    Simulation::Result test;
    
    
    test.simulate_debug(1.0e20, 5.0e20, 1000);
    test.plot_steps("output_file.txt");
    
    /*
    test.advance_step(0, 1.0e20, 9.0e20, 1000);
    test.report_step(0);
    test.regress_step(1, 1.0e20, 9.0e20, 1000);
    test.report_step(1);
    test.report_step(0);
    */
    
    return 0;
}
