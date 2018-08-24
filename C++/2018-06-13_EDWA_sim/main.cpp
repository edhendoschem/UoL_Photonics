#include "headers.h"
#include "utility_functions.h"
//Pending recent


//Pending long term
//Add function to calculate overlap scattering
//Current sim does not accept 1480 pump, incorporate that posibility
    
int main(int argc, char **argv)
{
    Simulation::Init_params a;
    Simulation::Result test;
    //std::cout<<"Ccr = "<<a.Ccr<<'\n';
    //std::cout<<"Cup = "<<a.Cup<<'\n';

    /*
    test.find_all_n(0, 1.0e20, 5.0e20, 1000);
    
    
    std::cout<<"1535, 21 = "<<Utility::return_cross_section(1535000, Cross_sec::Er_emi)<<'\n';
    std::cout<<"1535, 12 = "<<Utility::return_cross_section(1535000, Cross_sec::Er_abs)<<'\n';
    std::cout<<"dps_ind = "<<test.dPs_ind(0, 1533000, 1.0e-6)<<'\n';
    test.report_step(0, false);
    test.advance_signal(0, 1);
    test.report_step(1, false);
    */
    //test.simulate_debug();
    //test.save_data("output_file.txt", true);
    test.plot_data("output_file.txt", "plot_s_test.txt", true);
    
    return 0;
}
