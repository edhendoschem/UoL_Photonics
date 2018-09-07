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
    /*
    a.lsEr = 0;
    a.lsYb = 0;
    a.lpEr = 0;
    a.lpYb = 0;
    */

    std::cout<<a.lsEr<<", "<<a.lpEr<<", "<<a.lsYb<<", "<<a.lpYb<<'\n';
    
    a.NEr = 4.45e26;
    a.NYb = 4.45e26;
    a.recalculate_constants();
    Simulation::Result r{a};
    r.simulate(false); //true = warn of negative pop values
    r.save_spectral_data("spectral.txt",r.data.size()-1, true);
    //r.save_data("output.csv", false); //true = signal, ase powers in dBm, else mW
    //r.plot_data("output.csv","plot_script.txt"); //
    
    //std::vector<std::array<double, 4>> results = find_ratio(a);
    
    
 
    return 0;
}
