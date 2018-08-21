#ifndef UTILITY_FUNCTIONS_H_INCLUDED
#define UTILITY_FUNCTIONS_H_INCLUDED

#include <iostream>
#include <map>
#include "constants.h"
#include "cross_section.h"
#include "/home/eduardo/Desktop/Programming/cpp/2018-06-13_aux_maths/aux_maths.h"

//Contains all auxiliary functions used for the EDWA simulation
using wl_map = std::map<int, double>;

std::ostream& operator << (std::ostream& os, wl_map const& m);


namespace Utility
{

    //Converts wavelength in pm to frequency in Hz in free space (default) or in the medium
    double wl_to_freq(int const wl_, double n = 1.0) noexcept;
    
    //Converts frequency in Hz to wavelength in pm in free space (default) or in the medium
    int freq_to_wl(double const freq, double n = 1.0) noexcept;
    
    //Returns the cross-section according to the case
    double return_cross_section(int const wl, int const var) noexcept;
    
    //Returns the cross section adjusted, photon flux (1/s)
    double return_photon_flux(int const wl,double const P, 
                                 double const area, double const var) noexcept;
    
    //Converts a concentration related loss parameter from db/m to a number from 0 to 1. 
    //loss in units of db/m, length in m
    double return_conc_loss(double const loss, double const length, double const steps) noexcept;
    
    void update_wl_map_values(wl_map& m, std::vector<double> const& val) noexcept;

} //End of namespace Utility

#endif