#include "edwa_structs.h"

Simulation::Init_params::Init_params() noexcept
{
    A         = 2.0e-6 * 2.0e-6;   // 1.3 by 1.3 micron
    l         = 0.01;              // 1 cm
    NEr       = 4.45e26;           //Doping silica beyond... Chandrappan            
    NYb       = 4.45e26;           //Doping silica beyond... Chandrappan    
    A21       = 1.0 / (12.94e-3);  //Doping silica beyond... Chandrappan 
    A32       = (1.0e9);           //Finite-element... Di Pasquale
    A43       = (1.0e9);           //Finite-element... Di Pasquale
    A65       = 1.0 / (1.5e-3);    //Improved gain... Di Pasquale and Federighi
    S_start   = 1533000;           
    S_end     = 1533000;       
    n_signals = 1;     
    n_ASE     = 101;         
    steps     = 301;         
    step_size = l / static_cast<double>(steps-1);
    
    //Effect of pair induced... Di Pasquale and Federighi. linear interpolation
    Cup = 2.252632e-49 * NEr - 6.463158e-24;
    
    //Effect of pair induced...Assuming Dipole-Dipole interactions and neglecting Er-Er dipoles,
    //We get the following critical interaction distance
    double const var1_ {9.0 / (PI * PI * 16.0)};
    double const var2_ {Cup * (1.0 / A21)};
    double const R0_pow6 {var1_ * var2_ / NEr}; //Critical distance
    //Cross relaxation calculation
    double const RYb_Er {0.4e-9}; //Distance between erbium-ytterbium ions (0.4 nm) (Effect of pair...)
    double const var3_ {4.0 * PI / 3.0};
    double const var4_ {pow(RYb_Er, 3.0) * (1.0 / A65)};
    Ccr = var3_ * R0_pow6 / var4_;
    
    
    //Fill the channel spacing map
    int const step_size_ {100000 / (n_ASE-1)}; //Subdivide 1500-1600 nm in n_ASE-1 subintervals
    for (int i = 0; i < n_ASE; ++i)
    {
        int const curr_wl {1500000 + i * step_size_}; //current wavelength
        int const prev_wl {curr_wl - step_size_}; //previous wavelength
        double const spacing {Utility::wl_to_freq(prev_wl) - Utility::wl_to_freq(curr_wl)};
        channel_spacing.emplace(curr_wl, spacing);
        overlap.emplace(curr_wl, 0.4); //Becker et al chapter 6
    }
    
    overlap.emplace(976000, 0.64); //Becker et al chapter 6
    
    //Fill the signal map
    double const nm_step {(S_end - S_start)/static_cast<double>(n_signals)};
    
    for (int i = 0; i < n_signals; ++i)
    {
        double const curr_wl {S_start + static_cast<double>(i) * nm_step};
        double const pow {1.0e-6};                                      //1 microwatt or -30 dBm
        Ps0.emplace(curr_wl, pow);
    }

    //Fill the pump map
    Pp0_f.emplace(976000, 0.3);                                          //Forward pump of 300 mW
    Pp0_b.emplace(976000, 0.0);                                          //No backwards pump
    
    return;
}