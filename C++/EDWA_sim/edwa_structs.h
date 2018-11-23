#ifndef EDWA_STRUCTS_H_INCLUDED
#define EDWA_STRUCTS_H_INCLUDED

#include "utility_functions.h"

namespace Simulation
{
    
    struct Init_params
    {
        double width;                   //Width of the cross section in m
        double height;                  //Height of the cross section in m
        double A   {0.0};               //Area of the cross section in m^2
        double l   {0.0};               //Length of the waveguide
        double step_size {0.0};         //length/steps
        double NEr {0.0};               //Total Erbium ions/m^3
        double NYb {0.0};               //Total Ytterbium ions/m^3
        double A21 {0.0};               //Spontaneous transition rate from Erbium l2 to l1 in 1/s
        double A32 {0.0};               //Spontaneous transition rate from Erbium l3 to l2 in 1/s
        double A43 {0.0};               //Spontaneous transition rate from Erbium l4 to l3 in 1/s
        double A65 {0.0};               //Spontaneous transition rate from Erbium l6 to l5 in 1/s
        double Ccr {0.0};               //Cross relaxation coefficient from Yb to Er in m^3/s
        double Cup {0.0};               //Up conversion coefficient m^3/s
        double lpEr {0.0};              //Pump scattering loss due to erbium concentration
        double lsEr {0.0};              //Signal scatterng loss due to erbium concentration
        double lpYb {0.0};              //Pump scattering loss due to ytterbium concentration
        double lsYb {0.0};              //Signal scattering loss due to ytterbium concentration
        int    S_start   {0};           //Starting wavelength of the signals in pm
        int    S_end     {0};           //End wavelength of the signals pm
        int    n_signals {0};           //Total number of signals
        int    steps     {0};           //Number of steps/subdivision for the waveguide
        int    curr_step {0};
        wl_map channel_spacing;         //Channel spacing per wavelength of ASE signal
        wl_map Ps0;                      //Signal power by wavelength
        wl_map Pp0_f;                   //Initial Forward pump power by wavelength
        wl_map Pp0_b;                   //Initial Backward pump power by wavelength
        wl_map overlap;                 //Overlap between traverse modes and ions
        
        Init_params() noexcept; //Constructor, initializes the struct with default test data
        Init_params(Init_params const& p) = default;
        
        //Returns loss as db/m based on the concentration using linear interpolation of data
        //provided in Improved gain characteristics... Di Pasquale and M. Federighi
        void calculate_loss() noexcept; 
        
        //Recalculates dependent constants such as Cup and Ccr and loss
        void recalculate_constants() noexcept; 
    };
    
    //Contains the state of the amplifier at each step
    struct Step
    {
        double n1       {0.0};  //Erbium population in the ground state
        double n2       {0.0};  //Erbium population in the lasing level
        double n3       {0.0};  //Erbium population in the Non radiative relaxation level
        double n4       {0.0};  //Erbium population in the Up-conversion level
        double n5       {0.0};  //Ytterbium population in the ground state
        double n6       {0.0};  //Erbium Population in the excited state
        double lsEr     {0.0};  //Signal loss due to Erbium concentration (dB/m)
        double lsYb     {0.0};  //Signal loss due to Ytterbium concentration (dB/m)
        double lpEr     {0.0};  //Pump loss due to Erbium concentration (dB/m)
        double lpYb     {0.0};  //Pump loss due to Ytterbium concentration (dB/m)
        double residual {0.0};  //Residual of the solution
        unsigned int curr_step {0}; //Current step
        wl_map Pp_f;        //Forward pump power by wavelength (map<double, double>, wl in nm NOT m)
        wl_map Pp_b;        //Backwards pump power by wavelength
        wl_map Ps;          //Signals power
        wl_map PASE_f;      //Forward ASE power
        wl_map PASE_b;      //Backwards ASE power
        wl_map Pp_ASE_f;      //Forward Pump ASE power
        wl_map Pp_ASE_b;      //Backwards Pump ASE power
    };

} //End of namespace simulation


namespace Utility
{
    //Function that returns an array containing a vector with start iterators and a vector with
    //end iterators
    template <typename T>
    std::array<std::vector<T>, 2> return_iterators(Simulation::Step const& s)
    {
        //start iterators
        std::vector<T> start_it;
        start_it.emplace_back(s.Ps.begin());
        start_it.emplace_back(s.Pp_f.begin());
        start_it.emplace_back(s.Pp_b.begin());
        start_it.emplace_back(s.PASE_f.begin());
        start_it.emplace_back(s.PASE_b.begin());
        start_it.emplace_back(s.Pp_ASE_f.begin());
        start_it.emplace_back(s.Pp_ASE_b.begin());
        
        //end iterators
        std::vector<T> end_it;
        end_it.emplace_back(s.Ps.end());
        end_it.emplace_back(s.Pp_f.end());
        end_it.emplace_back(s.Pp_b.end());
        end_it.emplace_back(s.PASE_f.end());
        end_it.emplace_back(s.PASE_b.end());
        end_it.emplace_back(s.Pp_ASE_f.end());
        end_it.emplace_back(s.Pp_ASE_b.end());
        
        std::array<std::vector<T>, 2> output {start_it, end_it};
        return output;
    }
}



#endif /* EDWA_STRUCTS_H_INCLUDED */