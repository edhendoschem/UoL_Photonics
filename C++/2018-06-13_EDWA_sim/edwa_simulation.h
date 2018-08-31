#ifndef EDWA_SIMULATION_H_INCLUDED
#define EDWA_SIMULATION_H_INCLUDED

#include <vector>
#include <optional>
#include <fstream>
#include <string>
#include <regex>
#include "constants.h"
#include "cross_section.h"
#include "utility_functions.h"
#include "edwa_structs.h"


namespace Simulation
{
    class Result
    {
    public:
        //vars
        Simulation::Init_params p;
        std::vector<Simulation::Step> data;
        Simulation::Step start_step, end_step; //Steps stored for convergence comparison
        
        //Helper structs
        //Additional data needed in order to be able to use the static functions
        struct Data_pack
        {
            double W12, W13, W21, W65, W56;
            double Cup, Ccr;
            double NEr, NYb;
            double A21, A32, A43, A65;
            std::optional<Simulation::Init_params> params;
            std::optional<Simulation::Step> prev_z;
            std::optional<Simulation::Step> curr_z;
            std::optional<Simulation::Step> next_z;
        };
        
        //Abbreviations
        using Data_p = Simulation::Result::Data_pack;
        using u_int  = unsigned int;
        
        //Constructors
        Result(Simulation::Init_params const& initial_state) noexcept;
        Result() noexcept;
        
        //Functions
        static double Eq_1(std::vector<double> const& vars, Data_p const& p) noexcept; //dn1/dt
        static double Eq_2(std::vector<double> const& vars, Data_p const& p) noexcept; //dn2/dt
        static double Eq_3(std::vector<double> const& vars, Data_p const& p) noexcept; //dn3/dt
        static double Eq_4(std::vector<double> const& vars, Data_p const& p) noexcept; //dn5/dt
        static double Eq_5(std::vector<double> const& vars, Data_p const& p) noexcept; //NEr-n1-n2-n3-n4 = 0
        static double Eq_6(std::vector<double> const& vars, Data_p const& p) noexcept; //NYb - n5 - n6 = 0;
        
        //n finder (uses Newton's method to find the root of the equation system)
        void find_all_n(u_int const z, double const h, double const tol, unsigned int n_it) noexcept;
        
        //Pump, signal and ASE change
        //Individual case (wavelength wl_ in nm)
        double dPp_f_ind   (u_int const z, int const wl_, double const val) noexcept;
        double dPp_b_ind   (u_int const z, int const wl_, double const val) noexcept;
        double dPs_ind     (u_int const z, int const wl_, double const val) noexcept;
        double dPASE_f_ind (u_int const z, int const wl_, double const val) noexcept;
        double dPASE_b_ind (u_int const z, int const wl_, double const val) noexcept;
        double dPp_ASE_f_ind (u_int const z, int const wl_, double const val) noexcept;
        double dPp_ASE_b_ind (u_int const z, int const wl_, double const val) noexcept;
        
        //Simulation march
        void reset_start     ()                              noexcept; //Resets the starting parameters
        void reset_end       ()                              noexcept; //Resets the ending parameters
        void advance_signal  (u_int const z, int sign = 1) noexcept; //Advances 1 step the signal
        void advance_pump_f  (u_int const z, int sign = 1) noexcept; //Advances 1 step back pump
        void advance_pump_b  (u_int const z, int sign = 1) noexcept; //Advances 1 step forward pump
        void advance_ASE_f   (u_int const z, int sign = 1) noexcept; //Advances 1 step forward ASE
        void advance_ASE_b   (u_int const z, int sign = 1) noexcept; //Advances 1 step back ASE
        void advance_p_ASE_f (u_int const z, int sign = 1) noexcept; //Advances 1 step pump forward ASE
        void advance_p_ASE_b (u_int const z, int sign = 1) noexcept; //Advances 1 step pump back ASE
        
        
        //calculates all n in z and marches signal, pump and ASE to z+1
        void advance_step (u_int const z, double const h, double const tol, u_int n_it) noexcept; 
        void regress_step (u_int const z, double const h, double const tol, u_int n_it) noexcept;
        
        //Propagates forward and backwards to obtain the results
        void simulate(bool warn = false) noexcept;
        
        //Auxiliary Functions
        double  calculate_W(std::size_t const z, int const var) const noexcept;
        void    initialize_ASE() noexcept;
        void    report_step(u_int z, bool show_ASE = false) noexcept;
        void    save_data (std::string const filename_, bool dBm_units = false) noexcept;
        void    plot_data (std::string const data_file, std::string const plot_script) noexcept;
    };
    
    
    double df_dNEr(std::vector<double> const vars, Simulation::Init_params const p) noexcept;
    double df_dNYb(std::vector<double> const vars, Simulation::Init_params const p) noexcept; 
    std::vector<double> find_critical_point(Simulation::Init_params const p) noexcept;
    //Will attempt to find the optimal ratio and return an array with Ner, NYb,length and max gain
    std::vector<std::array<double, 4>> find_ratio(Simulation::Init_params const p) noexcept;
}

#endif  /* EDWA_SIMULATION_H_INCLUDED */