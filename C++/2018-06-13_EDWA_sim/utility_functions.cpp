#include "utility_functions.h"

std::ostream& operator << (std::ostream& os, wl_map const& m)
{
    for (auto i = m.begin(); i != m.end(); ++i)
    {
        os<<(i->first)<<", "<<(i->second)<<'\n';
    }
    
    return os;
}

double Utility::wl_to_freq(int const wl_, double const n) noexcept
{
    double const wl {static_cast<double>(wl_) * 1.0e-12}; //Convert from picometres to metres
    double freq {C / (wl * n)};
    return freq;
}


int Utility::freq_to_wl(double const freq, double const n) noexcept
{
    double wl_ {C / (freq * n)};
    int wl {static_cast<int> (wl_ * 1.0e12)};
    return wl;
}


double Utility::return_cross_section(int const wl, Cross_sec const var) noexcept
{
    
    switch (var)
    {
        //Erbium
        case Cross_sec::Er_abs:
        return Simulation::erbium_abs(wl);
        
        case Cross_sec::Er_emi:
        return Simulation::erbium_emi(wl);
        
        //Ytterbium
        case Cross_sec::Yb_abs:
        return Simulation::ytterbium_abs(wl);
        
        case Cross_sec::Yb_emi:
        return Simulation::ytterbium_emi(wl);
        
        default:
        {
            std::cout<<"Error n get_cross_section: invalid case\n";
            return -1.0;
        }
    }
}


double Utility::return_photon_flux(int const wl, double const P, 
                                      double const area, Cross_sec const var) noexcept
{
    double const sigma {Utility::return_cross_section(wl, var)};
    double const freq {Utility::wl_to_freq(wl)};
    double const flux {sigma * P / (H * freq * area)};
    return flux;
}


double Utility::return_conc_loss(double const loss, double const length, double const steps) noexcept
{
    double const stepsize {length / steps};
    double const out_pc {pow(10.0, (loss/10.0))};
    double const loss_pc {(1.0-out_pc) * stepsize};
    return loss_pc;
}


//val and wl_map must have the same number of elements
void Utility::update_wl_map_values(wl_map& m, std::vector<double> const& val) noexcept
{
    auto m_begin {m.begin()};
    auto m_end {m.end()};
    std::size_t count {0};
    
    for (m_begin; m_begin != m_end; ++m_begin)
    {
        m_begin->second = val[count];
        ++count;
    }
    
    return;
}


//Converts two power values into gain
double Utility::power_to_gain(double const start_pow, double const end_pow) noexcept
{
    return 10.0 * log10(end_pow/start_pow);
}


//Converts power in Watts to dBm
double Utility::power_to_dBm(double const power) noexcept
{
    return 10.0 * log10(power*1000.0);
}




