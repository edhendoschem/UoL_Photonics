#include "edwa_simulation.h"


//Constructor based on initial state
void Simulation::Result::initialize_ASE() noexcept
{
    for (Simulation::Step& S : data) 
    {
        int s_sz {150000 / (150)};
        for (int j = 0; j < 151; ++j)
        {
            int const wl  {1450000 + j*s_sz};
            int const wl2 {900000  + j*s_sz};
            S.PASE_b.emplace(wl, 0.0);
            S.PASE_f.emplace(wl, 0.0);
            S.Pp_ASE_b.emplace(wl2, 0.0);
            S.Pp_ASE_f.emplace(wl2, 0.0);
        }
    }
}


Simulation::Result::Result(Simulation::Init_params const& initial_state) noexcept
{
    p = initial_state;
    data.reserve(p.steps+1);
    for (auto i = 0; i < p.steps; ++i)
    {
        data.emplace_back(Simulation::Step{});
        data[i].Pp_f = p.Pp0_f;     //Initial forwards pump
        data[i].Pp_b = p.Pp0_b;     //Initial backwards pump
        data[i].Ps   = p.Ps0;       //Initial signal
        data[i].curr_step = i;      //Step number
        data[i].lsEr = p.lsEr;
        data[i].lpEr = p.lpEr;
        data[i].lsYb = p.lsYb;
        data[i].lpYb = p.lpYb;
    }
    
    initialize_ASE();           //Fills the ASE maps for every step with 0.0
    return;
}

//Default constructor
Simulation::Result::Result() noexcept
{
    Simulation::Init_params params_ {};
    p = params_;
    data.reserve(p.steps+1);
    for (auto i = 0; i < p.steps; ++i)
    {
        data.emplace_back(Simulation::Step{});
        data[i].Pp_f = p.Pp0_f;     //Initial forwards pump
        data[i].Pp_b = p.Pp0_b;     //Initial backwards pump
        data[i].Ps   = p.Ps0;       //Initial signal
        data[i].curr_step = i;      //Step number
        data[i].lsEr = p.lsEr;
        data[i].lpEr = p.lpEr;
        data[i].lsYb = p.lsYb;
        data[i].lpYb = p.lpYb;
    }
    
    initialize_ASE();           //Fills the ASE maps for every step with 0.0
    return;
}


double Simulation::Result::calculate_W(std::size_t const z, int const var) const noexcept
{
    double result {0.0};
    
    if (var == 56 || var == 65)
    {
        //Forward Pump contribution
        auto st1 {data[z].Pp_f.cbegin()}; //Iterate over forward pump in case there are more than one
        auto end1 {data[z].Pp_f.cend()};
        Cross_sec var_ {var == 56 ? Cross_sec::Yb_abs : Cross_sec::Yb_emi};
    
        for (st1; st1 != end1; ++st1)
        {
            int const wl1        {(st1->first)};
            if (wl1 == 1480000) continue; //Ignore 1480 nm pump
            double const pump_pow1  {st1->second};
            double const flux       {Utility::return_photon_flux(
                wl1,
                pump_pow1,
                p.A,
                var_
            )};
        
            result += flux;
        }
        
        //Backwards pump contribution
        auto st2 {data[z].Pp_b.cbegin()}; //Iterate over forward pump in case there are more than one
        auto end2 {data[z].Pp_b.cend()};
    
        for (st2; st2 != end2; ++st2)
        {
            int const wl2        {(st2->first)};
            if (wl2 == 1480000) continue; //Ignore 1480 nm wavelength
            double const pump_pow2  {st2->second};
            double const flux       {Utility::return_photon_flux(
                wl2,
                pump_pow2,
                p.A,
                var_
            )};
        
            result += flux;
        }
        
        
        //Forwards pump ASE contribution
        auto st3 {data[z].Pp_ASE_f.cbegin()}; //Iterate over forward pump in case there are more than one
        auto end3 {data[z].Pp_ASE_f.cend()};
    
        for (st3; st3 != end3; ++st3)
        {
            int const wl3        {(st3->first)};
            if (wl3 == 1480000) continue; //Ignore 1480 nm wavelength
            double const pump_pow3  {st3->second};
            double const flux       {Utility::return_photon_flux(
                wl3,
                pump_pow3,
                p.A,
                var_
            )};
        
            result += flux;
        }


        //Backwards pump ASE contribution
        auto st4 {data[z].Pp_ASE_b.cbegin()}; //Iterate over forward pump in case there are more than one
        auto end4 {data[z].Pp_ASE_b.cend()};
    
        for (st4; st4 != end4; ++st4)
        {
            int const wl4        {(st4->first)};
            if (wl4 == 1480000) continue; //Ignore 1480 nm wavelength
            double const pump_pow4  {st4->second};
            double const flux       {Utility::return_photon_flux(
                wl4,
                pump_pow4,
                p.A,
                var_
            )};
        
            result += flux;
        }
    
    } else {
        //Signal contribution
        auto st1 {data[z].Ps.cbegin()}; //Iterate over forward pump in case there are more than one
        auto end1 {data[z].Ps.cend()};
        Cross_sec var_ {(var == 12 || var == 13) ? Cross_sec::Er_abs : Cross_sec::Er_emi};
        
        for (st1; st1 != end1; ++st1)
        {
            if (var == 13) break;
            int const wl1        {(st1->first)};
            double const sig_pow1   {st1->second};
            double const flux       {Utility::return_photon_flux(
                wl1,
                sig_pow1,
                p.A,
                var_
            )};
        
            result += flux;
        }
        
        //Forwards ASE contribution
        auto st2 {data[z].PASE_f.cbegin()}; //Iterate over forward pump in case there are more than one
        auto end2 {data[z].PASE_f.cend()};
    
        for (st2; st2 != end2; ++st2)
        {
            if (var == 13) break;
            int const wl2        {(st2->first)};
            double const PASE_f_pow {st2->second};
            double const flux       {Utility::return_photon_flux(
                wl2,
                PASE_f_pow,
                p.A,
                var_
            )};
        
            result += flux;
        }
        
        //Backwards ASE contribution
        auto st3 {data[z].PASE_b.cbegin()}; //Iterate over forward pump in case there are more than one
        auto end3 {data[z].PASE_b.cend()};
    
        for (st3; st3 != end3; ++st3)
        {
            if (var == 13) break;
            int const wl3        {(st3->first)};
            double const PASE_b_pow {st3->second};
            double const flux       {Utility::return_photon_flux(
                wl3,
                PASE_b_pow,
                p.A,
                var_
            )};
        
            result += flux;
        }
        
        //Forward 1480 nm Pump contribution
        auto st4 {data[z].Pp_f.cbegin()}; 
        auto end4 {data[z].Pp_f.cend()};
    
        for (st4; st4 != end4; ++st4)
        {
            if (var == 13) break;
            int const wl4        {(st4->first)};
            if (wl4 != 1480000) continue; //Ignore wavelengths different from 1480 nm
            double const pump_pow4  {st4->second};
            double const flux       {Utility::return_photon_flux(
            wl4,
            pump_pow4,
            p.A,
            Cross_sec::Er_abs
            )};
        
            result += flux;
        }
            
            
        //Backwards 1480 pump contribution
        auto st5 {data[z].Pp_b.cbegin()}; 
        auto end5 {data[z].Pp_b.cend()};
    
        for (st5; st5 != end5; ++st5)
        {
            if (var == 13) break;
            int const wl5        {(st5->first)};
            if (wl5 != 1480000) continue; //Ignore wavelengths different from 1480 nm
            double const pump_pow5  {st5->second};
            double const flux       {Utility::return_photon_flux(
            wl5,
            pump_pow5,
            p.A,
            Cross_sec::Er_abs
            )};
        
            result += flux;
        }
        
        
        if (var == 13)
        {
            //Forward Pump contribution
            auto st6 {data[z].Pp_f.cbegin()}; 
            auto end6 {data[z].Pp_f.cend()};
    
            for (st6; st6 != end6; ++st6)
            {
                int const wl6        {(st6->first)};
                if (wl6 == 1480000) continue; //Ignore 1480 nm
                double const pump_pow6  {st6->second};
                double const flux       {Utility::return_photon_flux(
                    wl6,
                    pump_pow6,
                    p.A,
                    Cross_sec::Yb_abs
                )};
        
                result += flux;
            }
            
            
            //Backwards pump
            auto st7 {data[z].Pp_b.cbegin()}; 
            auto end7 {data[z].Pp_b.cend()};
    
            for (st7; st7 != end7; ++st7)
            {
                int const wl7        {(st7->first)};
                if (wl7 == 1480000) continue; //Ignore 1480 nm
                double const pump_pow7  {st7->second};
                double const flux       {Utility::return_photon_flux(
                    wl7,
                    pump_pow7,
                    p.A,
                    Cross_sec::Yb_abs
                )};
        
                result += flux;
            }
            
            //Pump forwards pump ASE contribution
            auto st8 {data[z].Pp_ASE_f.cbegin()}; 
            auto end8 {data[z].Pp_ASE_f.cend()};
    
            for (st8; st8 != end8; ++st8)
            {
                int const wl8        {(st8->first)};
                if (wl8 < 962) continue;
                if (wl8 > 990) break;
                double const pump_pow8  {st8->second};
                double const flux       {Utility::return_photon_flux(
                    wl8,
                    pump_pow8,
                    p.A,
                    Cross_sec::Yb_abs
                )};
        
                result += flux;
            }
            
            
            //Pump backwards pump ASE contribution
            auto st9 {data[z].Pp_ASE_b.cbegin()}; 
            auto end9 {data[z].Pp_ASE_b.cend()};
    
            for (st9; st9 != end9; ++st9)
            {
                int const wl9        {(st9->first)};
                if (wl9 < 962) continue;
                if (wl9 > 990) break;
                double const pump_pow9  {st9->second};
                double const flux       {Utility::return_photon_flux(
                    wl9,
                    pump_pow9,
                    p.A,
                    Cross_sec::Yb_abs
                )};
        
                result += flux;
            }
        
        } //End of if var == 13
    } //End of else
    
    return result;
}



//Individual calculations
double Simulation::Result::dPp_f_ind (u_int const z, int const wl_, double const val) noexcept
{
    double const n1   {data[z].n1};
    double const n5   {data[z].n5};
    double const n6   {data[z].n6};
    double const lpEr {data[z].lpEr};
    double const lpYb {data[z].lpYb};
    double const tot_loss {Utility::return_conc_loss(lpEr + lpYb)};
    double output     {0.0};
    
    int const wl {wl_};
    double const power {val};
    double const sigEr_13 {Utility::return_cross_section(wl, Cross_sec::Er_abs)};
    double const sigYb_56 {Utility::return_cross_section(wl, Cross_sec::Yb_abs)};
    double const sigYb_65 {Utility::return_cross_section(wl, Cross_sec::Yb_emi)};
    double const overlap {p.overlap.at(wl_)};
    output = power * overlap * (sigYb_65 * n6 - (sigEr_13 * n1 + sigYb_56 * n5)) - (tot_loss) * power;

    return output;
}


double Simulation::Result::dPp_b_ind (u_int const z, int const wl_, double const val) noexcept
{
    double const n1   {data[z].n1};
    double const n5   {data[z].n5};
    double const n6   {data[z].n6};
    double const lpEr {data[z].lpEr};
    double const lpYb {data[z].lpYb};
    double const tot_loss {Utility::return_conc_loss(lpEr + lpYb)};
    double output     {0.0};
    
    int const wl {wl_};
    double const power {val};
    double const sigEr_13 {Utility::return_cross_section(wl, Cross_sec::Er_abs)};
    double const sigYb_56 {Utility::return_cross_section(wl, Cross_sec::Yb_abs)};
    double const sigYb_65 {Utility::return_cross_section(wl, Cross_sec::Yb_emi)};
    double const overlap {p.overlap.at(wl_)};
    output = -power * overlap * (sigYb_65 * n6 - (sigEr_13 * n1 + sigYb_56 * n5)) + (tot_loss) * power;

    return output;
}

double Simulation::Result::dPs_ind (u_int const z, int const wl_, double const val) noexcept
{
    double const n1   {data[z].n1};
    double const n2   {data[z].n2};
    double const lsEr {data[z].lsEr};
    double const lsYb {data[z].lsYb};
    double const tot_loss {Utility::return_conc_loss(lsEr + lsYb)};
    double output     {0.0};
    
    int const wl {wl_};
    double const power {val};
    double const sigEr_21 {Utility::return_cross_section(wl, Cross_sec::Er_emi)};
    double const sigEr_12 {Utility::return_cross_section(wl, Cross_sec::Er_abs)};
    double const overlap {p.overlap.at(wl_)};
    
    output = power * overlap * (sigEr_21 * n2 - sigEr_12 * n1) - (tot_loss) * power;
    return output;                                          
} //End of dPs_ind


double Simulation::Result::dPASE_f_ind (u_int const z, int const wl_, double const val) noexcept
{
    double const n1 {data[z].n1};
    double const n2 {data[z].n2};
    double const lsEr {data[z].lsEr};
    double const lsYb {data[z].lsYb};
    double const tot_loss {Utility::return_conc_loss(lsEr + lsYb)};
    double output {0.0};
    
    int    const wl {wl_};
    double const power {val};
    double const sigEr_21 {Utility::return_cross_section(wl, Cross_sec::Er_emi)};
    double const sigEr_12 {Utility::return_cross_section(wl, Cross_sec::Er_abs)};
    double const v {Utility::wl_to_freq(wl, 1.0)};
    double const delta_v {p.channel_spacing.at(wl_)};
    double const overlap {p.overlap.at(wl_)};
    double const a {power * overlap * (sigEr_21 * n2 - sigEr_12 * n1)};
    double const b {2.0 * H * v * delta_v * overlap * sigEr_21 * n2};
    double const c {(tot_loss) * power};
    output = a + b - c;
    
    return output;
} //End dPASE_f_ind


double Simulation::Result::dPASE_b_ind (u_int const z, int const wl_, double const val) noexcept
{
    double const n1 {data[z].n1};
    double const n2 {data[z].n2};
    double const lsEr {data[z].lsEr};
    double const lsYb {data[z].lsYb};
    double const tot_loss {Utility::return_conc_loss(lsEr + lsYb)};
    double output {0.0};
    
    int    const wl {wl_};
    double const power {val};
    double const sigEr_21 {Utility::return_cross_section(wl, Cross_sec::Er_emi)};
    double const sigEr_12 {Utility::return_cross_section(wl, Cross_sec::Er_abs)};
    double const v {Utility::wl_to_freq(wl, 1.0)};
    double const delta_v {p.channel_spacing.at(wl_)};
    double const overlap {p.overlap.at(wl_)};
    double const a {power * overlap * (sigEr_21 * n2 - sigEr_12 * n1)};
    double const b {2.0 * H * v * delta_v * overlap * sigEr_21 * n2};
    double const c {(tot_loss) * power};
    output = -a - b + c;

    return output;
} //End of dPASE_b_ind


double Simulation::Result::dPp_ASE_f_ind (u_int const z, int const wl_, double const val) noexcept
{

    double const n1 {data[z].n1};
    double const n5 {data[z].n5};
    double const n6 {data[z].n6};
    double const lpEr {data[z].lpEr};
    double const lpYb {data[z].lpYb};
    double const tot_loss {Utility::return_conc_loss(lpEr + lpYb)};
    double output {0.0};
    
    int    const wl {wl_};
    double const power {val};
    double const sigYb_65 {Utility::return_cross_section(wl, Cross_sec::Yb_emi)};
    double const sigYb_56 {Utility::return_cross_section(wl, Cross_sec::Yb_abs)};
    double sigEr_13 {0.0};
    if ( wl >= 962 && wl <= 990)
    {
        sigEr_13 = Utility::return_cross_section(wl, Cross_sec::Er_abs);
    }
    
    double const v {Utility::wl_to_freq(wl, 1.0)};
    double const delta_v {p.channel_spacing.at(wl_)};
    double const overlap {p.overlap.at(wl_)};
    double const a {power * overlap * (sigYb_65 * n6 - (sigEr_13 * n1 + sigYb_56 * n5))};
    double const b {2.0 * H * v * delta_v * overlap * sigYb_65 * n6};
    double const c {(tot_loss) * power};
    output = a + b - c;
    
    return output;
} //end of dPp_ASE_f_ind


double Simulation::Result::dPp_ASE_b_ind (u_int const z, int const wl_, double const val) noexcept
{
    double const n1 {data[z].n1};
    double const n5 {data[z].n5};
    double const n6 {data[z].n6};
    double const lpEr {data[z].lpEr};
    double const lpYb {data[z].lpYb};
    double const tot_loss {Utility::return_conc_loss(lpEr + lpYb)};
    double output {0.0};
    
    int    const wl {wl_};
    double const power {val};
    double const sigYb_65 {Utility::return_cross_section(wl, Cross_sec::Yb_emi)};
    double const sigYb_56 {Utility::return_cross_section(wl, Cross_sec::Yb_abs)};
    double sigEr_13 {0.0};
    if ( wl >= 962 && wl <= 990)
    {
        sigEr_13 = Utility::return_cross_section(wl, Cross_sec::Er_abs);
    }
    
    double const v {Utility::wl_to_freq(wl, 1.0)};
    double const delta_v {p.channel_spacing.at(wl_)};
    double const overlap {p.overlap.at(wl_)};
    double const a {power * overlap * (sigYb_65 * n6 - (sigEr_13 * n1 + sigYb_56 * n5))};
    double const b {2.0 * H * v * delta_v * overlap * sigYb_65 * n6};
    double const c {(tot_loss) * power};
    output = -a - b + c;
    
    return output;
} //end of dPp_ASE_b_ind


//Prints information about this step
void Simulation::Result::report_step(u_int z, bool show_ASE) noexcept
{
    std::cout<<"Step "<<z<<" data:\n";
    std::cout<<"Residual = "<<data[z].residual<<'\n';
    std::cout<<"n6 = "<<data[z].n6<<'\n';
    std::cout<<"n5 = "<<data[z].n5<<'\n';
    std::cout<<"n4 = "<<data[z].n4<<'\n';
    std::cout<<"n3 = "<<data[z].n3<<'\n';
    std::cout<<"n2 = "<<data[z].n2<<'\n';
    std::cout<<"n1 = "<<data[z].n1<<'\n';
    std::cout<<"-----------------------\n";
    std::cout<<".........Ps.........\n"<<data[z].Ps<<'\n';
    std::cout<<".........Pp_f.........\n"<<data[z].Pp_f<<'\n';
    std::cout<<".........Pp_b.........\n"<<data[z].Pp_b<<'\n';
    if (show_ASE) std::cout<<".........PASE_f (1533 nm).........\n"<<data[z].PASE_f.at(1533000)<<'\n';
    if (show_ASE) std::cout<<".........PASE_b (1533 nm).........\n"<<data[z].PASE_b.at(1533000)<<'\n'; 
    /*
    if (show_ASE) std::cout<<".........PASE_f.........\n"<<data[z].PASE_f<<'\n';
    if (show_ASE) std::cout<<".........PASE_b.........\n"<<data[z].PASE_b<<'\n';
    */
    std::cout<<"=======================\n";
}


//Advances signal one step forwards or backwards
void Simulation::Result::advance_signal(u_int const z, int sign) noexcept
{
    double const h {p.step_size};
    //double const h {p.step_size/101};
    auto st_it {data[z].Ps.begin()}; //start iterator
    auto ed_it {data[z].Ps.end()};   //End iterator
    double const sign_d {static_cast<double>(sign)};
    
    //Starting values of Runge - Kutta
    for (st_it; st_it != ed_it; ++st_it)
    {
        int const wl {st_it->first};
        double P0 {st_it->second};
        double const k1 {h * sign_d * dPs_ind(z, wl, P0)};
        double const k2 {h * sign_d * dPs_ind(z, wl, (P0 + k1/2.0))};
        double const k3 {h * sign_d * dPs_ind(z, wl, (P0 + k2/2.0))};
        double const k4 {h * sign_d * dPs_ind(z, wl, (P0 + k3))};
        data[z+sign*1].Ps.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        if (data[z+sign*1].Ps.at(wl) < 0.0) data[z+sign*1].Ps.at(wl) = 0.0;
    }
    
    return;
}


//Advances forwad pump one step forwards or backwards
void Simulation::Result::advance_pump_f(u_int const z, int sign) noexcept
{
    double const h {p.step_size};
    auto st_it {data[z].Pp_f.begin()}; //start iterator
    auto ed_it {data[z].Pp_f.end()};   //End iterator
    double const sign_d {static_cast<double>(sign)};
    
    //Starting values of Runge - Kutta
    for (st_it; st_it != ed_it; ++st_it)
    {
        int const wl {st_it->first};
        if (wl != 1480000)
        {    
            double const P0 {st_it->second};
            double const k1 {h * sign_d * dPp_f_ind(z, wl, P0)};
            double const k2 {h * sign_d * dPp_f_ind(z, wl, (P0 + k1/2.0))};
            double const k3 {h * sign_d * dPp_f_ind(z, wl, (P0 + k2/2.0))};
            double const k4 {h * sign_d * dPp_f_ind(z, wl, (P0 + k3))};
            data[z+sign*1].Pp_f.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            if (data[z+sign*1].Pp_f.at(wl) < 0.0) data[z+sign*1].Pp_f.at(wl) = 0.0;
        } else {
            double const P0 {st_it->second};
            double const k1 {h * sign_d * dPs_ind(z, wl, P0)};
            double const k2 {h * sign_d * dPs_ind(z, wl, (P0 + k1/2.0))};
            double const k3 {h * sign_d * dPs_ind(z, wl, (P0 + k2/2.0))};
            double const k4 {h * sign_d * dPs_ind(z, wl, (P0 + k3))};
            data[z+sign*1].Pp_f.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            if (data[z+sign*1].Pp_f.at(wl) < 0.0) data[z+sign*1].Pp_f.at(wl) = 0.0;
        }
    }
    
    return;
}


//Advances backwards pump one step forward or backwards
void Simulation::Result::advance_pump_b(u_int const z, int sign) noexcept
{
    double const h {p.step_size};
    auto st_it {data[z].Pp_b.begin()}; //start iterator
    auto ed_it {data[z].Pp_b.end()};   //End iterator
    double const sign_d {static_cast<double>(sign)};
    
    //Starting values of Runge - Kutta
    for (st_it; st_it != ed_it; ++st_it)
    {
        int const wl {st_it->first};
        if (wl != 1480000)
        {    
            double const P0 {st_it->second};
            double const k1 {h * sign_d * dPp_b_ind(z, wl, P0)};
            double const k2 {h * sign_d * dPp_b_ind(z, wl, (P0 + k1/2.0))};
            double const k3 {h * sign_d * dPp_b_ind(z, wl, (P0 + k2/2.0))};
            double const k4 {h * sign_d * dPp_b_ind(z, wl, (P0 + k3))};
            data[z+sign*1].Pp_b.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            if (data[z+sign*1].Pp_b.at(wl) < 0.0) data[z+sign*1].Pp_b.at(wl) = 0.0;
        } else {
            double const P0 {st_it->second};
            double const k1 {h * sign_d * -1.0 * dPs_ind(z, wl, P0)};
            double const k2 {h * sign_d * -1.0 * dPs_ind(z, wl, (P0 + k1/2.0))};
            double const k3 {h * sign_d * -1.0 * dPs_ind(z, wl, (P0 + k2/2.0))};
            double const k4 {h * sign_d * -1.0 * dPs_ind(z, wl, (P0 + k3))};
            data[z+sign*1].Pp_b.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            if (data[z+sign*1].Pp_b.at(wl) < 0.0) data[z+sign*1].Pp_b.at(wl) = 0.0;
        }
    }
    
    return;
}



//Advances forward ASE one step forward or backwards
void Simulation::Result::advance_ASE_f(u_int const z, int sign) noexcept
{
    double const h {p.step_size};
    auto st_it {data[z].PASE_f.begin()}; //start iterator
    auto ed_it {data[z].PASE_f.end()};   //End iterator
    double const sign_d {static_cast<double>(sign)};
    
    //Starting values of Runge - Kutta
    for (st_it; st_it != ed_it; ++st_it)
    {
        int const wl {st_it->first};
        double const P0 {st_it->second};
        double const k1 {h * sign_d * dPASE_f_ind(z, wl, P0)};
        double const k2 {h * sign_d * dPASE_f_ind(z, wl, (P0 + k1/2.0))};
        double const k3 {h * sign_d * dPASE_f_ind(z, wl, (P0 + k2/2.0))};
        double const k4 {h * sign_d * dPASE_f_ind(z, wl, (P0 + k3))};
        data[z+sign*1].PASE_f.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        if (data[z+sign*1].PASE_f.at(wl) < 0.0) data[z+sign*1].PASE_f.at(wl) = 0.0;
    }
    
    return;
}


//Advances backward ASE one step forward or backwards
void Simulation::Result::advance_ASE_b(u_int const z, int sign) noexcept
{
    double const h {p.step_size};
    auto st_it {data[z].PASE_b.begin()}; //start iterator
    auto ed_it {data[z].PASE_b.end()};   //End iterator
    double const sign_d {static_cast<double>(sign)};
    
    //Starting values of Runge - Kutta
    for (st_it; st_it != ed_it; ++st_it)
    {
        int const wl {st_it->first};
        double const P0 {st_it->second};
        double const k1 {h * sign_d * dPASE_b_ind(z, wl, P0)};
        double const k2 {h * sign_d * dPASE_b_ind(z, wl, (P0 + k1/2.0))};
        double const k3 {h * sign_d * dPASE_b_ind(z, wl, (P0 + k2/2.0))};
        double const k4 {h * sign_d * dPASE_b_ind(z, wl, (P0 + k3))};
        data[z+sign*1].PASE_b.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        if (data[z+sign*1].PASE_b.at(wl) < 0.0) data[z+sign*1].PASE_b.at(wl) = 0.0;
    }
    
    return;
}


//Advances pump forward ASE one step forward or backwards
void Simulation::Result::advance_p_ASE_f(u_int const z, int sign) noexcept
{
    double const h {p.step_size};
    auto st_it {data[z].Pp_ASE_f.begin()}; //start iterator
    auto ed_it {data[z].Pp_ASE_f.end()};   //End iterator
    double const sign_d {static_cast<double>(sign)};
    
    //Starting values of Runge - Kutta
    for (st_it; st_it != ed_it; ++st_it)
    {
        int const wl {st_it->first};
        double const P0 {st_it->second};
        double const k1 {h * sign_d * dPp_ASE_f_ind(z, wl, P0)};
        double const k2 {h * sign_d * dPp_ASE_f_ind(z, wl, (P0 + k1/2.0))};
        double const k3 {h * sign_d * dPp_ASE_f_ind(z, wl, (P0 + k2/2.0))};
        double const k4 {h * sign_d * dPp_ASE_f_ind(z, wl, (P0 + k3))};
        data[z+sign*1].Pp_ASE_f.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        if (data[z+sign*1].Pp_ASE_f.at(wl) < 0.0) data[z+sign*1].Pp_ASE_f.at(wl) = 0.0;
    }
    
    return;
}


//Advances pump backward ASE one step forward or backwards
void Simulation::Result::advance_p_ASE_b(u_int const z, int sign) noexcept
{
    double const h {p.step_size};
    auto st_it {data[z].Pp_ASE_b.begin()}; //start iterator
    auto ed_it {data[z].Pp_ASE_b.end()};   //End iterator
    double const sign_d {static_cast<double>(sign)};
    
    //Starting values of Runge - Kutta
    for (st_it; st_it != ed_it; ++st_it)
    {
        int const wl {st_it->first};
        double const P0 {st_it->second};
        double const k1 {h * sign_d * dPp_ASE_b_ind(z, wl, P0)};
        double const k2 {h * sign_d * dPp_ASE_b_ind(z, wl, (P0 + k1/2.0))};
        double const k3 {h * sign_d * dPp_ASE_b_ind(z, wl, (P0 + k2/2.0))};
        double const k4 {h * sign_d * dPp_ASE_b_ind(z, wl, (P0 + k3))};
        data[z+sign*1].Pp_ASE_b.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        if (data[z+sign*1].Pp_ASE_b.at(wl) < 0.0) data[z+sign*1].Pp_ASE_b.at(wl) = 0.0;
    }
    
    return;
}


double Simulation::Result::Eq_1(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n1 {vars[0]};
    double const n2 {vars[1]};
    double const n3 {vars[2]};
    double const n4 {vars[3]};
    double const n6 {vars[5]};
    
    double const W12 {p.W12};
    double const W13 {p.W13};
    double const W21 {p.W21};
    double const t1  {-(W12 + W13) * n1};
    double const t2  {(p.A21 + W21) * n2};
    double const t3  {p.Cup * pow(n2, 2.0)};
    double const t4  {-p.Ccr * n1 * n4};
    double const t5  {p.Cup * pow(n3, 2.0)};
    double const t6  {-p.Ccr * n1 * n6};
    
    return t1 + t2 + t3 + t4 + t5 + t6;
}


double Simulation::Result::Eq_2(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n1 {vars[0]};
    double const n2 {vars[1]};
    double const n3 {vars[2]};
    double const n4 {vars[3]};
    
    double const W12 {p.W12};
    double const W21 {p.W21};
    double const t1  {W12 * n1 - (p.A21 + W21) * n2};
    double const t2  {p.A32 * n3 - 2.0 * p.Cup * pow(n2, 2.0)};
    double const t3  {2.0 * p.Ccr * n1 * n4};
    
    return t1 + t2 + t3;
}


double Simulation::Result::Eq_3(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n1 {vars[0]};
    double const n3 {vars[2]};
    double const n4 {vars[3]};
    double const n6 {vars[5]};
    
    double const W13 {p.W13};
    double const t1  {W13 * n1 - p.A32 * n3};
    double const t2  {p.A43 * n4 - 2.0 * p.Cup * pow(n3, 2.0)};
    double const t3  {p.Ccr * n1 * n6};
    
    return t1 + t2 + t3;
}


double Simulation::Result::Eq_4(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n1 {vars[0]};
    double const n5 {vars[4]};
    double const n6 {vars[5]};
    
    double const W56 {p.W56};
    double const W65 {p.W65};
    double const t1  {-W56 * n5 + (p.A65 + W65) * n6};
    double const t2  {p.Ccr * n1 * n6};
    
    return t1 + t2;
}


double Simulation::Result::Eq_5(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n1 {vars[0]};
    double const n2 {vars[1]};
    double const n3 {vars[2]};
    double const n4 {vars[3]};

    return p.NEr - n1 - n2 - n3 - n4;
}


double Simulation::Result::Eq_6(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n5 {vars[4]};
    double const n6 {vars[5]};
    return p.NYb - n5 - n6;
}


double Simulation::Result::dEq_1_dn1(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n4 {vars[3]};
    double const n6 {vars[5]};
    return -(p.W12 + p.W13) - p.Ccr * (n4 + n6);
}

double Simulation::Result::dEq_1_dn2(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n2 {vars[1]};
    return p.A21 + p.W21 + 2.0 * p.Cup * n2;
}

double Simulation::Result::dEq_1_dn3(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n3 {vars[2]};
    return 2.0 * p.Cup * n3;
}

double Simulation::Result::dEq_1_dn4(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n1 {vars[0]};
    return -p.Ccr * n1;
}

double Simulation::Result::dEq_1_dn5(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_1_dn6(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n1 {vars[0]};
    return -p.Ccr * n1;
}

double Simulation::Result::dEq_2_dn1(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n4 {vars[3]};
    return p.W12 + 2.0 * p.Ccr * n4;
}

double Simulation::Result::dEq_2_dn2(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n2 {vars[1]};    
    return -(p.A21 + p.W21) - 4.0 * p.Cup * n2;
}

double Simulation::Result::dEq_2_dn3(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return p.A32;
}

double Simulation::Result::dEq_2_dn4(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n1 {vars[0]};
    return 2.0 * p.Ccr * n1;
}

double Simulation::Result::dEq_2_dn5(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_2_dn6(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_3_dn1(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n6 {vars[5]};
    return p.W13 + p.Ccr * n6;
}

double Simulation::Result::dEq_3_dn2(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_3_dn3(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n3 {vars[2]};
    return -p.A32 - 4.0 * p.Cup * n3;
}

double Simulation::Result::dEq_3_dn4(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return p.A43;
}

double Simulation::Result::dEq_3_dn5(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_3_dn6(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n1 {vars[0]};
    return p.Ccr * n1;
}

double Simulation::Result::dEq_4_dn1(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n6 {vars[5]};
    return p.Ccr * n6;
}

double Simulation::Result::dEq_4_dn2(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_4_dn3(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_4_dn4(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_4_dn5(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return -p.W56;
}

double Simulation::Result::dEq_4_dn6(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n1 {vars[0]};
    return p.A65 + p.W65 + p.Ccr * n1;
}

double Simulation::Result::dEq_5_dn1(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return -1.0;
}

double Simulation::Result::dEq_5_dn2(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return -1.0;
}

double Simulation::Result::dEq_5_dn3(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return -1.0;
}

double Simulation::Result::dEq_5_dn4(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return -1.0;
}

double Simulation::Result::dEq_5_dn5(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_5_dn6(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_6_dn1(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_6_dn2(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_6_dn3(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_6_dn4(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return 0.0;
}

double Simulation::Result::dEq_6_dn5(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return -1.0;
}

double Simulation::Result::dEq_6_dn6(std::vector<double> const& vars, Data_p const& p) noexcept
{
    return -1.0;
}


void Simulation::Result::find_all_n(u_int const z, double const h, double const tol, u_int const n_it) noexcept
{
    Data_p p_;
    p_.W12 = calculate_W(z, 12);
    p_.W13 = calculate_W(z, 13);
    p_.W21 = calculate_W(z, 21);
    p_.W65 = calculate_W(z, 65);
    p_.W56 = calculate_W(z, 56);
    p_.Cup = p.Cup;
    p_.Ccr = p.Ccr;
    p_.NEr = p.NEr;
    p_.NYb = p.NYb;
    p_.A21 = p.A21;
    p_.A32 = p.A32;
    p_.A43 = p.A43;
    p_.A65 = p.A65;
    
    using f_ptr = double (*)(std::vector<double> const&, Data_p const&);
    std::vector<f_ptr> f {Simulation::Result::Eq_1, 
                          Simulation::Result::Eq_2,
                          Simulation::Result::Eq_3,
                          Simulation::Result::Eq_4,
                          Simulation::Result::Eq_5,
                          Simulation::Result::Eq_6};
    std::vector<f_ptr> df {Simulation::Result::dEq_1_dn1, 
                           Simulation::Result::dEq_2_dn1,
                           Simulation::Result::dEq_3_dn1,
                           Simulation::Result::dEq_4_dn1,
                           Simulation::Result::dEq_5_dn1,
                           Simulation::Result::dEq_6_dn1,
                           Simulation::Result::dEq_1_dn2,
                           Simulation::Result::dEq_2_dn2,
                           Simulation::Result::dEq_3_dn2,
                           Simulation::Result::dEq_4_dn2,
                           Simulation::Result::dEq_5_dn2,
                           Simulation::Result::dEq_6_dn2,
                           Simulation::Result::dEq_1_dn3,
                           Simulation::Result::dEq_2_dn3,
                           Simulation::Result::dEq_3_dn3,
                           Simulation::Result::dEq_4_dn3,
                           Simulation::Result::dEq_5_dn3,
                           Simulation::Result::dEq_6_dn3,
                           Simulation::Result::dEq_1_dn4,
                           Simulation::Result::dEq_2_dn4,
                           Simulation::Result::dEq_3_dn4,
                           Simulation::Result::dEq_4_dn4,
                           Simulation::Result::dEq_5_dn4,
                           Simulation::Result::dEq_6_dn4,
                           Simulation::Result::dEq_1_dn5,
                           Simulation::Result::dEq_2_dn5,
                           Simulation::Result::dEq_3_dn5,
                           Simulation::Result::dEq_4_dn5,
                           Simulation::Result::dEq_5_dn5,
                           Simulation::Result::dEq_6_dn5,
                           Simulation::Result::dEq_1_dn6,
                           Simulation::Result::dEq_2_dn6,
                           Simulation::Result::dEq_3_dn6,
                           Simulation::Result::dEq_4_dn6,
                           Simulation::Result::dEq_5_dn6,
                           Simulation::Result::dEq_6_dn6};
                          
    std::vector<double> x0(6, 0.0);
    if (z == 0 && 
        data[z].n1 == 0.0 &&
        data[z].n2 == 0.0 &&
        data[z].n3 == 0.0 &&
        data[z].n4 == 0.0 &&
        data[z].n5 == 0.0 &&
        data[z].n6 == 0.0)
    {
        //Initial guess
        x0[0] = 0.90*p.NEr;
        x0[1] = 0.90*p.NEr;
        x0[2] = 0.90*p.NEr;
        x0[3] = 0.90*p.NEr;
        x0[4] = 0.90*p.NYb;
        x0[5] = 0.90*p.NYb;
    } else if (z == 0)
    {
        //If the simulation has been run previously then use curr vals as initial guess
        x0[0] = data[z].n1 > -0.1 ? data[z].n1 : 0.90*p.NEr;
        x0[1] = data[z].n2 > -0.1 ? data[z].n2 : 0.90*p.NEr;
        x0[2] = data[z].n3 > -0.1 ? data[z].n3 : 0.90*p.NEr;
        x0[3] = data[z].n4 > -0.1 ? data[z].n4 : 0.90*p.NEr;
        x0[4] = data[z].n5 > -0.1 ? data[z].n5 : 0.90*p.NYb;
        x0[5] = data[z].n6 > -0.1 ? data[z].n6 : 0.90*p.NYb;
        
    } else if (z > 0)
    {
        //Use last values as starting values (should be different from zero if the simulation
        //forward iteration was run
        x0[0] = data[z-1].n1 > -0.1 ? data[z-1].n1 : 0.90*p.NEr;
        x0[1] = data[z-1].n2 > -0.1 ? data[z-1].n2 : 0.90*p.NEr;
        x0[2] = data[z-1].n3 > -0.1 ? data[z-1].n3 : 0.90*p.NEr;
        x0[3] = data[z-1].n4 > -0.1 ? data[z-1].n4 : 0.90*p.NEr;
        x0[4] = data[z-1].n5 > -0.1 ? data[z-1].n5 : 0.90*p.NYb;
        x0[5] = data[z-1].n6 > -0.1 ? data[z-1].n6 : 0.90*p.NYb;
     }
    
    //std::vector<double> output {Maths::jac_newton_method(f, x0, h, tol, n_it, p_)}; //old
    std::vector<double> output {Maths::jac_newton_method(f, x0, df, tol, n_it, p_)}; //new
    for (auto i =0; i < output.size(); ++i)
    {
        if (output[i] > (i <= 4 ? p.NEr:p.NYb)) 
        {
            output[i] = i <= 4 ? p.NEr:p.NYb;
            continue;
        }
        
        if (output[i] < -0.1)
        {
            output[i] = 0.0;
        }
    }
        
    data[z].n1 = output[0];
    data[z].n2 = output[1];
    data[z].n3 = output[2];
    data[z].n4 = output[3];
    data[z].n5 = output[4];
    data[z].n6 = output[5];
    data[z].residual = output[6];
    
    return;
}


void Simulation::Result::reset_start() noexcept
{
    //Reset the forwards ASE which should be 0 at the start
    auto ASE_f_start {data[0].PASE_f.begin()};
    auto ASE_f_end {data[0].PASE_f.end()};
    for (ASE_f_start; ASE_f_start != ASE_f_end; ++ASE_f_start)
    {
        ASE_f_start->second = 0.0;
    }
    
    //Reset the forwards pump ASE
    auto p_ASE_f_start {data[0].Pp_ASE_f.begin()};
    auto p_ASE_f_end {data[0].Pp_ASE_f.end()};
    for (p_ASE_f_start; p_ASE_f_start != p_ASE_f_end; ++p_ASE_f_start)
    {
        p_ASE_f_start->second = 0.0;
    }
    
    //Reset forwards pump to its starting value
    data[0].Pp_f = p.Pp0_f;
    
    //Reset signal to its starting value
    data[0].Ps = p.Ps0;
    
    //Reset NaN values if any
    auto ASE_b_start {data[0].PASE_b.begin()};
    auto ASE_b_end {data[0].PASE_b.end()};
    for (ASE_b_start; ASE_b_start != ASE_b_end; ++ASE_b_start)
    {
        if (std::isnan(ASE_b_start->second) || ASE_b_start->second < 0.0) 
            ASE_b_start->second = 0.0;
    }
    
    auto p_ASE_b_start {data[0].Pp_ASE_b.begin()};
    auto p_ASE_b_end {data[0].Pp_ASE_b.end()};
    for (p_ASE_b_start; p_ASE_b_start != p_ASE_b_end; ++p_ASE_b_start)
    {
        if (std::isnan(p_ASE_b_start->second) || p_ASE_b_start->second < 0.0) 
            p_ASE_b_start->second = 0.0;
    }
    
    
    auto backward_pump_st {data[data.size()-1].Pp_b.begin()};
    auto backward_pump_end {data[data.size()-1].Pp_b.end()};
    for (backward_pump_st; backward_pump_st != backward_pump_end; ++backward_pump_st)
    {
        if (std::isnan(backward_pump_st->second) || backward_pump_st->second < 0.0) 
            backward_pump_st->second = 0.0;
    }
    
} //End of reset start


void Simulation::Result::reset_end() noexcept
{
    //Reset the backwards ASE which should be 0 at the end
    auto ASE_b_start {data[data.size()-1].PASE_b.begin()};
    auto ASE_b_end {data[data.size()-1].PASE_b.end()};
    for (ASE_b_start; ASE_b_start != ASE_b_end; ++ASE_b_start)
    {
        ASE_b_start->second = 0.0;
    }
    
    //Reset the backwards pump ASE which should be 0 at the end
    auto p_ASE_b_start {data[data.size()-1].Pp_ASE_b.begin()};
    auto p_ASE_b_end {data[data.size()-1].Pp_ASE_b.end()};
    for (p_ASE_b_start; p_ASE_b_start != p_ASE_b_end; ++p_ASE_b_start)
    {
        p_ASE_b_start->second = 0.0;
    }
    
    
    //Reset backwards pump to its starting value
    data[data.size()-1].Pp_b = p.Pp0_b;
    
    //Reset NaN values if any
    auto ASE_f_start {data[data.size()-1].PASE_f.begin()};
    auto ASE_f_end {data[data.size()-1].PASE_f.end()};
    for (ASE_f_start; ASE_f_start != ASE_f_end; ++ASE_f_start)
    {
        if (std::isnan(ASE_f_start->second) || ASE_f_start->second < 0.0) 
            ASE_f_start->second = 0.0;
    }
    
    auto p_ASE_f_start {data[data.size()-1].Pp_ASE_f.begin()};
    auto p_ASE_f_end {data[data.size()-1].Pp_ASE_f.end()};
    for (p_ASE_f_start; p_ASE_f_start != p_ASE_f_end; ++p_ASE_f_start)
    {
        if (std::isnan(p_ASE_f_start->second) || p_ASE_f_start->second < 0.0) 
            p_ASE_f_start->second = 0.0;
    }
    
    auto signal_st {data[data.size()-1].Ps.begin()}; 
    auto signal_end {data[data.size()-1].Ps.end()};
    for (signal_st; signal_st != signal_end; ++signal_st)
    {
        if (std::isnan(signal_st->second) || signal_st->second < 0.0) 
            signal_st->second = 0.0;
    }
    
    auto forward_pump_st {data[data.size()-1].Pp_f.begin()};
    auto forward_pump_end {data[data.size()-1].Pp_f.end()};
    for (forward_pump_st; forward_pump_st != forward_pump_end; ++forward_pump_st)
    {
        if (std::isnan(forward_pump_st->second) || forward_pump_st->second < 0.0) 
            forward_pump_st->second = 0.0;
    }
    
} //End of reset end

void Simulation::Result::advance_step (u_int const z, double const h, double const tol, u_int n_it) noexcept
{
    find_all_n(z, h, tol, n_it);
    
    if (z < data.size()-1)
    {
        advance_signal  (z); //Advances 1 step the signal
        advance_pump_f  (z); //Advances 1 step back pump
        advance_pump_b  (z); //Advances 1 step forward pump
        advance_ASE_f   (z); //Advances 1 step forward ASE
        advance_ASE_b   (z); //Advances 1 step backward ASE
        //advance_p_ASE_f (z); //Advances 1 step forward ASE
        advance_p_ASE_b (z); //Advances 1 step backward ASE
    }

    return;
}


void Simulation::Result::regress_step (u_int const z, double const h, double const tol, u_int n_it) noexcept
{
    find_all_n(z, h, tol, n_it);
    
    if (z > 0)
    {
        advance_signal  (z, -1); //Advances 1 step the signal
        advance_pump_f  (z, -1); //Advances 1 step back pump
        advance_pump_b  (z, -1); //Advances 1 step forward pump
        advance_ASE_f   (z, -1); //Advances 1 step forward ASE
        advance_ASE_b   (z, -1); //Advances 1 step back ASE
        advance_p_ASE_f (z, -1); //Advances 1 step forward ASE
        advance_p_ASE_b (z, -1); //Advances 1 step back ASE
    }

    return;
}


void Simulation::Result::simulate(bool warn) noexcept
{
    double const NAvg {(p.NYb+p.NEr) / 2.0};
    double const h {1.0e-7 * NAvg};
    double const tol {5.0e-6 * NAvg};
    u_int n_it {1000};

    
    for (int n = 0; n < 3; ++n)
    {
        std::cout<<"===================cycle "<<n<<"===============\n";

        for (int i = 0; i < data.size(); ++i)
        {
            if (i % 20 == 0) std::cout<<"step = "<<i<<'\n';
            advance_step(i, h, tol, n_it);  
        
            if (warn             && (
                data[i].n1 < 0.0 ||
                data[i].n2 < 0.0 ||
                data[i].n3 < 0.0 ||
                data[i].n4 < 0.0 ||
                data[i].n5 < 0.0 ||
                data[i].n6 < 0.0))
            {
                std::cout<<"Warning in simulate: negative value found at step "<<i<<'\n';
                report_step(i, false);
            }
        } //End of fowards iteration
    
        //Reset PASE_b to 0 and Pp_b to Pp0_b
        reset_end();
    
        for (int i = data.size()-1; i >= 0 ; --i)
        {
        
            if (i % 20 == 0) std::cout<<"regress step = "<<i<<'\n';
            regress_step(i, h, tol, n_it);
        
            if (warn             && (
                data[i].n1 < 0.0 ||
                data[i].n2 < 0.0 ||
                data[i].n3 < 0.0 ||
                data[i].n4 < 0.0 ||
                data[i].n5 < 0.0 ||
                data[i].n6 < 0.0))
            {
                std::cout<<"Warning in simulate: negative value found at step "<<i<<'\n';
                report_step(i, true);
            }
        
        } //End of backwards iteration

        //Resets Ps to Ps0, PASE_f to 0, Pp_f to Pp0_f
        reset_start();
        
    } //End of Cycle
}


void Simulation::Result::save_data(std::string const filename, bool dBm_units) noexcept
{
    std::cout<<"Writing file: "<<filename<<'\n';
    std::ofstream file_handle {filename};
    std::string dBm {dBm_units ? "dBm" : "mW"};
    file_handle<<"Length,n1,n2,n3,n4,n5,n6,Ps ("+dBm+"),"
    "Gain,Pp_f (976),Pp_f (1480),Pp_b (976),Pp_b (1480),PASE_f,PASE_b,Pp_ASE_f,Pp_ASE_b\n";
    
    for (auto i = 0; i < data.size(); ++i)
    {
        double const i_d {static_cast<double>(i)};
        
        file_handle<<(p.step_size*i_d*100.0)
                   <<','
                   <<data[i].n1 * 100.0/p.NEr
                   <<','
                   <<data[i].n2 * 100.0/p.NEr
                   <<','
                   <<data[i].n3 * 100.0/p.NEr
                   <<','
                   <<data[i].n4 * 100.0/p.NEr
                   <<','
                   <<data[i].n5 * 100.0/p.NYb
                   <<','
                   <<data[i].n6 * 100.0/p.NYb
                   <<','
                   <<(dBm_units ? Utility::power_to_dBm(data[i].Ps.at(1533000)) : data[i].Ps.at(1533000) * 1000.0)
                   <<','
                   <<Utility::power_to_gain(data[0].Ps.at(1533000), data[i].Ps.at(1533000))
                   <<','
                   <<(data[i].Pp_f.count(976000)  > 0 ? data[i].Pp_f.at(976000) * 1000.0  : 0.0)
                   <<','
                   <<(data[i].Pp_f.count(1480000) > 0 ? data[i].Pp_f.at(1480000) * 1000.0 : 0.0)
                   <<','
                   <<(data[i].Pp_b.count(976000)  > 0 ? data[i].Pp_b.at(976000) * 1000.0  : 0.0)
                   <<','
                   <<(data[i].Pp_b.count(1480000) > 0 ? data[i].Pp_b.at(1480000) * 1000.0 : 0.0)
                   <<','
                   <<(dBm_units ? Utility::power_to_dBm(data[i].PASE_f.at(1533000)) : data[i].PASE_f.at(1533000) * 1000.0)
                   <<','
                   <<(dBm_units ? Utility::power_to_dBm(data[i].PASE_b.at(1533000)) : data[i].PASE_b.at(1533000) * 1000.0)
                   <<','
                   <<(dBm_units ? Utility::power_to_dBm(data[i].Pp_ASE_f.at(980000)) : data[i].Pp_ASE_f.at(980000) * 1000.0)
                   <<','
                   <<(dBm_units ? Utility::power_to_dBm(data[i].Pp_ASE_b.at(980000)) : data[i].Pp_ASE_b.at(980000) * 1000.0)
                   <<'\n';
    }
    
    file_handle.close();
    std::cout<<"Finished writing file\n";
    system("gnuplot plot_script.txt");
    
    return;
}


void Simulation::Result::plot_data (std::string const data_file, std::string const plot_script) noexcept
{
    std::ofstream script_out {plot_script};
    std::regex rx (R"(([a-zA-z0-9_]+)(.txt))"); //We have to add* which means 0 or more
    std::regex dBm(R"(dBm)");
    std::smatch matches {};
    std::smatch dBm_match{};
    std::regex_match(data_file, matches, rx);
    std::ifstream file {data_file};
    std::string val_to_search;
    getline(file,val_to_search);
    bool dBm_units {std::regex_match(val_to_search, dBm_match, dBm)};
    file.close();
    std::cout<<"Plotting data...\n";
    std::string name {matches[1]};
    
    std::string textfile {
    "set terminal pngcairo size 1200, 600\n"
    "set datafile separator \",\"\n"
    "unset arrow\n"
    "set view map\n"
    "set tics front\n"
    "set autoscale y\n"
    "set autoscale x\n"
    "filename = \""+data_file+"\"\n"
    "\n"
    "set output \""+name+"_plot_1.png\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \"Pop. Inversion (%)\"\n"
    "set title \"Population inversion\"\n"
    "plot filename using 1:2 title \'n1 ^4I_{15/2}\' lt rgb \"black\" with lines,\\\n"
    "filename using 1:3 title \'n2 ^4I_{13/2}\' lt rgb \"red\" with lines,\\\n"
    "filename using 1:4 title \'n3 ^4I_{11/2}\' lt rgb \"blue\" with lines,\\\n"
    "filename using 1:5 title \'n4 ^4I_{9/2}\' lt rgb \"green\" with lines,\\\n"
    "filename using 1:6 title \'n5 ^2F_{7/2}\' lt rgb \"orange\" with lines,\\\n"
    "filename using 1:7 title \'n6 ^2F_{5/2}\' lt rgb \"violet\" with lines\n"
    "\n"
    "set autoscale y\n"
    "set output \""+name+"_plot_2.png\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \"Pop. Inversion (%)\"\n"
    "set title \"Population inversion\"\n"
    "plot filename using 1:4 title \'n3 ^4I_{11/2}\' lt rgb \"blue\" with lines,\\\n"
    "filename using 1:5 title \'n4 ^4I_{9/2}\' lt rgb \"green\" with lines,\\\n"
    "\n"
    "set autoscale y\n"
    "set output \""+name+"_plot_3.png\"\n"
    "set title \"Signal at 1533 nm\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \""+(dBm_units ? "Power (dBm)" : "Power (mW)")+"\"\n"
    "plot filename using 1:8 title \'Signal\' with lines\n"
    "\n"
    "set output \""+name+"_plot_4.png\"\n"
    "set title \"Signal gain at 1533 nm\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \"Gain (dB)\"\n"
    "plot filename using 1:9 title \'Signal\' with lines\n"
    "\n"
    "set output \""+name+"_plot_5.png\"\n"
    "set title \"Forward and backwards pump\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \"Pump power (mW)\"\n"
    "plot filename using 1:10 title \'Pp_f 976 nm\' with lines,\\\n"
    "filename using 1:11 title \'Pp_f 1480 nm\' lt rgb \"red\" with lines,\\\n"
    "filename using 1:12 title \'Pp_b 976 nm\' with lines,\\\n"
    "filename using 1:13 title \'Pp_b 1480 nm\' lt rgb \"green\" with lines\n"
    "\n"
    "set output \""+name+"_plot_6.png\"\n"
    "set title \"Forward and backward ASE at 1533 nm\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \""+(dBm_units ? "Power (dBm)" : "Power (mW)")+"\"\n"
    "plot filename using 1:14 title \'PASE_f\' with lines,\\\n"
    "filename using 1:15 title \'PASE_b\' with lines\n"
    "\n"
    "set output \""+name+"_plot_7.png\"\n"
    "set title \"Forward and backward pump ASE at 980 nm\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \""+(dBm_units ? "Power (dBm)" : "Power (mW)")+"\"\n"
    "plot filename using 1:16 title \'PpASE_f\' with lines,\\\n"
    "filename using 1:17 title \'PpASE_b\' with lines\n"};
    script_out<<textfile;
    script_out.close();
    std::string const command {"gnuplot " + plot_script};
    system(command.data());
}


std::vector<std::array<double, 4>> Simulation::find_ratio(Simulation::Init_params const p) noexcept
{
    double const N_total {p.NYb + p.NEr};
    std::vector<std::array<double, 4>> output;
    for (auto i = 1; i < 11; ++i)
    {
        double const i_d {static_cast<double>(i)};
        std::array<double, 4> elem {0.0, 0.0, 0.0, 0.0};
        elem[1] = N_total * i_d/(2.0+i_d-1);
        elem[0] = N_total - elem[1];
        output.emplace_back(elem);
    }
    
    for (auto i = 1; i < 11; ++i)
    {
        double const i_d {static_cast<double>(i)};
        std::array<double, 4> elem {0.0, 0.0, 0.0, 0.0};
        elem[0] = N_total * i_d/(2.0+i_d-1);
        elem[1] = N_total - elem[0];
        output.emplace_back(elem);
    }
    
    for (auto i = 0; i < output.size(); ++i)
    {
        Simulation::Init_params p_ {p};
        p_.NEr = output[i][0];
        p_.NYb = output[i][1];
        p_.recalculate_constants();
        
        Simulation::Result r{p_};
        r.simulate();
        double highest_val {-1.0};
        
        for (auto j = 0; j < r.data.size(); ++j)
        {
            if (r.data[j].Ps.at(1533000) > highest_val) 
            {
                output[i][2] = static_cast<double>(j) * r.p.step_size*100.0;
                highest_val = r.data[j].Ps.at(1533000);
            }
        }
        std::cout<<"===============================highest_val = "<<highest_val<<'\n';
        output[i][3] = Utility::power_to_gain(r.data[0].Ps.at(1533000), highest_val);
    }
    //Save in a file steeam
    std::ofstream file_handle {"Ratios.csv"};
    file_handle<<"NEr,"<<"NYb,"<<"Ratio Er/Yb,"<<"z pos (cm)"<<"Max Gain,"<<","<<"Gain/cm\n";
    for (auto i : output)
    {
        file_handle<<i[0]<<","<<i[1]<<","<<i[2]<<","<<(i[0]/i[1])<<","<<i[3]<<","<<i[3]/i[2]<<'\n';
    }
    
    return output;
}


void Simulation::Result::save_spectral_data(std::string const filename, u_int const step, bool dBm_units) noexcept
{
    std::string name
    {filename[filename.size()-4] == '.' ? filename.substr(0, filename.size()-4) : filename.substr(0, filename.size())};
    std::string extension
    {filename[filename.size()-4] == '.' ? filename.substr(filename.size()-4, 4) : std::string{}};
    
    
    std::string new_name {name + "_step_" + std::to_string(step) + extension};
    std::cout<<"Saving spectral information: "<<new_name<<'\n';
    std::ofstream file_handle {new_name};
    std::string dBm {dBm_units ? "dBm" : "mW"};
    
    file_handle<<"Type, Wavelength (nm),Power ("<<dBm<<")\n";
    auto classify = [] (int k)
    {
        switch (k)
        {
            case 0:
            return "Ps";
            
            case 1:
            return "Pp_f";
            
            case 2:
            return "Pp_b";
            
            case 3:
            return "ASE_f";
            
            case 4:
            return "ASE_b";
            
            case 5:
            return "Pp_ASE_f";
            
            case 6:
            return "Pp_ASE_b";
        }
    };
    
    //0. Ps, 1. Pp_f, 2. Pp_b, 3. PASE_f, 4. PASE_b, 5. Pp_ASE_f, 6. Pp_ASE_b
    using it_type = decltype(std::cbegin(data[step].Ps));
    std::array<std::vector<it_type>, 2> iterators {Utility::return_iterators<it_type>(data[step])};
    
    for (auto i = 0; i < iterators[0].size(); ++i)
    {
        for (auto st_it = iterators[0][i]; st_it != iterators[1][i]; ++st_it)
        {
            file_handle<<classify(i)<<','<<st_it->first<<','
            <<(dBm_units ? Utility::power_to_dBm(st_it->second):st_it->second * 1000.0)<<'\n';
        }
    }
    
    file_handle.close();
    return;
}