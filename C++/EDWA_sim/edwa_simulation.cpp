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
            S.PASE_b.emplace(wl, 0.0);
            S.PASE_f.emplace(wl, 0.0);
        }
    }
}


Simulation::Result::Result(Simulation::Init_params const& initial_state) noexcept
{
    p = initial_state;
    data.reserve(p.steps);
    
    
    //Creating the first values
    data.emplace_back(Simulation::Step{});
        data[0].curr_length = 0.0;      //Length
        data[0].lsEr = p.lsEr;
        data[0].lpEr = p.lpEr;
        data[0].lsYb = p.lsYb;
        data[0].lpYb = p.lpYb;
    
    for (auto i = std::cbegin(p.Pp0_f); i != std::cend(p.Pp0_f); ++i) {
        int const wl = i->first;
        data[0].Pp_f[wl] = 0.0;
    }       

    for (auto i = std::cbegin(p.Pp0_b); i != std::cend(p.Pp0_b); ++i) {
        int const wl = i->first;
        data[0].Pp_b[wl] = 0.0;
    }

    for (auto i = std::cbegin(p.Ps0); i != std::cend(p.Ps0); ++i) {
        int const wl = i->first;
            data[0].Ps[wl] = 0.0;
    }

    for (auto i = 1; i < p.steps; ++i)
    {
        data.emplace_back(Simulation::Step{});
        data[i].curr_length = static_cast<double>(i) * p.step_size;      //Step number
        data[i].lsEr = p.lsEr;
        data[i].lpEr = p.lpEr;
        data[i].lsYb = p.lsYb;
        data[i].lpYb = p.lpYb;
        data[i].Pp_f = data[0].Pp_f;     //Initialize with forwards pump with 0
        data[i].Pp_b = data[0].Pp_b;     //Initialize backwards pump with 0
        data[i].Ps   = data[0].Ps;       //Initialize signal with 0
    }
    
    data[0].Pp_f = p.Pp0_f;            //Initial forwards pump
    data[p.steps-1].Pp_b = p.Pp0_b;    //Initial backwards pump
    data[0].Ps   = p.Ps0;              //Initial signal
    
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
            double const pump_pow2  {st2->second};
            double const flux       {Utility::return_photon_flux(
                wl2,
                pump_pow2,
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
        
        
        if (var == 13)
        {
            //Forward Pump contribution
            auto st6 {data[z].Pp_f.cbegin()}; 
            auto end6 {data[z].Pp_f.cend()};
    
            for (st6; st6 != end6; ++st6)
            {
                int const wl6        {(st6->first)};
                //Validity interval for Er absorption in this range 962 nm to 990 nm
                if (wl6 < 962000) continue;
                if (wl6 > 990000) break;
                double const pump_pow6  {st6->second};
                double const flux       {Utility::return_photon_flux(
                    wl6,
                    pump_pow6,
                    p.A,
                    var_
                )};
        
                result += flux;
            }
            
            
            //Backwards pump
            auto st7 {data[z].Pp_b.cbegin()}; 
            auto end7 {data[z].Pp_b.cend()};
    
            for (st7; st7 != end7; ++st7)
            {
                int const wl7        {(st7->first)};
                //Validity interval for Er absorption in this range 962 nm to 990 nm
                if (wl7 < 962000) continue;
                if (wl7 > 990000) break;
                double const pump_pow7  {st7->second};
                double const flux       {Utility::return_photon_flux(
                    wl7,
                    pump_pow7,
                    p.A,
                    var_
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
    double const lpEr {data[z].lpEr};
    double const lpYb {data[z].lpYb};
    double const tot_loss {Utility::return_conc_loss(lpEr + lpYb)};
    double output     {0.0};

    int const wl {wl_};
    double const power {val};
    double sigEr_13 {0.0};
    if (wl >= 962000 && wl <= 990000) sigEr_13 = Utility::return_cross_section(wl, Cross_sec::Er_abs);
    double const sigYb_56 {Utility::return_cross_section(wl, Cross_sec::Yb_abs)};
    double const overlap {p.overlap.at(wl_)};
    output = -power * overlap * (sigEr_13 * n1 + sigYb_56 * n5) - (tot_loss) * power;

    return output;
}


double Simulation::Result::dPp_b_ind (u_int const z, int const wl_, double const val) noexcept
{
    double const n1   {data[z].n1};
    double const n5   {data[z].n5};
    double const lpEr {data[z].lpEr};
    double const lpYb {data[z].lpYb};
    double const tot_loss {Utility::return_conc_loss(lpEr + lpYb)};
    double output     {0.0};
    
    int const wl {wl_};
    double const power {val};
    double sigEr_13 {0.0};
    if (wl >= 962000 && wl <= 990000) sigEr_13 = Utility::return_cross_section(wl, Cross_sec::Er_abs);
    double const sigYb_56 {Utility::return_cross_section(wl, Cross_sec::Yb_abs)};
    double const overlap {p.overlap.at(wl_)};
    output = power * overlap * (sigEr_13 * n1 + sigYb_56 * n5) + (tot_loss) * power;

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



//Advances signal one step forwards or backwards
void Simulation::Result::advance_signal(u_int const z, int const sign) noexcept
{
    double const h {p.step_size};
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
        data[z+sign].Ps.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        if (data[z+sign].Ps.at(wl) < 0.0) data[z+sign].Ps.at(wl) = 0.0;
    }
    
    return;
}


//Advances forwad pump one step forwards or backwards
void Simulation::Result::advance_pump_f(u_int const z, int const sign) noexcept
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
            data[z+sign].Pp_f.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            if (data[z+sign].Pp_f.at(wl) < 0.0) data[z+sign].Pp_f.at(wl) = 0.0;
        } else {
            //1480 nm pump is stored with signals
            double const P0 {st_it->second};
            double const k1 {h * sign_d * dPs_ind(z, wl, P0)};
            double const k2 {h * sign_d * dPs_ind(z, wl, (P0 + k1/2.0))};
            double const k3 {h * sign_d * dPs_ind(z, wl, (P0 + k2/2.0))};
            double const k4 {h * sign_d * dPs_ind(z, wl, (P0 + k3))};
            data[z+sign].Pp_f.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            if (data[z+sign].Pp_f.at(wl) < 0.0) data[z+sign].Pp_f.at(wl) = 0.0;
        }
    }
    
    return;
}


//Advances backwards pump one step forward or backwards
void Simulation::Result::advance_pump_b(u_int const z, int const sign) noexcept
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
            if (data[z+sign*1].Pp_b.at(wl) < 0.0) data[z+sign].Pp_b.at(wl) = 0.0;
        } else {
            //1480 nm pump is stored with signals
            double const P0 {st_it->second};
            double const k1 {h * sign_d * -1.0 * dPs_ind(z, wl, P0)};
            double const k2 {h * sign_d * -1.0 * dPs_ind(z, wl, (P0 + k1/2.0))};
            double const k3 {h * sign_d * -1.0 * dPs_ind(z, wl, (P0 + k2/2.0))};
            double const k4 {h * sign_d * -1.0 * dPs_ind(z, wl, (P0 + k3))};
            data[z+sign].Pp_b.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
            if (data[z+sign].Pp_b.at(wl) < 0.0) data[z+sign].Pp_b.at(wl) = 0.0;
        }
    }
    
    return;
}


//Advances forward ASE one step forward or backwards
void Simulation::Result::advance_ASE_f(u_int const z, int const sign) noexcept
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
        data[z+sign].PASE_f.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        
        if (data[z+sign].PASE_f.at(wl) < 0.0) data[z+sign].PASE_f.at(wl) = 0.0;
    }
    
    return;
}


//Advances backward ASE one step forward or backwards
void Simulation::Result::advance_ASE_b(u_int const z, int const sign) noexcept
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
        data[z+sign].PASE_b.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        if (data[z+sign].PASE_b.at(wl) < 0.0) data[z+sign].PASE_b.at(wl) = 0.0;
    }
    
    return;
}



double Simulation::Result::Eq_1(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n1 {vars[0]};
    double const n2 {vars[1]};
    double const n3 {vars[2]};
    double const n4 {p.NEr - vars[0] - vars[1] - vars[2]};
    double const n6 {vars[3]};
    
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
    double const n4 {p.NEr - vars[0] - vars[1] - vars[2]};
    
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
    double const n4 {p.NEr - vars[0] - vars[1] - vars[2]};
    double const n6 {vars[3]};
    
    double const W13 {p.W13};
    double const t1  {W13 * n1 - p.A32 * n3};
    double const t2  {p.A43 * n4 - 2.0 * p.Cup * pow(n3, 2.0)};
    double const t3  {p.Ccr * n1 * n6};
    
    return t1 + t2 + t3;
}


double Simulation::Result::Eq_4(std::vector<double> const& vars, Data_p const& p) noexcept
{
    double const n1 {vars[0]};
    double const n5 {p.NYb - vars[3]};
    double const n6 {vars[3]};
    
    double const W56 {p.W56};
    double const W65 {p.W65};
    double const t1  {-W56 * n5 + (p.A65 + W65) * n6};
    double const t2  {p.Ccr * n1 * n6};
    
    return t1 + t2;
}


void Simulation::Result::find_all_n(u_int const z, double const h, double const tol, u_int const n_it) noexcept
{
    
    
    Maths::Iter_params const it_p {h, tol, n_it}; 
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
    
    data[z].W12 = p_.W12;
    data[z].W13 = p_.W13;
    data[z].W21 = p_.W21;
    data[z].W65 = p_.W65;
    data[z].W56 = p_.W56;
    data[z].A21 = p_.A21;
    data[z].A32 = p_.A32;
    data[z].A43 = p_.A43;
    data[z].A65 = p_.A65;
    
    
    
    
    using f_ptr = double (*)(std::vector<double> const&, Data_p const&);
    std::vector<f_ptr> f {Simulation::Result::Eq_1, 
                          Simulation::Result::Eq_2,
                          Simulation::Result::Eq_3,
                          Simulation::Result::Eq_4};
    
    //std::array<double, 4> init_g {0.20, 0.70, 0.10, 0.45};                      
    std::vector<double> x0(4, 0.0);
    
    if (first_run) {
        x0[0] = 0.017 * p.NEr;  
        x0[1] = 0.980 * p.NEr;
        x0[2] = 0.002 * p.NEr;
        x0[3] = 0.010 * p.NYb;
        first_run = false;
    } else if (z == 0) {
        x0[0] = data[z].n1;  
        x0[1] = data[z].n2;
        x0[2] = data[z].n3;
        x0[3] = data[z].n6;
    } else {
        x0[0] = data[z-1].n1;  
        x0[1] = data[z-1].n2;
        x0[2] = data[z-1].n3;
        x0[3] = data[z-1].n6;
    }
    
    try {
    
        std::vector<double> output {Maths::newton_method(f, x0, it_p, p_)};

        data[z].n1 = output[0];
        data[z].n2 = output[1];
        data[z].n3 = output[2];
        data[z].n4 = p.NEr - output[0] - output[1] - output[2];
        data[z].n5 = p.NYb - output[3];
        data[z].n6 = output[3];
    
    } 
    catch (Maths::Error& e) 
    {
        //Throw any other type of error that is not dimension shrinking
        if (e.et != Maths::Error_type::zero_value_determinant) {
            throw(e);
        }
        
        //Will attempt to recalculate assuming that the populations n5 and n6 do not change (the 
        //Jacobian's last row is 0) by redifining the Jacobian to exclude the last row
        std::vector<f_ptr> f2 {Simulation::Result::Eq_1, 
                          Simulation::Result::Eq_2,
                          Simulation::Result::Eq_3};
        std::vector<double> x02(3, 0.0);
        if (first_run) {
        x02[0] = 0.017 * p.NEr;  
        x02[1] = 0.980 * p.NEr;
        x02[2] = 0.002 * p.NEr;
        first_run = false;
    } else if (z == 0) {
        x02[0] = data[z].n1;  
        x02[1] = data[z].n2;
        x02[2] = data[z].n3;
    } else {
        x02[0] = data[z-1].n1;  
        x02[1] = data[z-1].n2;
        x02[2] = data[z-1].n3;
    }
    
    std::vector<double> output2 {Maths::newton_method(f, x02, it_p, p_)};
    
    data[z].n1 = output2[0];
    data[z].n2 = output2[1];
    data[z].n3 = output2[2];
    data[z].n4 = p.NEr - output2[0] - output2[1] - output2[2];
    data[z].n5 = data[z-1].n5;
    data[z].n6 = data[z-1].n6;
    
    }//End of catch
    
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
    
} //End of reset start


void Simulation::Result::reset_end() noexcept
{
    //Reset the backwards ASE which should be 0 at the end
    auto last_elem {data.size() - 1};
    auto ASE_b_start {data[last_elem].PASE_b.begin()};
    auto ASE_b_end {data[last_elem].PASE_b.end()};
    for (ASE_b_start; ASE_b_start != ASE_b_end; ++ASE_b_start)
    {
        ASE_b_start->second = 0.0;
    }
    
    
    //Reset backwards pump to its starting value
    data[last_elem].Pp_b = p.Pp0_b;
    
    //Reset NaN values if any
    auto ASE_f_start {data[last_elem].PASE_f.begin()};
    auto ASE_f_end {data[last_elem].PASE_f.end()};
    for (ASE_f_start; ASE_f_start != ASE_f_end; ++ASE_f_start)
    {
        if (std::isnan(ASE_f_start->second) || ASE_f_start->second < 0.0) 
            ASE_f_start->second = 0.0;
    }
    
    
    auto signal_st {data[last_elem].Ps.begin()}; 
    auto signal_end {data[last_elem].Ps.end()};
    for (signal_st; signal_st != signal_end; ++signal_st)
    {
        if (std::isnan(signal_st->second) || signal_st->second < 0.0) 
            signal_st->second = 0.0;
    }
    
    
} //End of reset end

void Simulation::Result::advance_step (u_int const z, 
                                       double const h, 
                                       double const tol, 
                                       u_int n_it, 
                                       bool enable_ASE) noexcept
{
    find_all_n(z, h, tol, n_it);
    
    auto last_elem {data.size() - 1};
    
    if (z < last_elem)
    {
        advance_signal  (z); //Advances 1 step the signal
        advance_pump_f  (z); //Advances 1 step back pump
        advance_pump_b  (z); //Advances 1 step forward pump
        if (enable_ASE)
        {
            advance_ASE_f   (z); //Advances 1 step forward ASE
            advance_ASE_b   (z); //Advances 1 step backward ASE
        }
    }

    return;
}


void Simulation::Result::regress_step (u_int const z, 
                                       double const h, 
                                       double const tol, 
                                       u_int n_it, 
                                       bool enable_ASE) noexcept
{
    find_all_n(z, h, tol, n_it);
    
    if (z > 0)
    {
        advance_signal  (z, -1); //Advances 1 step the signal
        advance_pump_f  (z, -1); //Advances 1 step back pump
        advance_pump_b  (z, -1); //Advances 1 step forward pump
        if (enable_ASE)
        {
            advance_ASE_f   (z, -1); //Advances 1 step forward ASE
            advance_ASE_b   (z, -1); //Advances 1 step backward ASE
        }
    }

    return;
}


void Simulation::Result::simulate(float& report, bool const warn, bool const enable_ASE) noexcept
{
    
    double const NAvg {(p.NYb + p.NEr) / 2.0};
    double const h {1.0e-7};
    double const tol {1.0e-7 * NAvg};
    u_int n_it {1000};
    int n_loops {3};
    float curr_tot {0.0};
    float total {static_cast<float>(n_loops) * static_cast<float>(data.size()) * 2.0f};
    
    for (int n = 0; n < n_loops; ++n)
    {
        for (int i = 0; i < data.size(); ++i)
        {
            curr_tot += 1.0;
            report = curr_tot / total;
            advance_step(i, h, tol, n_it, enable_ASE);  

        } //End of fowards iteration
    
        //Reset PASE_b to 0 and Pp_b to Pp0_b
        reset_end();
    
        for (int i = data.size()-1; i >= 0 ; --i)
        {
            curr_tot += 1.0;
            report = curr_tot / total;
            regress_step(i, h, tol, n_it, enable_ASE);
        
        } //End of backwards iteration

        //Resets Ps to Ps0, PASE_f to 0, Pp_f to Pp0_f
        reset_start();
      
    } //End of Cycle
    
    for (auto i = 0; i < data.size(); ++i) {
        log(i);
    }
} //End of simulate





void Simulation::Result::save_data(std::string_view filename, bool const dBm_units, int const s_wl, int const p_wl_1, int const p_wl_2) noexcept
{
    std::cout<<"Entered save data\n";
    std::cout<<"s_wl = "<<s_wl<<'\n';
    std::cout<<"p_wl_1 = "<<p_wl_1<<'\n';
    std::cout<<"p_wl_2 = "<<p_wl_2<<'\n';
    boost::filesystem::path curr_path{boost::filesystem::current_path()};
    boost::filesystem::path folder {"output"};
    curr_path /= folder;
    if (!boost::filesystem::exists(curr_path)) boost::filesystem::create_directories(curr_path);
    std::regex rx ("([a-zA-Z0-9_]+)(.*)(\w*)");
    std::smatch s;
    std::string const filename_str {filename.data()};
    std::regex_match(filename_str, s, rx);
    std::string file_name {s[1]};
    curr_path /= boost::filesystem::path(file_name);
    if (!boost::filesystem::exists(curr_path)) boost::filesystem::create_directories(curr_path);
    file_name = file_name + std::string(".csv");
    curr_path /= boost::filesystem::path(file_name);
    
    std::ofstream file_handle {curr_path.string()};
    
    std::string dBm {dBm_units ? "dBm" : "mW"};
    std::string s_wl_s {std::to_string(s_wl/1000)};
    std::string p_wl_1_s {std::to_string(p_wl_1/1000)};
    std::string p_wl_2_s {std::to_string(p_wl_2/1000)};
    file_handle<<"Length,n1,n2,n3,n4,n5,n6,Ps @ "+s_wl_s+" ("+dBm+"),"
    "Gain,Pp_f ("+p_wl_1_s+"),Pp_f ("+p_wl_2_s+"),Pp_b ("+p_wl_1_s+"),Pp_b ("+p_wl_2_s+"),PASE_f,PASE_b,"
    "W12,W13,W21,W65,W56,A21,A32,A43,A65,NEr,NYb,Cup,Ccr\n";

    for (auto i = 0; i < data.size(); ++i)
    {
        double const i_d {static_cast<double>(i)};
        
        auto Pp_f_wl_1 = p_wl_1 >= 1420000 ? data[i].Ps : data[i].Pp_f;
        auto Pp_f_wl_2 = p_wl_2 >= 1420000 ? data[i].Ps : data[i].Pp_f;
        auto Pp_b_wl_1 = p_wl_1 >= 1420000 ? data[i].Ps : data[i].Pp_b;
        auto Pp_b_wl_2 = p_wl_2 >= 1420000 ? data[i].Ps : data[i].Pp_b;
        
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
                   <<(dBm_units ? Utility::power_to_dBm(data[i].Ps.at(s_wl)) : data[i].Ps.at(s_wl) * 1000.0)
                   <<','
                   <<Utility::power_to_gain(data[0].Ps.at(s_wl), data[i].Ps.at(s_wl))
                   <<','
                   <<(Pp_f_wl_1.count(p_wl_1) > 0 ? Pp_f_wl_1.at(p_wl_1) * 1000.0 : 0.0)
                   <<','
                   <<(Pp_f_wl_2.count(p_wl_2) > 0 ? Pp_f_wl_2.at(p_wl_2) * 1000.0 : 0.0)
                   <<','
                   <<(Pp_b_wl_1.count(p_wl_1) > 0 ? Pp_b_wl_1.at(p_wl_1) * 1000.0 : 0.0)
                   <<','
                   <<(Pp_b_wl_2.count(p_wl_2) > 0 ? Pp_b_wl_2.at(p_wl_2) * 1000.0 : 0.0)
                   <<','
                   <<(dBm_units ? Utility::power_to_dBm(data[i].PASE_f.at(s_wl)) : data[i].PASE_f.at(s_wl) * 1000.0)
                   <<','
                   <<(dBm_units ? Utility::power_to_dBm(data[i].PASE_b.at(s_wl)) : data[i].PASE_b.at(s_wl) * 1000.0)
                   <<','
                   <<data[i].W12
                   <<','
                   <<data[i].W13
                   <<','
                   <<data[i].W21
                   <<','
                   <<data[i].W65
                   <<','
                   <<data[i].W56
                   <<','
                   <<data[i].A21
                   <<','
                   <<data[i].A32
                   <<','
                   <<data[i].A43
                   <<','
                   <<data[i].A65
                   <<','
                   <<p.NEr
                   <<','
                   <<p.NYb
                   <<','
                   <<p.Cup
                   <<','
                   <<p.Ccr
                   <<'\n';
    }

    file_handle.close();
    std::cout<<"Exited save_data\n";
    return;
}


void Simulation::Result::plot_data (std::string_view data_file_, int const s_wl, int const p_wl_1, int const p_wl_2) noexcept
{
    std::string s_wl_s {std::to_string(s_wl / 1000)};
    std::string p_wl_1_s {std::to_string(p_wl_1 / 1000)};
    std::string p_wl_2_s {std::to_string(p_wl_2 / 1000)};
    boost::filesystem::path curr_path{boost::filesystem::current_path()};
    boost::filesystem::path folder {"output"};
    curr_path /= folder;
    if (!boost::filesystem::exists(curr_path)) boost::filesystem::create_directories(curr_path);
    std::regex rx (R"(([a-zA-Z0-9_]+)?(\.\w+)?)");
    std::regex dBm(R"((dBm))");
    std::smatch matches {};
    std::smatch dBm_match{};
    std::string const data_file_str {data_file_.data()};
    std::regex_match(data_file_str, matches, rx);
    std::string data_file {matches[1]};
    curr_path /= boost::filesystem::path(data_file);
    boost::filesystem::path curr_path_copy {curr_path};
    data_file = data_file + std::string(".csv");
    curr_path_copy /= boost::filesystem::path(data_file);
    data_file = curr_path_copy.string();
    std::string filename {matches[1]};
    if (!boost::filesystem::exists(curr_path)) boost::filesystem::create_directories(curr_path);
    curr_path /= boost::filesystem::path(filename);
    std::string plot_script {curr_path.string() + std::string("_script") + std::string(".txt")};
    std::ofstream script_out {plot_script};
    std::ifstream file {data_file};
    std::string val_to_search;
    getline(file, val_to_search);
    bool dBm_units {std::regex_search(val_to_search, dBm_match, dBm)};
    file.close();
    std::string name {matches[1]};
    curr_path_copy = curr_path_copy.parent_path() / boost::filesystem::path(name);
    name = curr_path_copy.string();
    
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
    "set title \"Signal at "+s_wl_s+" nm\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \""+(dBm_units ? "Power (dBm)" : "Power (mW)")+"\"\n"
    "plot filename using 1:8 title \'Signal\' with lines\n"
    "\n"
    "set output \""+name+"_plot_4.png\"\n"
    "set title \"Signal gain at "+s_wl_s+" nm\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \"Gain (dB)\"\n"
    "plot filename using 1:9 title \'Signal\' with lines\n"
    "\n"
    "set output \""+name+"_plot_5.png\"\n"
    "set title \"Forward and backwards pump\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \"Pump power (mW)\"\n"
    "plot filename using 1:10 title \'Pp_f "+p_wl_1_s+" nm\' with lines,\\\n"
    "filename using 1:11 title \'Pp_f "+p_wl_2_s+" nm\' lt rgb \"red\" with lines,\\\n"
    "filename using 1:12 title \'Pp_b "+p_wl_1_s+" nm\' with lines,\\\n"
    "filename using 1:13 title \'Pp_b "+p_wl_2_s+"\' lt rgb \"green\" with lines\n"
    "\n"};
    
    script_out<<textfile;
    script_out.close();
    std::string const command {"gnuplot " + plot_script};
    system(command.data());
}


void Simulation::Result::save_spectral_data(std::string_view filename, u_int const step, bool const dBm_units) noexcept
{
    std::string name
    {filename[filename.size()-4] == '.' ? filename.substr(0, filename.size()-4) : filename.substr(0, filename.size())};
    std::string extension
    {filename[filename.size()-4] == '.' ? filename.substr(filename.size()-4, 4) : std::string{".csv"}};
    
    
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
            file_handle<<classify(i)<<','<<st_it->first/1000<<','
            <<(dBm_units ? Utility::power_to_dBm(st_it->second):st_it->second * 1000.0)<<'\n';
        }
    }
    
    file_handle.close();
    return;
}


