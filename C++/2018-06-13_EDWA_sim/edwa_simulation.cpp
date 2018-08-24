#include "edwa_simulation.h"

//Using abbreviations


//Constructor based on initial state
void Simulation::Result::initialize_ASE() noexcept
{
    for (Simulation::Step& S : data) 
    {
        int step_size {100000 / (p.n_ASE-1)};
        for (int j = 0; j < p.n_ASE; ++j)
        {
            int const wl {1500000 + j*step_size};
            S.PASE_b.emplace(wl, 0.0);
            S.PASE_f.emplace(wl, 0.0);
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
        
        //Bakwards ASE contribution
        auto st3 {data[z].PASE_b.cbegin()}; //Iterate over forward pump in case there are more than one
        auto end3 {data[z].PASE_b.cend()};
    
        for (st3; st3 != end3; ++st3)
        {
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
            auto st4 {data[z].Pp_f.cbegin()}; 
            auto end4 {data[z].Pp_f.cend()};
    
            for (st4; st4 != end4; ++st4)
            {
                int const wl4        {(st4->first)};
                double const pump_pow1  {st4->second};
                double const flux       {Utility::return_photon_flux(
                    wl4,
                    pump_pow1,
                    p.A,
                    Cross_sec::Yb_abs
                )};
        
                result += flux;
            }
            
            auto st5 {data[z].Pp_b.cbegin()}; 
            auto end5 {data[z].Pp_b.cend()};
    
            for (st5; st5 != end5; ++st5)
            {
                int const wl5        {(st5->first)};
                double const pump_pow1  {st5->second};
                double const flux       {Utility::return_photon_flux(
                    wl5,
                    pump_pow1,
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
    double const lpEr {data[z].lpEr};
    double const lpYb {data[z].lpYb};
    double output     {0.0};
    
    int const wl {wl_};
    double const power {val};
    double const sigEr_13 {Utility::return_cross_section(wl, Cross_sec::Er_abs)};
    double const sigYb_56 {Utility::return_cross_section(wl, Cross_sec::Yb_abs)};
    double const overlap {p.overlap.at(wl_)};
    output = -power * overlap * (sigEr_13 * n1 + sigYb_56 * n5) - (lpEr + lpYb) * power;
    
    return output;
}


double Simulation::Result::dPp_b_ind (u_int const z, int const wl_, double const val) noexcept
{
    double const n1   {data[z].n1};
    double const n5   {data[z].n5};
    double const lpEr {data[z].lpEr};
    double const lpYb {data[z].lpYb};
    double output     {0.0};
    
    int const wl {wl_};
    double const power {val};
    double const sigEr_13 {Utility::return_cross_section(wl, Cross_sec::Er_abs)};
    double const sigYb_56 {Utility::return_cross_section(wl, Cross_sec::Yb_abs)};
    double const overlap {p.overlap.at(wl_)};
    output = power * overlap * (sigEr_13 * n1 + sigYb_56 * n5) + (lpEr + lpYb) * power;
    
    return output;                                          
}


double Simulation::Result::dPs_ind (u_int const z, int const wl_, double const val) noexcept
{
    double const n1   {data[z].n1};
    double const n2   {data[z].n2};
    double const lsEr {data[z].lsEr};
    double const lsYb {data[z].lsYb};
    double output     {0.0};
    
    int const wl {wl_};
    double const power {val};
    double const sigEr_21 {Utility::return_cross_section(wl, Cross_sec::Er_emi)};
    double const sigEr_12 {Utility::return_cross_section(wl, Cross_sec::Er_abs)};
    double const overlap {p.overlap.at(wl_)};
    output = power * overlap * (sigEr_21 * n2 - sigEr_12 * n1) - (lsEr + lsYb) * power;
    return output;                                          
}


double Simulation::Result::dPASE_f_ind (u_int const z, int const wl_, double const val) noexcept
{
    double const n1 {data[z].n1};
    double const n2 {data[z].n2};
    double const lsEr {data[z].lsEr};
    double const lsYb {data[z].lsYb};
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
    double const c {(lsEr + lsYb) * power};
    output = a + b - c;
    
    return output;
}


double Simulation::Result::dPASE_b_ind (u_int const z, int const wl_, double const val) noexcept
{
    double const n1 {data[z].n1};
    double const n2 {data[z].n2};
    double const lsEr {data[z].lsEr};
    double const lsYb {data[z].lsYb};
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
    double const c {(lsEr + lsYb) * power};
    output = -a - b + c;

    return output;
}


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
void Simulation::Result::advance_signal(u_int const z, double sign) noexcept
{
    double const h {p.step_size};
    //double const h {p.step_size/101};
    auto st_it {data[z].Ps.begin()}; //start iterator
    auto ed_it {data[z].Ps.end()};   //End iterator
    int sign_ {static_cast<int>(sign)};
    
    //Starting values of Runge - Kutta
    for (st_it; st_it != ed_it; ++st_it)
    {
        int const wl {st_it->first};
        double P0 {st_it->second};
        //for (auto i = 0; i < 101; ++i)
        //{
        //std::cout<<"P0 = "<<P0<<", wl = "<<wl<<'\n';
        double const k1 {h * sign * dPs_ind(z, wl, P0)};
        double const k2 {h * sign * dPs_ind(z, wl, (P0 + k1/2.0))};
        double const k3 {h * sign * dPs_ind(z, wl, (P0 + k2/2.0))};
        double const k4 {h * sign * dPs_ind(z, wl, (P0 + k3))};
        //P0 += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        //std::cout<<"k1 = "<<k1<<", k2 = "<<k2<<", k3 = "<<k3<<", k4 = "<<k4<<'\n';
        //}
        
        data[z+sign_*1].Ps.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        //data[z+sign_*1].Ps.at(wl) = P0;
        
        //std::cout<<"h = "<<h<<'\n';
    }
    
    return;
}


//Advances forwad pump one step forwards or backwards
void Simulation::Result::advance_pump_f(u_int const z, double sign) noexcept
{
    double const h {p.step_size};
    auto st_it {data[z].Pp_f.begin()}; //start iterator
    auto ed_it {data[z].Pp_f.end()};   //End iterator
    int sign_ {static_cast<int>(sign)};
    
    //Starting values of Runge - Kutta
    for (st_it; st_it != ed_it; ++st_it)
    {
        int const wl {st_it->first};
        double const P0 {st_it->second};
        double const k1 {h * sign * dPp_f_ind(z, wl, P0)};
        double const k2 {h * sign * dPp_f_ind(z, wl, (P0 + k1/2.0))};
        double const k3 {h * sign * dPp_f_ind(z, wl, (P0 + k2/2.0))};
        double const k4 {h * sign * dPp_f_ind(z, wl, (P0 + k3))};
        data[z+sign_*1].Pp_f.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    
    return;
}


//Advances backwards pump one step forward or backwards
void Simulation::Result::advance_pump_b(u_int const z, double sign) noexcept
{
    double const h {p.step_size};
    auto st_it {data[z].Pp_b.begin()}; //start iterator
    auto ed_it {data[z].Pp_b.end()};   //End iterator
    int sign_ {static_cast<int>(sign)};
    
    //Starting values of Runge - Kutta
    for (st_it; st_it != ed_it; ++st_it)
    {
        int const wl {st_it->first};
        double const P0 {st_it->second};
        double const k1 {h * sign * dPp_b_ind(z, wl, P0)};
        double const k2 {h * sign * dPp_b_ind(z, wl, (P0 + k1/2.0))};
        double const k3 {h * sign * dPp_b_ind(z, wl, (P0 + k2/2.0))};
        double const k4 {h * sign * dPp_b_ind(z, wl, (P0 + k3))};
        data[z+sign_*1].Pp_b.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    
    return;
}


//Advances forward ASE one step forward or backwards
void Simulation::Result::advance_ASE_f(u_int const z, double sign) noexcept
{
    double const h {p.step_size};
    auto st_it {data[z].PASE_f.begin()}; //start iterator
    auto ed_it {data[z].PASE_f.end()};   //End iterator
    int sign_ {static_cast<int>(sign)};
    //Starting values of Runge - Kutta
    for (st_it; st_it != ed_it; ++st_it)
    {
        int const wl {st_it->first};
        double const P0 {st_it->second};
        double const k1 {h * sign * dPASE_f_ind(z, wl, P0)};
        double const k2 {h * sign * dPASE_f_ind(z, wl, (P0 + k1/2.0))};
        double const k3 {h * sign * dPASE_f_ind(z, wl, (P0 + k2/2.0))};
        double const k4 {h * sign * dPASE_f_ind(z, wl, (P0 + k3))};
        data[z+sign_*1].PASE_f.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }
    
    return;
}


//Advances backward ASE one step forward or backwards
void Simulation::Result::advance_ASE_b(u_int const z, double sign) noexcept
{
    double const h {p.step_size};
    auto st_it {data[z].PASE_b.begin()}; //start iterator
    auto ed_it {data[z].PASE_b.end()};   //End iterator
    int sign_ {static_cast<int>(sign)};
    
    //Starting values of Runge - Kutta
    for (st_it; st_it != ed_it; ++st_it)
    {
        int const wl {st_it->first};
        double const P0 {st_it->second};
        double const k1 {h * sign * dPASE_b_ind(z, wl, P0)};
        double const k2 {h * sign * dPASE_b_ind(z, wl, (P0 + k1/2.0))};
        double const k3 {h * sign * dPASE_b_ind(z, wl, (P0 + k2/2.0))};
        double const k4 {h * sign * dPASE_b_ind(z, wl, (P0 + k3))};
        
        data[z+sign_*1].PASE_b.at(wl) = P0 + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
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
    } else if (z == 0 &&
        data[z].n1 > 0.0 &&
        data[z].n2 > 0.0 &&
        data[z].n3 > 0.0 &&
        data[z].n4 > 0.0 &&
        data[z].n5 > 0.0 &&
        data[z].n6 > 0.0)
    {
        //If the simulation has been run previously then use curr vals as initial guess
        x0[0] = data[z].n1;
        x0[1] = data[z].n2;
        x0[2] = data[z].n3;
        x0[3] = data[z].n4;
        x0[4] = data[z].n5;
        x0[5] = data[z].n6;
        
    } else if (z > 0)
    {
        //Use last values as starting values (should be different from zero if the simulation
        //forward iteration was run
        x0[0] = data[z-1].n1;
        x0[1] = data[z-1].n2;
        x0[2] = data[z-1].n3;
        x0[3] = data[z-1].n4;
        x0[4] = data[z-1].n5;
        x0[5] = data[z-1].n6;
    }
    
    std::vector<double> output {Maths::jac_newton_method(f, x0, h, tol, n_it, p_)};
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
    
    //Reset forwards pump to its starting value
    data[0].Pp_f = p.Pp0_f;
    
    //Reset signal to its starting value
    data[0].Ps = p.Ps0;
}


void Simulation::Result::reset_end() noexcept
{
    //Reset the backwards ASE which should be 0 at the end
    auto ASE_b_start {data[data.size()-1].PASE_b.begin()};
    auto ASE_b_end {data[data.size()-1].PASE_b.end()};
    for (ASE_b_start; ASE_b_start != ASE_b_end; ++ASE_b_start)
    {
        ASE_b_start->second = 0.0;
    }
    
    //Reset backwards pump to its starting value
    data[data.size()-1].Pp_b = p.Pp0_b;
    
}

void Simulation::Result::advance_step (u_int const z, double const h, double const tol, u_int n_it) noexcept
{
    find_all_n(z, h, tol, n_it);
    if (z < data.size()-1)
    {
        advance_signal (z); //Advances 1 step the signal
        advance_pump_f (z); //Advances 1 step back pump
        advance_pump_b (z); //Advances 1 step forward pump
        advance_ASE_f  (z); //Advances 1 step forward ASE
        advance_ASE_b  (z); //Advances 1 step back ASE
    }

    return;
}


void Simulation::Result::regress_step (u_int const z, double const h, double const tol, u_int n_it) noexcept
{
    find_all_n(z, h, tol, n_it);
    if (z > 0)
    {
        advance_signal (z, -1.0); //Advances 1 step the signal
        advance_pump_f (z, -1.0); //Advances 1 step back pump
        advance_pump_b (z, -1.0); //Advances 1 step forward pump
        advance_ASE_f  (z, -1.0); //Advances 1 step forward ASE
        advance_ASE_b  (z, -1.0); //Advances 1 step back ASE
    }

    return;
}


Matrix<double> Simulation::Result::delta_it(Simulation::Step const& a, Simulation::Step const& b) noexcept
{
    Matrix<double> delta{5, 1};//Ps, Pp_f, Pp_b, PASE_f, PASE_b
    //Gives an array with all start iterators and end iterators
    std::array<std::vector<decltype(a.Ps.begin())>, 2> a_iterators {Utility::return_iterators<decltype(a.Ps.begin())>(a)};
    std::array<std::vector<decltype(b.Ps.begin())>, 2> b_iterators {Utility::return_iterators<decltype(b.Ps.begin())>(b)};
    double a_val {0.0};
    Matrix<double> a_vals{5,1};
    double b_val {0.0};
    Matrix<double> b_vals{5,1};
    
    for (auto j = 0; j < a_iterators[0].size(); ++j)
    {
        for (auto it = a_iterators[0][j]; it != a_iterators[1][j]; ++it)
        {
            a_val += it->second;
        }
        
        a_vals[j] = a_val;
        a_val = 0.0;
        
        
        for (auto it = b_iterators[0][j]; it != b_iterators[1][j]; ++it)
        {
            b_val += it->second;
        }
        
        b_vals[j] = b_val;
        b_val = 0;
        
        delta[j] = Maths::abs_val(a_vals[j] - b_vals[j]);
    }
    
    return delta;
} //End of delta_it


void Simulation::Result::simulate() noexcept
{
    double const NAvg {(p.NYb+p.NEr) / 2.0};
    double const h {1.0e-7 * NAvg};
    double const tol {1.0e-6 * NAvg};
    u_int n_it {1000};
    
    for (int n = 0; n < 2; ++n)
    {
        std::cout<<"===================cycle "<<n<<"===============\n";
        for (int i = 0; i < data.size(); ++i)
        {
            if (i % 20 == 0) std::cout<<"step = "<<i<<'\n';
            advance_step(i, h, tol, n_it);
        
            if (i == 0)
            {
                start_step = data[i];
            }
        
            if (i == data.size()-1)
            {
                end_step = data[i];
            }
        
            if (data[i].n1 < 0.0 ||
                data[i].n2 < 0.0 ||
                data[i].n3 < 0.0 ||
                data[i].n4 < 0.0 ||
                data[i].n5 < 0.0 ||
                data[i].n6 < 0.0)
            {
                std::cout<<"Warning in simulate: negative value found at step "<<i<<'\n';
                report_step(i, true);
            }
        } //End of fowards iteration
    
        //Reset PASE_b to 0 and Pp_b to Pp0_b
        reset_end();
    
        for (int i = data.size()-1; i >= 0 ; --i)
        {
        
            if (i % 20 == 0) std::cout<<"regress step = "<<i<<'\n';
            regress_step(i, h, tol, n_it);
        
            if (i == 0)
            {
                start_step = data[i];
            }
        
            if (i == data.size()-1)
            {
                end_step = data[i];
            }
        
            if (data[i].n1 < 0.0 ||
                data[i].n2 < 0.0 ||
                data[i].n3 < 0.0 ||
                data[i].n4 < 0.0 ||
                data[i].n5 < 0.0 ||
                data[i].n6 < 0.0)
            {
                std::cout<<"Warning in simulate: negative value found at step "<<i<<'\n';
                report_step(i, true);
            }
        
        } //End of backwards iteration

        //Resets Ps to Ps0, PASE_f to 0, Pp_f to Pp0_f
        reset_start();
    }//End of cycle
}


void Simulation::Result::save_data(std::string const filename, bool dBm_units) noexcept
{
    std::cout<<"Writing file: "<<filename<<'\n';
    std::ofstream file_handle {filename};
    
    file_handle<<"Length,n1,n2,n3,n4,n5,n6,Ps,Gain, Pp_f,Pp_b,PASE_f,PASE_b\n";
    
    for (auto i = 0; i < data.size(); ++i)
    {
        double const i_d {static_cast<double>(i)};
        
        file_handle<<(p.step_size*i_d*100.0)
                   <<','
                   <<data[i].n1
                   <<','
                   <<data[i].n2
                   <<','
                   <<data[i].n3
                   <<','
                   <<data[i].n4
                   <<','
                   <<data[i].n5
                   <<','
                   <<data[i].n6
                   <<','
                   <<(dBm_units ? Utility::power_to_dBm(data[i].Ps.at(1533000)) : data[i].Ps.at(1533000) * 100.0)
                   <<','
                   <<Utility::power_to_gain(data[0].Ps.at(1533000), data[i].Ps.at(1533000) * 100.0)
                   <<','
                   <<data[i].Pp_f.at(976000) * 100.0
                   <<','
                   <<data[i].Pp_b.at(976000) * 100.0
                   <<','
                   <<(dBm_units ? Utility::power_to_dBm(data[i].PASE_f.at(1533000)) : data[i].PASE_f.at(1533000) * 100.0)
                   <<','
                   <<(dBm_units ? Utility::power_to_dBm(data[i].PASE_b.at(1533000)) : data[i].PASE_b.at(1533000) * 100.0)
                   <<'\n';
    }
    
    file_handle.close();
    std::cout<<"Finished writing file\n";
    system("gnuplot plot_script.txt");
    
    return;
}


void Simulation::Result::plot_data (std::string const data_file, std::string const plot_script, bool dBm_units) noexcept
{
    std::ofstream script_out {plot_script};
    std::regex rx (R"(([a-zA-z0-9_]+)(.txt))"); //We have to add* which means 0 or more
    std::smatch matches {};
    std::regex_match(data_file, matches, rx);
    std:cout<<"Plotting data...\n";
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
    "set ylabel \"Pop. density (ions/m^3)\"\n"
    "set title \"Population densities\"\n"
    "plot filename using 1:2 title \'n1 ^4I_{15/2}\' with lines,\\\n"
    "filename using 1:3 title \'n2 ^4I_{15/2}\' with lines,\\\n"
    "filename using 1:4 title \'n3 ^4I_{13/2}\' lt rgb \"black\" with lines,\\\n"
    "filename using 1:5 title \'n4 ^4I_{9/2}\' with lines,\\\n"
    "filename using 1:6 title \'n5 ^2F_{7/2}\' with lines,\\\n"
    "filename using 1:7 title \'n6 ^2F_{5/2}\' with lines\n"
    "\n"
    "set yrange[0:1e23]\n"
    "set output \""+name+"_plot_2.png\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \"Pop. density (ions/m^3)\"\n"
    "set title \"Population densities\"\n"
    "plot filename using 1:2 title \'n1 ^4I_{15/2}\' with lines,\\\n"
    "filename using 1:3 title \'n2 ^4I_{15/2}\' with lines,\\\n"
    "filename using 1:4 title \'n3 ^4I_{13/2}\' lt rgb \"black\" with lines,\\\n"
    "filename using 1:5 title \'n4 ^4I_{9/2}\' with lines,\\\n"
    "filename using 1:6 title \'n5 ^2F_{7/2}\' with lines,\\\n"
    "filename using 1:7 title \'n6 ^2F_{5/2}\' with lines\n"
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
    "set title \"Forward and backwards pump change at 976 nm\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \"Pump power (mW)\"\n"
    "plot filename using 1:10 title \'Pp_f\' with lines,\\\n"
    "filename using 1:11 title \'Pp_b\' with lines\n"
    "\n"
    "set output \""+name+"_plot_6.png\"\n"
    "set title \"Forward and backward ASE at 1533 nm\"\n"
    "set xlabel \"Length (cm)\"\n"
    "set ylabel \""+(dBm_units ? "Power (dBm)" : "Power (mW)")+"\"\n"
    "plot filename using 1:12 title \'PASE_f\' with lines,\\\n"
    "filename using 1:13 title \'PASE_b\' with lines\n"};

    script_out<<textfile;
    script_out.close();
    std::string const command {"gnuplot " + plot_script};
    system(command.data());
}


double Simulation::Result::find_critical_point() noexcept
{
    
}
