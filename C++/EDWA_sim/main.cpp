#include "headers.h"
//Pending
//Fix scattering by erbium and ytterbium
//Double check up conversion and cross relaxation
//Fix plotting of 1480 nm pump
//Add overlap columns to the GUI
//Add function to calculate overlap


//Helper struct for Parallel_idx
struct Trio
{
    float first {0.0f};
    bool second {false};
    std::string third {std::string()};
};

//Enables to check which index is available to represent in the ImGui::ProgressBar
struct Parallel_idx
{
    std::mutex m;
    std::vector<Trio> idx;
    
    //Initializes by checking the number of threads and creating a vector of that size
    Parallel_idx()
    {
        std::size_t ts {std::thread::hardware_concurrency()};
        for (auto i = 0; i < ts; ++i)
        idx.emplace_back(Trio());
    }
    
    //Deleted copy constructors to ensure only references are used
    Parallel_idx(Parallel_idx const&) = delete;
    Parallel_idx& operator = (Parallel_idx const&) = delete;
    
    //Returns an optional containing size_t in case there are no available indexes to assign
    std::optional<std::size_t> assign_idx(std::string_view s)
    {
        //Locks the vector to prevent other threads from accessing it
        std::lock_guard<std::mutex> lk(m);
        std::optional<std::size_t> output;
        
        for (std::size_t i = 0;  i < idx.size(); ++i)
        {
            if (idx[i].second == false)
            {
                //If the index is assigned, reset the value, set to true and save the profile name
                idx[i].first  = 0.0f;
                idx[i].second = true;
                idx[i].third  = s;
                *output = i;
                break;
            }
            
        }
        
        return output;
    }

    //Set to false the corresponding variable to signal the index is no longer in use
    void release_idx(std::size_t i)
    {
        std::lock_guard<std::mutex> lk(m);
        idx[i].second = false;
        return;
    }
};


//This structure helps convert the input from the gui to Init_params required by the simulation
struct Merge_helper
{
    std::array<double, 30> temp_doubles;
    std::array<bool, 30> temp_bools;
    wl_map temp_signals;
    wl_map temp_f_pump;
    wl_map temp_b_pump;
    
    //Initializes Merge_helper with the default values of Init_params
    void init_temp(Simulation::Init_params& p)
    {
        
        temp_doubles[0] = p.width * 1.0e6;
        temp_doubles[1] = p.height * 1.0e6;
        temp_doubles[2] = p.l * 1.0e2;
        temp_doubles[3] = p.NEr;
        temp_doubles[4] = p.NYb;
        temp_doubles[5] = (1.0/p.A21) * 1.0e3;
        temp_doubles[6] = (1.0/p.A32) * 1.0e3;
        temp_doubles[7] = (1.0/p.A43) * 1.0e3;
        temp_doubles[8] = (1.0/p.A65) * 1.0e3;
        temp_doubles[9] = p.Cup;
        temp_doubles[10] = p.Ccr;
        temp_doubles[11] = p.lpEr;
        temp_doubles[12] = p.lsEr;
        temp_doubles[13] = p.lpYb;
        temp_doubles[14] = p.lsYb;
        
        for (auto i = 15; i < temp_doubles.size(); ++i)
        {
            temp_doubles[i] = 0.0;
        }
        
        //Default values for plotting
        temp_doubles[21] = 1533.0;
        temp_doubles[22] = 976.0;
        temp_doubles[23] = 1480.0;
        
        temp_bools[0] = true;
        temp_bools[1] = true;
        for (auto i = 2; i < temp_bools.size(); ++i)
        {
            temp_bools[i] = false;
        }
        
        return;
    }
    
    //Coonverts and merges the values into the Init_params format
    void merge_with_sim(Simulation::Init_params& p) const
    {
        p.width = temp_doubles[0] * 1.0e-6;
        p.height = temp_doubles[1] * 1.0e-6;
        p.set_length((temp_doubles[2] * 1.0e-2));
        p.NEr = temp_doubles[3];
        p.NYb = temp_doubles[4];
        p.A21 = 1.0/(temp_doubles[5] * 1.0e-3);
        p.A32 = 1.0/(temp_doubles[6] * 1.0e-3);
        p.A43 = 1.0/(temp_doubles[7] * 1.0e-3);
        p.A65 = 1.0/(temp_doubles[8] * 1.0e-3);
        p.Ps0 = temp_signals;
        p.Pp0_f = temp_f_pump;
        p.Pp0_b = temp_b_pump;
        p.Cup = temp_doubles[9];
        p.Ccr = temp_doubles[10];
        p.lpEr = temp_doubles[11];
        p.lsEr = temp_doubles[12]; 
        p.lpYb = temp_doubles[13];
        p.lsYb = temp_doubles[14];
        
        return;
    }
}; //End of Merge_helper struct


//Creates a vector of Merge_helpers with the desired number of elements
std::vector<Merge_helper> create_helpers(Simulation::Init_params& p, std::size_t const n)
{
    std::vector<Merge_helper> output;
    
    for (auto i = 0; i < n; ++i)
    {
        Merge_helper m;
        m.init_temp(p);
        output.emplace_back(m);
    }
    
    return output;
}

//Function from ImGuiDemo.cpp to show help markers
static void ShowHelpMarker(const char* desc)
{
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered())
    {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}


//This function will fill the columns in the gui with wavelength, power and type, reading from a
//wl_map
void log_map(wl_map& m, std::string s)
{
    for (auto& i : m)
    {
        std::stringstream ss;
        ss<<(i.first / 1000);
        ImGui::Text(ss.str().c_str()); ImGui::NextColumn();
        ss.str(std::string());
        ss.clear();
        ss<<(i.second * 1000);
        ImGui::Text(ss.str().c_str()); ImGui::NextColumn();
        ImGui::Text(s.c_str()); ImGui::NextColumn();
    }
}


//Will place as mane progress bars as threads available for concurrency
void multi_progress_bar(Parallel_idx& f, int curr_idx_)
{
    std::size_t threads(f.idx.size());
    ImGui::PushItemWidth(ImGui::GetWindowContentRegionWidth());
    for (auto i = 0; i < threads; ++i)
    {
        char inner_buf[32];
        sprintf(inner_buf, "Running %s %.2f/%.0f%s", f.idx[i].third.c_str(), f.idx[i].first*100.0f, 100.0f, "%");
        ImGui::ProgressBar(f.idx[i].first, ImVec2(0.0f,0.0f), inner_buf);
    }
    
    ImGui::PopItemWidth();
    
    return;
}


//Struct with overloaded () operator, used to help launch threads with the simulation parameters
struct Launch_sim
{
    Simulation::Init_params p_; 
    Merge_helper m; 
    std::string name; 
    std::size_t n_threads;
    
    Launch_sim(Simulation::Init_params const&p, Merge_helper const& m, 
               std::string const& name, std::size_t const n_threads):
               p_{p}, m{m}, name{name}, n_threads{n_threads} {}
    
    void operator() (std::atomic<std::size_t>& av_threads, Parallel_idx& progress)
    {
        --av_threads;
        std::size_t curr_t;
        std::optional<std::size_t> opt_t {progress.assign_idx(name)};
        curr_t = *opt_t;
        Simulation::Result r{p_};
        
        wl_map temp_signals;
        wl_map temp_f_pump;
        wl_map temp_b_pump;
        progress.idx[curr_t].first = 0.0f;
        r.simulate(progress.idx[curr_t].first, false, m.temp_bools[1]);
        int const s_wl   {static_cast<int>(m.temp_doubles[21]) * 1000};
        int const p_wl_1 {static_cast<int>(m.temp_doubles[22]) * 1000};
        int const p_wl_2 {static_cast<int>(m.temp_doubles[23]) * 1000};
        r.save_data(name, m.temp_bools[3], s_wl, p_wl_1, p_wl_2);
        r.plot_data(name, s_wl, p_wl_1, p_wl_2);
        progress.release_idx(curr_t);
        ++av_threads;
        return;
    }
    
    //Default copy constructor and copy assignment to ensure copies of the variables are distributed
    //to the threads
    Launch_sim& operator = (Launch_sim const&) = default;
    Launch_sim(Launch_sim const&) = default;
};

//Helper enum for the combine_wl_map
enum class S_type{Signal, F_pump, B_pump};
//Adds a signal or pump value to the current profile

void combine_wl_map(std::vector<Merge_helper>& m_vec, int& vec_idx, int& curr_idx, S_type t)
{
    if (t == S_type::B_pump)
    {
        int temp_wl {static_cast<int>(m_vec[curr_idx].temp_doubles[19] * 1000.0)};
        if (!(m_vec[curr_idx].temp_doubles[19] < 850.0 || m_vec[curr_idx].temp_doubles[19] > 1150.0))
        {
            double temp_pow {m_vec[curr_idx].temp_doubles[20] / 1000.0};
            m_vec[vec_idx].temp_b_pump.emplace(std::pair(temp_wl, temp_pow));
        }
    } else if (t == S_type::Signal)
    {
        int temp_wl {static_cast<int>(m_vec[curr_idx].temp_doubles[15] * 1000.0)};
        if (!(m_vec[curr_idx].temp_doubles[15] < 1450.0 || m_vec[curr_idx].temp_doubles[15] > 1600.0))
        {
            double temp_pow {m_vec[curr_idx].temp_doubles[16] / 1000.0};
            m_vec[vec_idx].temp_signals.emplace(std::pair(temp_wl, temp_pow));
        }
    } else 
    {
        int temp_wl {static_cast<int>(m_vec[curr_idx].temp_doubles[17] * 1000.0)};
        if (!(m_vec[curr_idx].temp_doubles[17] < 850.0 || m_vec[curr_idx].temp_doubles[17] > 1150.0))
        {
            double temp_pow {m_vec[curr_idx].temp_doubles[18] / 1000.0};
            m_vec[vec_idx].temp_f_pump.emplace(std::pair(temp_wl, temp_pow));
        }
    }
    
    return;
 }

 
int main(int argc, char **argv)
{
    /*
    //delete afterwards
    auto print_vec = [] (std::vector<double> const& v) {
        auto count = 0;
        for (auto i : v)
        {
            std::cout<<"output["<<count<<"] = "<<i<<'\n';
            ++count;
        }
    
    };
    
    auto return_data_pack = [] (Simulation::Result const& r, std::size_t const z) {
        Simulation::Result::Data_pack pack0;
        pack0.W12 = r.calculate_W(z, 12);
        pack0.W13 = r.calculate_W(z, 13);
        pack0.W21 = r.calculate_W(z, 21);
        pack0.W65 = r.calculate_W(z, 65);
        pack0.W56 = r.calculate_W(z, 56);
        pack0.Cup = r.p.Cup;
        pack0.Ccr = r.p.Ccr;
        pack0.NEr = r.p.NEr;
        pack0.NYb = r.p.NYb;
        pack0.A21 = r.p.A21;
        pack0.A32 = r.p.A32;
        pack0.A43 = r.p.A43;
        pack0.A65 = r.p.A65;
        
        return pack0;
    };
    
    auto integrate_step_data = [] (Simulation::Result& r, 
                                   std::size_t const z, 
                                   std::vector<double> const& output) 
    {
        r.data[z].n1 = output[0];
        r.data[z].n2 = output[1];
        r.data[z].n3 = output[2];
        r.data[z].n6 = output[3];
        r.data[z].n4 = r.p.NEr - output[0] - output[1] - output[2];
        r.data[z].n5 = r.p.NYb - output[3];
    };
    
    auto advance_everything = [] (Simulation::Result& r,
                                  std::size_t const z,
                                  int const sign)
    {
        r.advance_signal(z, sign);
        r.advance_pump_f(z, sign);
        r.advance_pump_b(z, sign);
        //r.advance_ASE_f(z, sign);
        //r.advance_ASE_b(z, sign);
    };
    
    auto return_x0 = [] (Simulation::Result const& r, std::size_t const z) {
        std::vector<double> output {r.data[z].n1,
                                    r.data[z].n2,
                                    r.data[z].n3,
                                    r.data[z].n6
        };
        
        return output;
    };
    
    auto iterate = [&return_data_pack, 
                    &return_x0, 
                    &print_vec,
                    &integrate_step_data,
                    &advance_everything] 
                    (Simulation::Result& r, 
                    std::size_t const& max_it, 
                    Maths::Iter_params<double> const& it_p, std::string const& filename,
                    int sign = 1) 
    {
        
        std::ofstream file_handle {filename};
        std::vector<double (*)(std::vector<double> const&, Simulation::Result::Data_pack const&)> f 
                          {Simulation::Result::Eq_1, 
                          Simulation::Result::Eq_2,
                          Simulation::Result::Eq_3,
                          Simulation::Result::Eq_4};
        
        Simulation::Result::Data_p pack0 {return_data_pack(r, 0)};
        std::vector<double> x0 {0.015*r.p.NEr, 0.98*r.p.NEr, 0.001*r.p.NEr, 0.01*r.p.NYb};    
        std::vector<double> output {Maths::newton_method(f, x0, it_p, pack0)};
    
        integrate_step_data(r, 0, output); //integrates the results from output to the step
        advance_everything(r, 0, 1); //Advances the pump, signal, ASE, etc
    
        file_handle<<"NEr = "<<r.p.NEr<<'\n'<<"NYb = "<<r.p.NYb<<'\n'<<"Length = "<<r.p.l
                    <<'\n'<<"Steps = "<<r.p.steps<<'\n';
        file_handle<<'\n'<<"step,curr_length,n1,n2,n3,n4,n5,n6,Ps(1533),Pp_f(976),Pp_b(976),PASE_f(1533),PASE_b(1533)\n";
            
        auto back_pump = r.p.Pp0_b;
        auto forward_pump = r.p.Pp0_f;
        auto signal = r.p.Ps0;
        
        for (auto j = 0; j < 3; ++j)
        {
        //Forward advance
        for (auto i = 1; i < max_it; ++i) {
            
            auto pack1 {return_data_pack(r, i)};
            auto x1 {return_x0(r, i-1)}; //Uses values of the given step
    
            std::vector<double> output1 {Maths::newton_method(f, x1, it_p, pack1)};

            integrate_step_data(r, i, output1); //integrates the results from output to the step
            if (i == max_it-1) break;
            advance_everything(r, i, sign); //Advances the pump, signal, ASE, etc
        }
      
        r.data[max_it-1].Pp_b = back_pump;
        //Backwards advance
        for (auto i = max_it - 1; i != 0; --i) {
            auto pack1 {return_data_pack(r, i)};
            auto x1 {return_x0(r, i-1)}; //Uses values of the given step
    
            std::vector<double> output1 {Maths::newton_method(f, x1, it_p, pack1)};

            integrate_step_data(r, i, output1); //integrates the results from output to the step
            if (i == 0) break;
            advance_everything(r, i, -1 * sign); //Advances the pump, signal, ASE, etc
        }
        
        r.data[0].Pp_f = forward_pump;
        r.data[0].Ps = signal;
        
        //Forward advance
        for (auto i = 1; i < max_it; ++i) {
            
            auto pack1 {return_data_pack(r, i)};
            auto x1 {return_x0(r, i-1)}; //Uses values of the given step
    
            std::vector<double> output1 {Maths::newton_method(f, x1, it_p, pack1)};

            integrate_step_data(r, i, output1); //integrates the results from output to the step
            if (i == max_it-1) break;
            advance_everything(r, i, sign); //Advances the pump, signal, ASE, etc
        }
        
        r.data[max_it-1].Pp_b = back_pump;
        
        
        } //End of cycle loop


        for (auto i = 0; i < max_it; ++i)
        {
                        auto step   = r.data[i].curr_step;
            auto curr_l = static_cast<double>(r.data[i].curr_step) * r.p.step_size;
            auto n1     = r.data[i].n1 * 100.0 / r.p.NEr;
            auto n2     = r.data[i].n2 * 100.0 / r.p.NEr;
            auto n3     = r.data[i].n3 * 100.0 / r.p.NEr;
            auto n4     = r.data[i].n4 * 100.0 / r.p.NEr;
            auto n5     = r.data[i].n5 * 100.0 / r.p.NYb;
            auto n6     = r.data[i].n6 * 100.0 / r.p.NYb;
            auto Ps     = r.data[i].Ps[1533000];
            auto Pp_f   = r.data[i].Pp_f[976000];
            auto Pp_b   = r.data[i].Pp_b[976000];
            auto PASE_f = r.data[i].PASE_f[1533000];
            auto PASE_b = r.data[i].PASE_b[1533000];
            file_handle<<step<<','<<curr_l<<','<<n1<<','<<n2<<','<<n3<<','<<n4<<','<<n5<<','<<n6
                        <<','<<Ps<<','<<Pp_f<<','<<Pp_b<<','<<PASE_f<<','<<PASE_b<<'\n';
        }
        
        file_handle.close();
    }; //End of iterate
    */
    /*
    //aquiaqui
    Simulation::Init_params p{};
    p.set_steps(1000);
    p.set_length(0.05);
    p.NEr = 4.45e26;
    p.NYb = 4.45e26;
    p.recalculate_constants();
    p.Ps0[1533000] = 1e-6;
    p.Pp0_f[976000] = 0.8;
    p.Pp0_b[976000] = 0.8;
    Simulation::Result r{p};
    float dummy_float{0.0};
    
    double const h = 1.0e-7;
    double const tol = (p.NEr + p.NYb)/2.0 * 1.0e-7;
    //double const tol = 1.0e20;
    std::size_t const n_it = 1000;
    Maths::Iter_params it_p {h, tol, n_it};
    
    //iterate(r, p.steps, it_p, "testout.csv", 1);
   
    
    r.simulate(dummy_float, false, true);
    r.save_data("test_data", false, 1533000, 976000, 1480000);
    r.plot_data("test_data", 1533000, 976000, 1480000);
    r.save_spectral_data("spectral_data", r.data.size()-1, false);
    */

    //dearimgui GUI
    std::size_t n_threads {std::thread::hardware_concurrency()}; //Number of threads
    std::atomic<std::size_t> available_threads {n_threads}; //Threads not in use
    //Support 12 different profiles
    char const * buf[] {"profile_0", "profile_1", "profile_2", "profile_3", "profile_4",
                              "profile_5", "profile_6", "profile_7", "profile_8", "profile_9",
                              "profile_10", "profile_11"};
    Simulation::Init_params p{};
    bool run_all {false};
    bool copy_all {true};
    bool monitor_launched {false};
    bool ready {false};
    int p_steps {301};
    
    //Create size 12 Merge_helper (same size as buf)
    std::vector<Merge_helper> m_vec {create_helpers(p, 12)};
    int curr_idx {0};
    int it_idx{0};
    Merge_helper m = m_vec[0];
    Parallel_idx progress {Parallel_idx()};

    //Creates an 800 x 800 window in which the gui window is going to be displayed
    sf::RenderWindow window(sf::VideoMode(800, 800), "");
    window.setVerticalSyncEnabled(true);
    ImGui::SFML::Init(window);
 
    sf::Color bgColor; //Default colour scheme
    char windowTitle[255] = "EDWA simulation window";
    ImVec2 button_size(100, 50); //Size for some of the ImGui::Button used

    window.setTitle(windowTitle);
    window.resetGLStates(); // call it if you only draw ImGui. Otherwise not needed.
    sf::Clock deltaClock;
    
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            ImGui::SFML::ProcessEvent(event);
 
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }
 
        ImGui::SFML::Update(window, deltaClock.restart());
        
        ImGui::Begin("EDWA Model"); // begin window
    
        if (ImGui::CollapsingHeader("Starting parameters"))
        {
 
            ImGui::Combo("Selected parameters", &curr_idx, buf, m_vec.size());

            
            if (ImGui::TreeNode("Waveguide Dimensions"))
            {
                ImGui::SameLine(); ShowHelpMarker("test help");
                ImGui::PushItemWidth(50);
                ImGui::InputDouble("Width (um)", &m_vec[curr_idx].temp_doubles[0], 0.0f, 0.0f, "%.2f");
                ImGui::InputDouble("Height (um)", &m_vec[curr_idx].temp_doubles[1], 0.0f, 0.0f, "%.2f");
                ImGui::InputDouble("Length (cm)", &m_vec[curr_idx].temp_doubles[2], 0.0f, 0.0f, "%.2f");
                ImGui::PopItemWidth();
                
                if (ImGui::Button("Apply to all profiles##0"))
                {
                    for (auto i = 0; i < m_vec.size(); ++i)
                    {
                        if (i != curr_idx)
                        {
                            m_vec[i].temp_doubles[0] = m_vec[curr_idx].temp_doubles[0];
                            m_vec[i].temp_doubles[1] = m_vec[curr_idx].temp_doubles[1];
                            m_vec[i].temp_doubles[2] = m_vec[curr_idx].temp_doubles[2];
                        }
                    }
                }
                
                ImGui::TreePop();
                ImGui::Separator();
            }
            
            if (ImGui::TreeNode("Rare earth content"))
            {
                ImGui::SameLine(); ShowHelpMarker("test help");
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("Erbium (ions/m^3)", &m_vec[curr_idx].temp_doubles[3], 0.0f, 0.0f, "%e");
                ImGui::InputDouble("Ytterbium (ions/m^3)", &m_vec[curr_idx].temp_doubles[4], 0.0f, 0.0f, "%e");
                ImGui::PopItemWidth();
                
                if (ImGui::Button("Apply to all profiles##1"))
                {
                    for (auto i = 0; i < m_vec.size(); ++i)
                    {
                        if (i != curr_idx)
                        {
                            m_vec[i].temp_doubles[3] = m_vec[curr_idx].temp_doubles[3];
                            m_vec[i].temp_doubles[4] = m_vec[curr_idx].temp_doubles[4];
                        }
                    }
                }
                
                ImGui::TreePop();
                ImGui::Separator();
            }
            
            if (ImGui::TreeNode("Lifetime information"))
            {
                ImGui::SameLine(); ShowHelpMarker("test help");
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("t21 (ms)", &m_vec[curr_idx].temp_doubles[5], 0.0f, 0.0f, "%e");
                ImGui::InputDouble("t32 (ms)", &m_vec[curr_idx].temp_doubles[6], 0.0f, 0.0f, "%e");
                ImGui::InputDouble("t43 (ms)", &m_vec[curr_idx].temp_doubles[7], 0.0f, 0.0f, "%e");
                ImGui::InputDouble("t65 (ms)", &m_vec[curr_idx].temp_doubles[8], 0.0f, 0.0f, "%e");
                ImGui::PopItemWidth();
                
                if (ImGui::Button("Apply to all profiles##2"))
                {
                    for (auto i = 0; i < m_vec.size(); ++i)
                    {
                        if (i != curr_idx)
                        {
                            m_vec[i].temp_doubles[5] = m_vec[curr_idx].temp_doubles[5];
                            m_vec[i].temp_doubles[6] = m_vec[curr_idx].temp_doubles[6];
                            m_vec[i].temp_doubles[7] = m_vec[curr_idx].temp_doubles[7];
                            m_vec[i].temp_doubles[8] = m_vec[curr_idx].temp_doubles[8];
                        }
                    }
                }
                
                ImGui::TreePop();
                ImGui::Separator();
            }
            
            if (ImGui::TreeNode("Other information"))
            {
                ImGui::SameLine(); ShowHelpMarker("test help");
                ImGui::PushItemWidth(100);
                ImGui::InputInt("Steps", &p_steps);
                ImGui::InputDouble("Cup (m^3/s)", &m_vec[curr_idx].temp_doubles[9], 0.0f, 0.0f, "%e");
                ImGui::InputDouble("Ccr (m^3/s)", &m_vec[curr_idx].temp_doubles[10], 0.0f, 0.0f, "%e");
                ImGui::NewLine();
                ImGui::Text("Erbium concentration scattering");
                ImGui::InputDouble("Pump scattering (dB/m)", &m_vec[curr_idx].temp_doubles[11], 0.0f, 0.0f, "%.2f");
                ImGui::InputDouble("Signal scattering (dB/m)", &m_vec[curr_idx].temp_doubles[12], 0.0f, 0.0f, "%.2f");
                ImGui::Text("Ytterbium concentration scattering");
                ImGui::InputDouble("Pump scattering (dB/m)##2", &m_vec[curr_idx].temp_doubles[13], 0.0f, 0.0f, "%.2f");
                ImGui::InputDouble("Signal scattering (dB/m)##2", &m_vec[curr_idx].temp_doubles[14], 0.0f, 0.0f, "%.2f");
                ImGui::PopItemWidth();
                
                if (ImGui::Button("Apply to all profiles##3"))
                {
                    for (auto i = 0; i < m_vec.size(); ++i)
                    {
                        if (i != curr_idx)
                        {
                            m_vec[i].temp_doubles[9] = m_vec[curr_idx].temp_doubles[9];
                            m_vec[i].temp_doubles[10] = m_vec[curr_idx].temp_doubles[10];
                            m_vec[i].temp_doubles[11] = m_vec[curr_idx].temp_doubles[11];
                            m_vec[i].temp_doubles[12] = m_vec[curr_idx].temp_doubles[12];
                            m_vec[i].temp_doubles[13] = m_vec[curr_idx].temp_doubles[13];
                            m_vec[i].temp_doubles[14] = m_vec[curr_idx].temp_doubles[14];
                        }
                    }
                }

                ImGui::TreePop();
                ImGui::Separator();
            }
            
            if (ImGui::TreeNode("Signals & Pump"))
            {
                ImGui::Text("Notes");
                ImGui::SameLine();
                ShowHelpMarker("Signal allowed ranges: 1450 nm to 1600 nm\nPump allowed ranges: 850 nm to 1150 nm\nFor 1480 nm pumping, add it as a signal\nThe 1480 pump direction will be the same as the signal");
                ImGui::NewLine();
                ImGui::Text("Wavelength (nm)"); ImGui::SameLine(150); ImGui::Text("Power (mW)");
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w1", &m_vec[curr_idx].temp_doubles[15], 0.0f, 0.0f, "%.2f");
                ImGui::PopItemWidth();
                ImGui::SameLine(150);
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w2", &m_vec[curr_idx].temp_doubles[16], 0.0f, 0.0f, "%e");
                ImGui::PopItemWidth();
                ImGui::SameLine();
                if (ImGui::Button("Add signal"))
                {
                    if (copy_all)
                    {
                        for (auto i = 0; i < m_vec.size(); ++i)
                        {
                            combine_wl_map(m_vec, i, curr_idx, S_type::Signal);
                        }
                    } else {
                        combine_wl_map(m_vec, curr_idx, curr_idx, S_type::Signal);
                    }
                }
                ImGui::SameLine();
                if (ImGui::Button("Clear all##1"))
                {
                    if (copy_all)
                    {
                        for (auto i = 0; i < m_vec.size(); ++i)
                        {
                            m_vec[i].temp_signals.clear();
                        }
                    } else {
                        m_vec[curr_idx].temp_signals.clear();
                    }
                }
                
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w3", &m_vec[curr_idx].temp_doubles[17], 0.0f, 0.0f, "%.2f");
                ImGui::PopItemWidth();
                ImGui::SameLine(150);
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w4", &m_vec[curr_idx].temp_doubles[18], 0.0f, 0.0f, "%e");
                ImGui::PopItemWidth();
                ImGui::SameLine();
                
                if (ImGui::Button("Add Forward pump"))
                {
                    if (copy_all)
                    {
                        for (auto i = 0; i < m_vec.size(); ++i)
                        {
                            combine_wl_map(m_vec, i, curr_idx, S_type::F_pump);
                        }
                    } else {
                        combine_wl_map(m_vec, curr_idx, curr_idx, S_type::F_pump);
                    }
                }
                
                ImGui::SameLine();
                
                if (ImGui::Button("Clear all##2"))
                {
                    if (copy_all)
                    {
                        for (auto i = 0; i < m_vec.size(); ++i)
                        {
                            m_vec[i].temp_f_pump.clear();
                        }
                    } else {
                        m_vec[curr_idx].temp_f_pump.clear();
                    }
                }
                
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w5", &m_vec[curr_idx].temp_doubles[19], 0.0f, 0.0f, "%.2f");
                ImGui::PopItemWidth();
                ImGui::SameLine(150);
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w6", &m_vec[curr_idx].temp_doubles[20], 0.0f, 0.0f, "%e");
                ImGui::PopItemWidth();
                ImGui::SameLine();
                
                if (ImGui::Button("Add backward pump"))
                {
                    if (copy_all)
                    {
                        for (auto i = 0; i < m_vec.size(); ++i)
                        {
                            combine_wl_map(m_vec, i, curr_idx, S_type::B_pump);
                        }
                    } else {
                        combine_wl_map(m_vec, curr_idx, curr_idx, S_type::B_pump);
                    }
                }
                
                ImGui::SameLine();
                
                if (ImGui::Button("Clear all##3"))
                {
                    if (copy_all)
                    {
                        for (auto i = 0; i < m_vec.size(); ++i)
                        {
                            m_vec[i].temp_b_pump.clear();
                        }
                    } else {
                        m_vec[curr_idx].temp_b_pump.clear();
                    }
                    
                }
                
                ImGui::Checkbox("Apply to all profiles##4", &copy_all);
                ImGui::Separator();
                ImGui::Columns(3, "##test columns");
                ImGui::Text("Wavelength (nm)"); ImGui::NextColumn();
                ImGui::Text("Power (mW)"); ImGui::NextColumn();
                ImGui::Text("Type"); ImGui::NextColumn();
                ImGui::Separator();
                
                if (m_vec[curr_idx].temp_signals.size() != 0)
                {
                    log_map(m_vec[curr_idx].temp_signals, "Signal");
                }
                
                if (m_vec[curr_idx].temp_f_pump.size() != 0)
                {
                    log_map(m_vec[curr_idx].temp_f_pump, "Forward pump");
                }
                
                if (m_vec[curr_idx].temp_b_pump.size() != 0)
                {
                    log_map(m_vec[curr_idx].temp_b_pump, "Backward pump");
                }
                
                ImGui::Columns(1);
                
                ImGui::TreePop();
                ImGui::Separator();
            }
            
            ImGui::NewLine();
            ImGui::Checkbox("Recalculate Cup and Ccr", &m_vec[curr_idx].temp_bools[0]);
            ImGui::SameLine(); ShowHelpMarker("Will use the default equations to calculate the values. I.E. will overwrite user provided values\n");
            ImGui::NewLine();
        } //End of starting parameters
        
        
        if (ImGui::CollapsingHeader("Simulation"))
        {

            ImGui::Text("Plotting options");
            ImGui::PushItemWidth(100);
            ImGui::Text("Signal wavelength to plot");
            ImGui::InputDouble("##Signal wavelength to plot", &m_vec[curr_idx].temp_doubles[21], 0.0f, 0.0f, "%.2f");
            ImGui::SameLine();
            if (ImGui::Button("Apply to all profiles##5"))
            {
                for (auto i = 0; i < m_vec.size(); ++i)
                {
                    if (i == curr_idx) continue;
                    m_vec[i].temp_doubles[21] = m_vec[curr_idx].temp_doubles[21];
                }
            }
        
            ImGui::Text("First pump wavelength to plot");
            ImGui::InputDouble("##First pump wavelength to plot", &m_vec[curr_idx].temp_doubles[22], 0.0f, 0.0f, "%.2f");
            ImGui::SameLine();
            if (ImGui::Button("Apply to all profiles##6"))
            {
                for (auto i = 0; i < m_vec.size(); ++i)
                {
                    if (i == curr_idx) continue;
                    m_vec[i].temp_doubles[22] = m_vec[curr_idx].temp_doubles[22];
                }
            }
        
            ImGui::Text("Second pump wavelength to plot");
            ImGui::InputDouble("##Second pump wavelength to plot", &m_vec[curr_idx].temp_doubles[23], 0.0f, 0.0f, "%.2f");
            ImGui::SameLine();
            if (ImGui::Button("Apply to all profiles##7"))
            {
                for (auto i = 0; i < m_vec.size(); ++i)
                {
                    if (i == curr_idx) continue;
                    m_vec[i].temp_doubles[23] = m_vec[curr_idx].temp_doubles[23];
                }
            }
        
            ImGui::PopItemWidth();
            ImGui::Checkbox("Enable signal ASE", &m_vec[curr_idx].temp_bools[1]);
            ImGui::Checkbox("Save in dBm units", &m_vec[curr_idx].temp_bools[3]);
            
            ImGui::NewLine();
            ImGui::SameLine(ImGui::GetWindowWidth() * 0.4f);
            

            
            if (ImGui::Button("Run", button_size))
            {
                if (run_all)
                {
                    auto launch_monitor = [&it_idx, &monitor_launched, &buf, m_vec, n_threads]
                    (std::atomic<std::size_t>& available_threads, Parallel_idx& progress)
                    {
                        for (it_idx; it_idx < m_vec.size();)
                        {
                            if (available_threads > 1)
                            {
                                Simulation::Init_params p_{};
                                m_vec[it_idx].merge_with_sim(p_);
                                if (m_vec[it_idx].temp_bools[0]) p_.recalculate_constants();
                                std::string name(buf[it_idx]);
                                Launch_sim launcher {p_, m_vec[it_idx], name, n_threads};
                                std::thread t(launcher, std::ref(available_threads), std::ref(progress));
                                t.detach();
                                
                                if (it_idx == m_vec.size() - 1)
                                {
                                    it_idx = 0;
                                    monitor_launched = false;
                                    return;
                                }
                                
                                ++it_idx;
                            } else {
                                std::this_thread::sleep_for(std::chrono::milliseconds(500));
                            }
                        }
                        
                        return;
                    }; //end of launch_monitor
                    
                    if (!monitor_launched)
                    {
                        monitor_launched = true;
                        std::thread tm(launch_monitor, std::ref(available_threads), std::ref(progress));
                        tm.detach();
                    }
                } else {
                    if (available_threads > 0)
                    {
                        Simulation::Init_params p_{};
                        p_.set_steps(p_steps);
                        m_vec[curr_idx].merge_with_sim(p_);
                        if (m_vec[curr_idx].temp_bools[0]) p_.recalculate_constants();
                        std::string name(buf[curr_idx]);
                        Launch_sim launcher {p_, m_vec[curr_idx], name, n_threads};
                        std::thread t(launcher, std::ref(available_threads), std::ref(progress));
                        t.detach();
                    }
                }
            }
            

            ImGui::Checkbox("Run all profiles", &run_all);
            multi_progress_bar(progress, curr_idx);
            
        }//End of simulation
        
        ImGui::End(); // end window
        
        //ImGui::ShowDemoWindow();         //demo windows for reference
        window.clear(bgColor); // fill background with default color
        ImGui::SFML::Render(window);
        window.display();
    } //end of while (window.isOpen())
 
    ImGui::SFML::Shutdown();
    

    
    return 0;
}
