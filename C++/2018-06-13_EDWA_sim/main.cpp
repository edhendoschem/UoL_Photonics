#include "headers.h"
#include "utility_functions.h"
//Pending recent
//Pending long term
//Add overlap columns to the GUI
//Add function to calculate overlap scattering
//Investigate power in dbm for signal

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
        
        for (auto i = 11; i < temp_doubles.size(); ++i)
        {
            temp_doubles[i] = 0.0;
        }
        
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
        p.l = temp_doubles[2] * 1.0e-2;
        p.NEr = temp_doubles[3];
        p.NYb = temp_doubles[4];
        p.A21 = 1.0/(temp_doubles[5] * 1.0e-3);
        p.A32 = 1.0/(temp_doubles[6] * 1.0e-3);
        p.A43 = 1.0/(temp_doubles[7] * 1.0e-3);
        p.A65 = 1.0/(temp_doubles[8] * 1.0e-3);
        p.step_size = p.l / static_cast<double>(p.steps-1);
        p.Ps0 = temp_signals;
        p.Pp0_f = temp_f_pump;
        p.Pp0_b = temp_b_pump;
        p.Cup = temp_doubles[9];
        p.Ccr = temp_doubles[10];
        
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
        progress.idx[curr_t].first = 0.0f;
        r.simulate(progress.idx[curr_t].first, false, m.temp_bools[1], m.temp_bools[2]);
        r.save_data(name, m.temp_bools[3]);
        r.plot_data(name);
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
        int temp_wl {static_cast<int>(m_vec[curr_idx].temp_doubles[15] * 1000.0)};
        if (!(m_vec[curr_idx].temp_doubles[15] < 850.0 || m_vec[curr_idx].temp_doubles[15] > 1150.0))
        {
            double temp_pow {m_vec[curr_idx].temp_doubles[16] / 1000.0};
            m_vec[vec_idx].temp_b_pump.emplace(std::pair(temp_wl, temp_pow));
        }
    } else if (t == S_type::Signal)
    {
        int temp_wl {static_cast<int>(m_vec[curr_idx].temp_doubles[11] * 1000.0)};
        if (!(m_vec[curr_idx].temp_doubles[11] < 1450.0 || m_vec[curr_idx].temp_doubles[11] > 1600.0))
        {
            double temp_pow {m_vec[curr_idx].temp_doubles[12] / 1000.0};
            m_vec[vec_idx].temp_signals.emplace(std::pair(temp_wl, temp_pow));
        }
    } else 
    {
        int temp_wl {static_cast<int>(m_vec[curr_idx].temp_doubles[13] * 1000.0)};
        if (!(m_vec[curr_idx].temp_doubles[13] < 850.0 || m_vec[curr_idx].temp_doubles[13] > 1150.0))
        {
            double temp_pow {m_vec[curr_idx].temp_doubles[14] / 1000.0};
            m_vec[vec_idx].temp_f_pump.emplace(std::pair(temp_wl, temp_pow));
        }
    }
    
    return;
}

int main(int argc, char **argv)
{
    std::size_t n_threads {std::thread::hardware_concurrency()}; //Number of threads
    std::atomic<std::size_t> available_threads {n_threads}; //Threads not in use
    //Support 12 different profiles
    char const * buf[] {"profile_0", "profile_1", "profile_2", "profile_3", "profile_4",
                              "profile_5", "profile_6", "profile_7", "profile_8", "profile_9",
                              "profile_10", "profile_11"};
    Simulation::Init_params p{};
    bool run_all {true};
    bool copy_all {true};
    bool monitor_launched {false};
    bool ready {false};
    
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
    char windowTitle[255] = "Test window";
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
                ImGui::InputInt("Steps", &p.steps);
                ImGui::InputDouble("Cup (m^3/s)", &m_vec[curr_idx].temp_doubles[9], 0.0f, 0.0f, "%e");
                ImGui::InputDouble("Ccr (m^3/s)", &m_vec[curr_idx].temp_doubles[10], 0.0f, 0.0f, "%e");
                ImGui::PopItemWidth();
                
                if (ImGui::Button("Apply to all profiles##3"))
                {
                    for (auto i = 0; i < m_vec.size(); ++i)
                    {
                        if (i != curr_idx)
                        {
                            m_vec[i].temp_doubles[9] = m_vec[curr_idx].temp_doubles[9];
                            m_vec[i].temp_doubles[10] = m_vec[curr_idx].temp_doubles[10];
                        }
                    }
                }

                ImGui::TreePop();
                ImGui::Separator();
            }
            
            if (ImGui::TreeNode("Signals & Pump"))
            {
                ImGui::Text("Wavelength (nm)"); ImGui::SameLine(150); ImGui::Text("Power (mW)");
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w1", &m_vec[curr_idx].temp_doubles[11], 0.0f, 0.0f, "%.2f");
                ImGui::PopItemWidth();
                ImGui::SameLine(150);
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w2", &m_vec[curr_idx].temp_doubles[12], 0.0f, 0.0f, "%e");
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
                ImGui::InputDouble("##w3", &m_vec[curr_idx].temp_doubles[13], 0.0f, 0.0f, "%.2f");
                ImGui::PopItemWidth();
                ImGui::SameLine(150);
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w4", &m_vec[curr_idx].temp_doubles[14], 0.0f, 0.0f, "%e");
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
                ImGui::InputDouble("##w5", &m_vec[curr_idx].temp_doubles[15], 0.0f, 0.0f, "%.2f");
                ImGui::PopItemWidth();
                ImGui::SameLine(150);
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w6", &m_vec[curr_idx].temp_doubles[16], 0.0f, 0.0f, "%e");
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
            ImGui::Checkbox("Enable signal ASE", &m_vec[curr_idx].temp_bools[1]);
            ImGui::Checkbox("Enable pump ASE", &m_vec[curr_idx].temp_bools[2]);
            ImGui::Checkbox("Save in dBm units", &m_vec[curr_idx].temp_bools[3]);
            ImGui::Checkbox("Run all profiles", &run_all);
            ImGui::NewLine();
            ImGui::SameLine(ImGui::GetWindowWidth() * 0.4f);
            

            
            if (ImGui::Button("Run", button_size))
            {
                if (run_all)
                {//aquiaqui
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
                                    break;
                                }
                                
                                ++it_idx;
                            } else {
                                std::this_thread::sleep_for(std::chrono::milliseconds(500));
                            }
                        }
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
                        m_vec[curr_idx].merge_with_sim(p_);
                        if (m_vec[curr_idx].temp_bools[0]) p_.recalculate_constants();
                        std::string name(buf[curr_idx]);
                        Launch_sim launcher {p_, m_vec[curr_idx], name, n_threads};
                        std::thread t(launcher, std::ref(available_threads), std::ref(progress));
                        t.detach();
                    }
                }
            }
            
            multi_progress_bar(progress, curr_idx);
            
        }//End of simulation
        
        ImGui::End(); // end window
        
        //ImGui::ShowDemoWindow();         //demo windows for reference
        window.clear(bgColor); // fill background with default color
        ImGui::SFML::Render(window);
        window.display();
    }
 
    ImGui::SFML::Shutdown();


    
    return 0;
}
