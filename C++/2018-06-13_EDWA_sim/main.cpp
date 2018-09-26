#include "headers.h"
#include "utility_functions.h"
//Pending recent
//Pending long term
//Add function to calculate overlap scattering
//Investigate power in dbm for signal

struct Merge_helper
{
    std::array<double, 30> temp_doubles;
    std::array<bool, 30> temp_bools;
    wl_map temp_signals;
    wl_map temp_f_pump;
    wl_map temp_b_pump;
    
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
        
        for (auto i = 9; i < temp_doubles.size(); ++i)
        {
            temp_doubles[i] = 0.0;
        }
        
        temp_bools[0] = true;
        
        return;
    }
    
    
    void merge_with_sim(Simulation::Init_params& p)
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
    
        return;
    }
}; //End of Merge_helper struct


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


int main(int argc, char **argv)
{
    Simulation::Init_params p{};
    
    Merge_helper m;
    m.init_temp(p);
    
    sf::RenderWindow window(sf::VideoMode(800, 800), "");
    window.setVerticalSyncEnabled(true);
    ImGui::SFML::Init(window);
 
    sf::Color bgColor;
    char windowTitle[255] = "Test window";
    ImVec2 button_size(100, 50);

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
            if (ImGui::TreeNode("Waveguide Dimensions"))
            {
                ImGui::SameLine(); ShowHelpMarker("test help");
                ImGui::PushItemWidth(50);
                ImGui::InputDouble("Width (um)", &m.temp_doubles[0], 0.0f, 0.0f, "%.2f");
                ImGui::InputDouble("Height (um)", &m.temp_doubles[1], 0.0f, 0.0f, "%.2f");
                ImGui::InputDouble("Length (cm)", &m.temp_doubles[2], 0.0f, 0.0f, "%.2f"); //3 decimals
                ImGui::PopItemWidth();
                ImGui::TreePop();
                ImGui::Separator();
            }
            
            if (ImGui::TreeNode("Rare earth content"))
            {
                ImGui::SameLine(); ShowHelpMarker("test help");
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("Erbium (ions/m^3)", &m.temp_doubles[3], 0.0f, 0.0f, "%e");
                ImGui::InputDouble("Ytterbium (ions/m^3)", &m.temp_doubles[4], 0.0f, 0.0f, "%e");
                ImGui::PopItemWidth();
                ImGui::TreePop();
                ImGui::Separator();
            }
            
            if (ImGui::TreeNode("Lifetime information"))
            {
                ImGui::SameLine(); ShowHelpMarker("test help");
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("t21 (ms)", &m.temp_doubles[5], 0.0f, 0.0f, "%e");
                ImGui::InputDouble("t32 (ms)", &m.temp_doubles[6], 0.0f, 0.0f, "%e");
                ImGui::InputDouble("t43 (ms)", &m.temp_doubles[7], 0.0f, 0.0f, "%e");
                ImGui::InputDouble("t65 (ms)", &m.temp_doubles[8], 0.0f, 0.0f, "%e");
                ImGui::PopItemWidth();
                ImGui::TreePop();
                ImGui::Separator();
            }
            
            if (ImGui::TreeNode("Other information"))
            {
                ImGui::SameLine(); ShowHelpMarker("test help");
                ImGui::PushItemWidth(100);
                ImGui::InputInt("Steps", &p.steps);
                ImGui::InputDouble("Cup (m^3/s)", &p.Cup, 0.0f, 0.0f, "%e");
                ImGui::InputDouble("Ccr (m^3/s)", &p.Ccr, 0.0f, 0.0f, "%e");
                ImGui::PopItemWidth();
                ImGui::TreePop();
                ImGui::Separator();
            }
            
            if (ImGui::TreeNode("Signals & Pump"))
            {
                ImGui::Text("Wavelength (nm)"); ImGui::SameLine(150); ImGui::Text("Power (mW)");
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w1", &m.temp_doubles[9], 0.0f, 0.0f, "%.2f");
                ImGui::PopItemWidth();
                ImGui::SameLine(150);
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w2", &m.temp_doubles[10], 0.0f, 0.0f, "%e");
                ImGui::PopItemWidth();
                ImGui::SameLine();
                if (ImGui::Button("Add signal"))
                {
                    int temp_wl {static_cast<int>(m.temp_doubles[9] * 1000.0)};
                    if (!(m.temp_doubles[9] < 1450.0 || m.temp_doubles[9] > 1600.0))
                    {
                        double temp_pow {m.temp_doubles[10] / 1000.0};
                        m.temp_signals.emplace(std::pair(temp_wl, temp_pow));
                    }
                }
                ImGui::SameLine();
                if (ImGui::Button("Clear all##1"))
                {
                    m.temp_signals.clear();
                }
                
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w3", &m.temp_doubles[11], 0.0f, 0.0f, "%.2f");
                ImGui::PopItemWidth();
                ImGui::SameLine(150);
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w4", &m.temp_doubles[12], 0.0f, 0.0f, "%e");
                ImGui::PopItemWidth();
                ImGui::SameLine();
                
                if (ImGui::Button("Add Forward pump"))
                {
                    int temp_wl {static_cast<int>(m.temp_doubles[11] * 1000.0)};
                    if (!(m.temp_doubles[11] < 850.0 || m.temp_doubles[11] > 1150.0))
                    {
                        double temp_pow {m.temp_doubles[12] / 1000.0};
                        m.temp_f_pump.emplace(std::pair(temp_wl, temp_pow));
                    }
                }
                
                ImGui::SameLine();
                
                if (ImGui::Button("Clear all##2"))
                {
                    m.temp_f_pump.clear();
                }
                
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w5", &m.temp_doubles[13], 0.0f, 0.0f, "%.2f");
                ImGui::PopItemWidth();
                ImGui::SameLine(150);
                ImGui::PushItemWidth(100);
                ImGui::InputDouble("##w6", &m.temp_doubles[14], 0.0f, 0.0f, "%e");
                ImGui::PopItemWidth();
                ImGui::SameLine();
                
                if (ImGui::Button("Add backward pump"))
                {
                    int temp_wl {static_cast<int>(m.temp_doubles[13] * 1000.0)};
                    if (!(m.temp_doubles[13] < 850.0 || m.temp_doubles[13] > 1150.0))
                    {
                        double temp_pow {m.temp_doubles[14] / 1000.0};
                        m.temp_b_pump.emplace(std::pair(temp_wl, temp_pow));
                    }
                }
                
                ImGui::SameLine();
                
                if (ImGui::Button("Clear all##3"))
                {
                    m.temp_b_pump.clear();
                }
                
                ImGui::Separator();
                ImGui::Columns(3, "test columns");
                ImGui::Text("Wavelength (nm)"); ImGui::NextColumn();
                ImGui::Text("Power (mW)"); ImGui::NextColumn();
                ImGui::Text("Type"); ImGui::NextColumn();
                ImGui::Separator();
                
                
                if (m.temp_signals.size() != 0)
                {
                    log_map(m.temp_signals, "Signal");
                }
                
                if (m.temp_f_pump.size() != 0)
                {
                    log_map(m.temp_f_pump, "Forward pump");
                }
                
                if (m.temp_b_pump.size() != 0)
                {
                    log_map(m.temp_b_pump, "Backward pump");
                }
                
                ImGui::Columns(1);
                
                ImGui::TreePop();
                ImGui::Separator();
            }
            
            ImGui::Checkbox("Recalculate Cup and Ccr", &m.temp_bools[0]);
            ImGui::SameLine(); ShowHelpMarker("Will use the default equations to calculate the values. I.E. will overwrite user provided values\n");
            ImGui::NewLine();
            ImGui::SameLine(ImGui::GetWindowWidth() * 0.4f);
            if (ImGui::Button("Update", button_size))
            {
                m.merge_with_sim(p);
                if (m.temp_bools[0]) p.recalculate_constants();
            }
            
        } //End of starting parameters
        
        
        if (ImGui::CollapsingHeader("Simulation"))
        {
            ImGui::Checkbox("Enable signal ASE", &m.temp_bools[1]);
            ImGui::Checkbox("Enable pump ASE", &m.temp_bools[2]);
            ImGui::Checkbox("Save in dBm units", &m.temp_bools[3]);
            static char buf1[64] = ""; ImGui::InputText("Data name", buf1, 64);
            static char buf2[64] = ""; ImGui::InputText("Plot script name", buf2, 64);

            ImGui::NewLine();
            ImGui::SameLine(ImGui::GetWindowWidth() * 0.4f);

            if (ImGui::Button("Run", button_size))
            {
                Simulation::Result r{p};
                r.simulate(false, m.temp_bools[1], m.temp_bools[2]);
                ImGui::ProgressBar(r.p.curr_step, ImVec2(0.0f,0.0f));
                r.save_data(std::string(buf1), m.temp_bools[3]);
                r.plot_data(std::string(buf1), std::string(buf2));
            }
        
        }//End of simulation
        
        ImGui::End(); // end window
        
        ImGui::ShowDemoWindow();         // Background color edit
        window.clear(bgColor); // fill background with default color
        ImGui::SFML::Render(window);
        window.display();
    }
 
    ImGui::SFML::Shutdown();


    
    return 0;
}
