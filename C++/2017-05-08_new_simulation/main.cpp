#include "definitions.h"

using namespace std;

int main()
{
//1D grid test
/*
    Grid1DTM a {801, 420, 1550e-9, 1.0, 20, 4, 0}; //Max time, length, wavelength, source node, source type

    a.enable_lossy_termination();
    a.show_params();
    a.add_object(100, 160, 9.0, 1.0, 0, 0); //start, finish, eps_rel, mu_rel, sigma, sigma_mag
    a.show_params();
    a.set_dispersion(30);
    a.set_time_delay(10.0);
    a.save_information();
    for (unsigned long long i = 0; i < a.maximum_time(); ++i) {
        //if (i % 5 == 0) {a.save_state();}
        a.advance_simulation();
    }
*/

/*
    double courant = 1.0 / sqrt(2.0);
    unsigned long long duration = 301;
    Grid1DTM a {duration, 171, 1550e-9, courant, 20, 0, 0}; //Max time, length, wavelength, courant, ppw, source node, source type
    a.enable_lossy_termination();
    Params b {
        1550.0e-9, //Wavelength
        duration, //Max_time
        251, //Size_x
        251, //Size_y
        0, //Source type
        {75, 75}, //TFSF start
        {175, 175},  //TFSF end
        0           //TFSF choice, 0 completely closed, 1 open at the right hand side

    };

    Grid2DTMZ c {b};
    Point s {51, 75}; //Starting point of rectangle
    Point f {200, 150}; //Ending point of rectangle
    Point center {2, 125};
    //c.add_diel_srectangle(s, f, 4.0, 0, 1.0, 0); //Start, end, eps_rel, sigma, mu_rel, sigma_m
    //c.add_diel_circle(center, 45, 4.0, 0, 1.0, 0); //center, radius, eps_rel, sigma, mu_rel, sigma_m
    //c.add_pec_circle(center, 25); //center, radius
    //c.add_pec_rectangle({130,50}, {132, 150}); //start, finish

    for (unsigned long long time = 0; time < a.return_maximum_time(); ++time) {

        if (time % 10 == 0) {
      //          cout<<time<<'\n';
                c.save_state();
                a.save_state();
        }
        c.advance(a);
    }
*/
//2D grid test
/*
    double courant = 1.0 / sqrt(2.0);
    //double courant = 1.0;
    unsigned long long duration = 101; //Simulation duration
    unsigned long long source_loc = 0; //1D source loc
    double ppw = 28; //Set non default ppw
    double delay = (ppw - source_loc) / courant; //this is to make sure it starts at 0
    Grid1DTM a {duration, 101, 1550e-9, courant, ppw, source_loc, 0}; //Max time, length, wavelength, courant, ppw, source node, source type
    //a.set_time_delay(30);
    a.start_source_at_0(); //Makes functions start at 0
    Params b {
        1550.0e-9, //Wavelength
        duration, //Max_time
        ppw, //ppw
        201, //Size_x
        201, //Size_y
        0, //Source type
        {50, 50}, //TFSF start
        {150, 150},  //TFSF end
        0           //TFSF choice, 0 completely closed, 1 open at the right hand side

    };

    Grid2DTMZ c {b};
    Point s {75, 75}; //Starting point of rectangle
    Point f {125, 125}; //Ending point of rectangle
    Point center {100, 100};
    //c.add_diel_rectangle(s, f, 4.0, 0, 1.0, 0); //Start, end, eps_rel, sigma, mu_rel, sigma_m
    //c.add_diel_circle(center, 45, 4.0, 0, 1.0, 0); //center, radius, eps_rel, sigma, mu_rel, sigma_m
    //c.add_pec_circle(center, 50); //center, radius
    //c.add_pec_rectangle({130,50}, {132, 150}); //start, finish

    for (unsigned long long time = 0; time < a.return_maximum_time(); ++time) {

        if (time % 10 == 0) {
                cout<<time<<'\n';
                c.save_state();
                //a.save_state();
        }
        c.advance(a);
        //a.advance_simulation_hs();

    }
    cout<<"Plotting...\n";
    system("gnuplot plot_script_2d");
*/
//3D grid test
    double courant = 1.0 / sqrt(3.0);
    //double courant = 1.0;
    unsigned long long duration = 301; //Simulation duration
    double ppw = 22; //Set non default ppw

    Params b {
        1550.0e-9, //Wavelength
        duration, //Max_time
        ppw, //ppw
        51, //Size_x
        51, //Size_y
        51, //Size_z
        0, //Source type
        {12, 12, 12}, //TFSF start
        {24, 24, 24},  //TFSF end
        0, //TFSF radius
        0           //TFSF choice, 0 completely closed, 1 open at the right hand side

    };

    Grid3D c {b};
    //Point s {75, 75}; //Starting point of rectangle
    //Point f {125, 125}; //Ending point of rectangle
    //Point center {100, 100};
    //c.add_diel_rectangle(s, f, 4.0, 0, 1.0, 0); //Start, end, eps_rel, sigma, mu_rel, sigma_m
    //c.add_diel_circle(center, 45, 4.0, 0, 1.0, 0); //center, radius, eps_rel, sigma, mu_rel, sigma_m
    //c.add_pec_circle(center, 50); //center, radius
    //c.add_pec_rectangle({130,50}, {132, 150}); //start, finish


    for (unsigned long long time = 0; time < duration; ++time) {

        if (time % 10 == 0) {
                cout<<time<<'\n';
                c.save_state();
                c.zx_cross_section(25);
                c.xy_cross_section(25);
                c.zy_cross_section(25);
                //a.save_state();
        }
        c.update_test();
    }

    //cout<<"Plotting...\n";
    system("gnuplot plot_script_3d_e_xy");
    system("gnuplot plot_script_3d_e_zx");
    system("gnuplot plot_script_3d_e_zy");
    return 0;
}

