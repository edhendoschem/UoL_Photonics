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
    Grid1DTM a {duration, 123, 1550e-9, courant, ppw, source_loc, 0}; //Max time, length, wavelength, courant, ppw, source node, source type
    //a.set_time_delay(30);
    a.start_source_at_0(); //Makes functions start at 0
    a.enable_lossy_termination();
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
    //Pending convert TFSF to parallel processing
    //convert_from_metres(double value, double wavelength, double eps_rel_max, double mu_rel_max, double ppw = 28);
    double courant = 1.0 / sqrt(3.0);
    unsigned long long duration = 301; //Simulation duration
    double ppw = 20; //Set non default ppw
    unsigned long long source_loc = 0; //1D source loc

    Grid1DTM a {duration, 53, 1550e-9, courant, 20, source_loc, 0}; //Max time, length, wavelength, courant, ppw, source node, source type
    a.enable_lossy_termination();
    //Note: Grid1DTM must be 3 points larger than TFSF or 23 if enable_lossy_termination() is set
    a.start_source_at_0(); //Makes functions start at 0

    Params b {
        1550.0e-9, //Wavelength
        duration, //Max_time
        ppw, //ppw
        50, //Size_x
        50, //Size_y
        50, //Size_z
        0, //Source type
        {10, 10, 10}, //TFSF start
        {40, 40, 40},  //TFSF end
        0, //TFSF radius
        0           //TFSF choice, 0 completely closed, 1 open at the right hand side

    };

    Grid3D c {b};
    Point centre {15,25,25};
    Point centre2 {22,25,25};
    Point s {40, 12, 12}; //Starting point of rectangle
    Point f {50, 37, 37}; //Ending point of rectangle
    //c.add_pec_cylinder(centre, 7, 20, 'x'); //centre, radius, length, axis alignment
    //c.add_diel_cylinder(centre, 7, 35, 'x', 9.0, 0, 1, 0); //centre, radius, length, axis alignment, eps_r, sigma, mu_r, sigma_m
    //c.add_pec_rectangle(s, f); //start, finish
    //c.add_pec_rectangle_centre(centre,  2, 26, 26); //centre, xlen, ylen, zlen
    //Point center {100, 100};
    //c.add_diel_rectangle(s, f, 4.0, 0, 1.0, 0); //Start, end, eps_rel, sigma, mu_rel, sigma_m
    //c.add_diel_rectangle_centre(centre, 10, 10, 10, 4.0, 0, 1.0, 0); //centre, xlen, ylen, zlen, eps_rel, sigma, mu_rel,
                                                                        //sigma_m
    //c.add_pec_sphere(centre, 7);
    //c.add_diel_sphere(centre, 7, 4, 0, 1.0, 0); //centre, radius, eps_rel, sigma, mu_rel, sigma_m
    //c.add_diel_rectangle_centre(centre2,  7, 14, 14, 1.0, 0, 1.0, 0);
    c.show_params();
/*
    for (unsigned long long time = 0; time < duration; ++time) {
        if (time % 10 == 0) {
        //if (129 < time && time < 171) {
                cout<<"time = "<<time<<'\n';

                c.save_state();
                c.xz_cross_section(25, 'a');
                c.xy_cross_section(25, 'a');
                c.yz_cross_section(25, 'a');

                c.parallel_save_all(25, 'a');
                //a.save_state();
        }
        c.advance_simulation(a);
    }
*/

    cout<<"Creating plots...\n";

    vector<thread> plots;
    //plots.push_back(thread(&system, "gnuplot plot_script_3d_e"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e_xy"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e_xz"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e_yz"));
    for (thread &t : plots) {
        t.join();
    }

    //convert_from_metres(distance in metres, wavelength, eps_r, mu_r, ppw)
    cout<<"conversion to nodes= "<<convert_from_metres(1.55e-6, 1.55e-6, 2.56, 1.0, 20)<<'\n';

    //convert_from_nodes(nodes, wavelength, eps_r, mu_r, ppw, courant)
    //cout<<"conversion to metres = "<<convert_from_nodes(64.0, 1.55e-6, 2.56, 1.0, 18)<<'\n';


    return 0;
}

