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


//3D grid simulation
    { //start of scope
    //Pending convert TFSF to parallel processing
    //convert_from_metres(double value, double wavelength, double eps_rel_max, double mu_rel_max, double ppw = 28);
    double courant = 1.0 / sqrt(3.0);
    unsigned long long duration = 1001; //Simulation duration
    double ppw = 10; //Set non default ppw
    unsigned long long source_loc = 0; //1D source loc

    Grid1DTM a {duration, 40, 1550e-9, courant, 16, source_loc, 0}; //Max time, length, wavelength, courant, ppw, source node, source type
    a.enable_lossy_termination();
    //Note: Grid1DTM must be 3 points larger than TFSF or 23 if enable_lossy_termination() is set
    a.start_source_at_0(); //Makes functions start at 0

    Params b {
        1550.0e-9, //Wavelength
        duration, //Max_time
        ppw, //ppw
        140, //Size_x
        592, //Size_y
        592, //Size_z
        0, //Source type
        {4, 4, 4}, //TFSF start
        {20, 584, 584},  //TFSF end
        0, //TFSF radius
        0           //TFSF choice, 0 completely closed, 1 open at the right hand side

    };

    Grid3D c {b};
    Point centre {16,296,296}; //optic fiber
    Point centre2 {108,296,296}; //waveguide and bulk glass covering top of the waveguide
    Point centre3 {108,296,161}; //bulk glass without covering the top of waveguide
    Point centre4 {22, 296, 296}; //cladding

    c.add_diel_rectangle_centre(centre2, 1, 1, 1, 2.56, 0, 1, 0); //increase the ppw to the max value
    c.add_diel_rectangle_centre(centre3, 48, 584, 312, 2.1025, 0, 1, 0); //add the bulk glass without covering top of the waveguide
    //c.add_diel_rectangle_centre(centre2, 48, 584, 584, 2.1025, 0, 1, 0); //add the bulk glass covering the top of the waveguide
    c.add_diel_rectangle_centre(centre2, 48, 42, 42, 2.56, 0, 1, 0); //add the waveguide

    c.add_diel_cylinder(centre4, 292, 42, 'x', 2.085136, 0, 1, 0); //cladding
    c.add_diel_cylinder(centre, 44, 48, 'x', 2.09525625, 0, 1, 0); //fiber core= centre, radius, length, axis alignment, eps_r, sigma, mu_r, sigma_m

    c.show_params();

    vector<thread> proc;
    for (unsigned long long time = 0; time < duration; ++time) {

        if (time % 10 == 0) {
                cout<<"time = "<<time<<'\n';
                if (time != duration - 1) {
                    proc.push_back(thread(&Grid3D::yz_cross_section, &c, 40, '1'));
                    proc.push_back(thread(&Grid3D::yz_cross_section, &c, 64, '2'));
                    proc.push_back(thread(&Grid3D::yz_cross_section, &c, 84, '3'));
                    proc.push_back(thread(&Grid3D::yz_cross_section, &c, 108, '4'));
                    for (thread &t : proc) {
                        t.join();
                    }
                    proc.clear();

                    proc.push_back(thread(&Grid3D::xy_cross_section, &c, 296, 'a'));
                    proc.push_back(thread(&Grid3D::xz_cross_section, &c, 296, 'a'));
                    for (thread &t : proc) {
                        t.join();
                    }
                    proc.clear();

                } else {

                    proc.push_back(thread(&Grid3D::yz_cross_section, &c, 40, '1'));
                    proc.push_back(thread(&Grid3D::yz_cross_section, &c, 64, '2'));
                    proc.push_back(thread(&Grid3D::yz_cross_section, &c, 84, '3'));
                    proc.push_back(thread(&Grid3D::yz_cross_section, &c, 108, '4'));
                    for (thread &t : proc) {
                        t.join();
                    }
                    proc.clear();

                    proc.push_back(thread(&Grid3D::xy_cross_section, &c, 296, 'a'));
                    proc.push_back(thread(&Grid3D::xz_cross_section, &c, 296, 'a'));
                    for (thread &t : proc) {
                        t.join();
                    }
                    proc.clear();
                    c.save_state();
                } //End of else

        } //End of if time%10
        c.parallel_advance_simulation(a);
    } //End of for loop


    } //End of scope to free memory


    cout<<"Creating plots...\n";

    vector<thread> plots;
    //plots.push_back(thread(&system, "gnuplot plot_script_3d_e"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e_yz-1"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e_yz-2"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e_yz-3"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e_yz-4"));
    for (thread &t : plots) {
        t.join();
    }
    plots.clear();

    plots.push_back(thread(&system, "gnuplot plot_script_3d_e_xy"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e_xz"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_h"));
    for (thread &t : plots) {
        t.join();
    }



/*
//3D grid testing area
    //convert_from_metres(double value, double wavelength, double eps_rel_max, double mu_rel_max, double ppw = 28);
    double courant = 1.0 / sqrt(3.0);
    unsigned long long duration = 1001; //Simulation duration
    double ppw = 10; //Set non default ppw
    unsigned long long source_loc = 0; //1D source loc

    Grid1DTM a {duration, 32, 1550e-9, courant, 10, source_loc, 0}; //Max time, length, wavelength, courant, ppw, source node, source type
    a.enable_lossy_termination();
    //Note: Grid1DTM must be 3 points larger than TFSF or 23 if enable_lossy_termination() is set
    a.start_source_at_0(); //Makes functions start at 0

    Params b {
        1550.0e-9, //Wavelength
        duration, //Max_time
        ppw, //ppw
        100, //Size_x
        30, //Size_y
        30, //Size_z
        0, //Source type
        {4, 4, 4}, //TFSF start
        {14, 26, 26},  //TFSF end
        0, //TFSF radius
        0           //TFSF choice, 0 completely closed, 1 open at the right hand side

    };

    Grid3D c {b};
    Point test = c.find_TFSF_centre();
    cout<<test<<'\n';
    cout<<"radius = "<<c.find_TFSF_radius();
    Point centre {10,15,15};
    Point centre2 {266,279,279};
    Point s {13, 0, 0}; //Starting point of rectangle
    Point f {14, 30, 30}; //Ending point of rectangle
    //c.add_pec_cylinder(centre, 7, 20, 'x'); //centre, radius, length, axis alignment
    c.add_pec_rectangle(s, f);

    c.add_diel_cylinder(centre, 9, 10, 'x', 1.1, 0, 1, 0); //centre, radius, length, axis alignment, eps_r, sigma, mu_r, sigma_m


    //c.add_diel_rectangle_centre(centre2, 86, 74, 74, 2.56, 0, 1.0, 0); //geometric centre, xlen, ylen, zlen, eps_rel,
                                                                        //sigma, mu_rel, sigma_m
    //c.add_diel_cylinder(centre, 98, 172, 'x', 2.0856, 0, 1, 0); //centre, radius, length, axis alignment, eps_r, sigma, mu_r, sigma_m
    //c.add_diel_cylinder(centre, 88, 172, 'x', 2.1025, 0, 1, 0); //centre, radius, length, axis alignment, eps_r, sigma, mu_r, sigma_m


    //c.add_pec_rectangle(s, f); //start, finish
    //c.add_pec_rectangle_centre(centre,  2, 26, 26); //centre, xlen, ylen, zlen
    //Point center {100, 100};
    //c.add_diel_rectangle(s, f, 4.0, 0, 1.0, 0); //Start, end, eps_rel, sigma, mu_rel, sigma_m
    //c.add_pec_sphere(centre, 7);
    //c.add_diel_sphere(centre, 7, 4, 0, 1.0, 0); //centre, radius, eps_rel, sigma, mu_rel, sigma_m
    //c.add_diel_rectangle_centre(centre2,  7, 14, 14, 1.0, 0, 1.0, 0);
    c.show_params();


    for (unsigned long long time = 0; time < duration; ++time) {
        if (time % 10 == 0) {
        //if (129 < time && time < 171) {
                cout<<"time = "<<time<<'\n';
                c.parallel_save_all(13, 15, 15, 'a');
                c.yz_cross_section(14, 'b');
                c.yz_cross_section(16, 'c');
        }

        c.parallel_advance_simulation(a);
    }


    cout<<"Creating plots...\n";

    vector<thread> plots;
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e_xy"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e_xz"));
    plots.push_back(thread(&system, "gnuplot plot_script_3d_e_yz"));
    for (thread &t : plots) {
        t.join();
    }
*/
    return 0;
}

