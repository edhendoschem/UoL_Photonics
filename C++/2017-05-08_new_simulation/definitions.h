/*--------------File containing all definitions, #1--------------*/
#ifndef DEFINITIONS_H_INCLUDED
#define DEFINITIONS_H_INCLUDED

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream> //File streams
#include <iomanip>
#include <complex>

//Multithreading
#include <thread> //Note for clang the linker must have -lpthread flag

//Constants
const double PI = 3.14159265358979; //Pi constant
const double eps_0 = 8.854187e-12; //Permittivity of vacuum/free space
const double mu_0 = 4.0 * PI * 1.0e-7; //Permeability of free space
const double c0 = 299792458.0; //Speed of light in vacuum m/s
const double imp0 = mu_0 * c0; //Characteristic impedance of free space Ohms
const unsigned int NUM_THREADS = std::thread::hardware_concurrency(); //Detects maximum number of threads

using namespace std;


class Flat_vec {
public:
    //Non-default constructor 1D vector
    Flat_vec(unsigned long long lengthh);
    //Non-default constructor 2D vector
    Flat_vec(unsigned long long lengthh, unsigned long long widthh);
    //Non-default constructor 3D vector
    Flat_vec(unsigned long long lengthh, unsigned long long widthh,unsigned long long depthh);
    //get index function
    unsigned long long get_ind_3D(unsigned long long i, unsigned long long j, unsigned long long k);
    unsigned long long get_ind_2D(unsigned long long i, unsigned long long j);
    //get size function
    unsigned long long size() const {return vec_size;}
    //Display vector's dimensions
    void display_dimensions() {cout<<"length = "<<length<<'\t'<<"width = "<<width<<'\t'<<"depth = "<<depth<<'\n';}

    //Operator overloading for access
    double &operator () (unsigned long long i, unsigned long long j, unsigned long long k) {return vec[get_ind_3D(i, j, k)];}
    double &operator () (unsigned long long i, unsigned long long j) {return vec[get_ind_2D(i, j)];}
    double &operator () (unsigned long long i) {return vec[i];}
private:
    vector<double> vec;
    unsigned long long length, width, depth, vec_size;
};


class Grid1DTM {
public:

    //Non-default constructor with source parameters
    Grid1DTM(unsigned long long max_timee, unsigned long long lengthh, double wavelengthh, double courantt,
             double ppww, unsigned long long source_nodee, int typee);
    //Display current params
    void show_params() const;
    void change_courant(double new_val);
    void change_dx(double new_val);
    void change_dt(double new_val);
    void change_ppw(double new_val);
    void place_source(unsigned long long node);
    void advance_simulation(); //advances the simulation 1 time unit
    void advance_simulation_hs(); //advances simulation 1 time unit using hard source
    void save_state(); //Save the vectors for plotting
    void enable_lossy_termination(); //Lossy termination in the last 20 nodes for insertion in 2D TFSF
    unsigned long long return_maximum_time() const; //Returns max_time
    unsigned long long return_current_time() const; //Returns current time
    void add_object(unsigned long long start, unsigned long long finish, double eps_rel, double mu_rel, double sigma,
                        double sigma_mag);
    void set_time_delay(double val);
    void set_dispersion(double val);
    void start_source_at_0();
    void save_information();
    double return_Ez(unsigned long long i);
    double return_Hy(unsigned long long i);


private:
    Flat_vec Ez, c1ez, c2ez;
    Flat_vec Hy, c1hy, c2hy;
    Flat_vec Ez_left, Ez_right;
    unsigned long long length; //Grid size
    unsigned long long curr_time, max_time; //Time parameters
    unsigned long long source_node; //where is the source going to be placed
    double A, B,C, D, time_delay = 30, dispersion = 10; //can be any number but ideally integers since they are in delta t units
    int type; //source type
    //Courant number, points per wavelength, wavelength (meters), spatial step, time step
    double courant, ppw, wavelength, dx, dt;
    //Stores the highest eps_rel and mu_rel and minimum wavelength
    double min_wavelength;
    double source(); //calculates the value of the source for the current time
    void update_magnetic(); //Updates the magnetic field vector
    void update_electric(); //Updates the electric field vector
    void update_magnetic_hs(); //Update magnetic field with hard source
    void update_electric_hs(); //update electric field with a hard source
    void abc_left(); //Apply second order absorbing boundary to the left
    void abc_right(); //Apply second order absorbing boundary to the right
    void update_abc(); //First order absorbing boundary
};


class Grid1DTE {
public:

    //Non-default constructor with source parameters
    Grid1DTE(unsigned long long max_timee, unsigned long long lengthh, double wavelengthh, double courantt,
             double ppww, unsigned long long source_nodee, int typee);
    //Display current params
    void show_params() const;
    void change_courant(double new_val);
    void change_dx(double new_val);
    void change_dt(double new_val);
    void change_ppw(double new_val);
    void place_source(unsigned long long node);
    void select_source(int type = 0); //0 = Harmonic, 1 = ?
    void advance_simulation(); //advances the simulation 1 time unit
    void advance_simulation_hs(); //ADvances simulation 1 time unit using hard source
    void save_state(); //Save the vectors for plotting
    void enable_lossy_termination(); //Lossy termination in the last 20 nodes for insertion in 2D TFSF
    unsigned long long return_maximum_time() const; //Returns max_time
    unsigned long long return_current_time() const; //Returns current time
    void add_object(unsigned long long start, unsigned long long finish, double eps_rel, double mu_rel, double sigma,
                        double sigma_mag);
    void set_time_delay(double val);
    void start_source_at_0();
    void set_dispersion(double val);
    void save_information();
    double return_Hz(unsigned long long i);
    double return_Ey(unsigned long long i);

private:
    Flat_vec Hz, c1hz, c2hz;
    Flat_vec Ey, c1ey, c2ey;
    Flat_vec Ey_left, Ey_right;
    unsigned long long length; //Grid size
    unsigned long long curr_time, max_time; //Time parameters
    unsigned long long source_node; //where is the source going to be placed
    double A, B, C, D, time_delay = 30, dispersion = 10; //can be any number but ideally integers since they are in delta t units
    int type; //source type
    //Courant number, points per wavelength, wavelength (meters), spatial step, time step
    double courant, ppw, wavelength, dx, dt;
    //Stores the highest eps_rel and mu_rel and minimum wavelength
    double min_wavelength;
    double source(); //calculates the value of the source for the current time
    void update_magnetic(); //Updates the magnetic field vector
    void update_electric(); //Updates the electric field vector
    void update_magnetic_hs(); //Updates the magnetic field vector with hard source
    void update_electric_hs(); //Updates the electric field vector with hard source7
    void abc_left(); //Apply second order absorbing boundary to the left
    void abc_right(); //Apply second order absorbing boundary to the right
    void update_abc(); //First order absorbing boundary
};


//Helper struct for setting points in the grid, only takes positive values as the grid has no negative values
struct Point {
    unsigned long long x, y, z;
    //2D constructor
    Point(unsigned long long xx, unsigned long long yy) : x {xx}, y {yy} {};
    //3D constructor
    Point(unsigned long long xx, unsigned long long yy, unsigned long long zz) : x {xx}, y {yy}, z{zz} {};

};

//Output operator for Point
ostream &operator << (ostream &os, Point &a);

struct Params {
    //2D constructor
    Params(double wavelengthh, unsigned long long max_timee, double ppww, unsigned long long size_xx,
           unsigned long long size_yy, int source_typee, Point tfsf_startt, Point tfsf_endd, int choicee);
    //3D constructor
    Params(double wavelengthh, unsigned long long max_timee, double ppww, unsigned long long size_xx,
           unsigned long long size_yy, unsigned long long size_zz, int source_typee, Point tfsf_startt,
           Point tfsf_endd, unsigned long long tfsf_radiuss, int choicee);

    unsigned long long size_x, size_y, size_z, max_time, tfsf_radius;
    Point tfsf_start, tfsf_end;
    int source_type, choice;
    double wavelength, ppw;
};

//Mode TMz
class Grid2DTMZ {
public:
    Grid2DTMZ(Params in);
    void show_params();
    void advance(Grid1DTM &aux_grid);
    unsigned long long return_maximum_time();
    unsigned long long return_current_time();
    Point size();
    void save_state();
    void add_diel_rectangle(Point start, Point finish, double eps_r, double sigma = 0, double mu_r = 1, double sigma_m = 0);
    void add_diel_circle(Point center, unsigned long long radius, double eps_r, double sigma = 0, double mu_r = 1,
                         double sigma_m = 0);
    void add_pec_circle(Point center, unsigned long long radius);
    void add_pec_rectangle(Point start, Point finish);

private:
    Flat_vec Hx, Hy, Ez, c1hx, c2hx, c1hy, c2hy, c1ez, c2ez, Ez_top, Ez_bottom, Ez_left, Ez_right;
    unsigned long long curr_time = 0, max_time, size_x, size_y;
    Point tfsf_start, tfsf_end;
    int source_type, choice = 0;
    double time_delay, dispersion = 10, ppw, wavelength, min_wavelength, dx, dt, courant;
    double cour_prime, A, B, C, D; //Second order absorbing boundary condition coefficients

    bool bounds_check(Point a); //returns true if point is inside boundary. Note it depends on unsigned wraparound behaviour
    void update_magnetic();
    void update_electric();
    //Absorbing boundary conditions helper functions,
    //Note: Lossy dielectric terminations may not work properly as that would require more vectors to keep track of the
    //electric and magnetic conductivity
    void abc_top();
    void abc_bottom();
    void abc_left();
    void abc_right();
    void apply_abc(); //applies second order absorbing boundary condition
    //TFSF function
    void apply_TFSF(Grid1DTM &aux_grid);
    void apply_open_right_TFSF(Grid1DTM &aux_grid);


};

//Helper function to convert metres to number of nodes
//default ppw set for 2d simulation, note min ppw is 10, 14, 18 for 1D, 2D and 3D respectively
unsigned long long convert_from_metres(double value, double wavelength, double eps_rel_max, double mu_rel_max,
                                       double ppw = 28);


//Mode TMz
class Grid2DTEZ {
public:
    Grid2DTEZ(Params in);
    void show_params();
    void advance(Grid1DTE &aux_grid);
    unsigned long long return_maximum_time();
    unsigned long long return_current_time();
    Point size();
    void save_state();
    void add_diel_rectangle(Point start, Point finish, double eps_r, double sigma = 0, double mu_r = 1, double sigma_m = 0);
    void add_diel_circle(Point center, unsigned long long radius, double eps_r, double sigma = 0, double mu_r = 1,
                         double sigma_m = 0);
    void add_pec_circle(Point center, unsigned long long radius);
    void add_pec_rectangle(Point start, Point finish);

private:
    Flat_vec Hz, Ey, Ex, c1hz, c2hz, c1ey, c2ey, c1ex, c2ex, Ex_top, Ex_bottom, Ey_left, Ey_right;
    unsigned long long curr_time = 0, max_time, size_x, size_y;
    Point tfsf_start, tfsf_end;
    int source_type, choice = 0;
    double time_delay = 30, dispersion = 10, ppw, wavelength, min_wavelength, dx, dt, courant = 1.0 / sqrt(2.0);
    double cour_prime, A, B, C, D; //Second order absorbing boundary condition coefficients

    bool bounds_check(Point a); //returns true if point is inside boundary. Note it depends on unsigned wraparound behaviour
    void update_magnetic();
    void update_electric();
    //Absorbing boundary conditions helper functions,
    //Note: Lossy dielectric terminations may not work properly as that would require more vectors to keep track of the
    //electric and magnetic conductivity
    void abc_top();
    void abc_bottom();
    void abc_left();
    void abc_right();
    void apply_abc(); //applies second order absorbing boundary condition
    //TFSF function
    void apply_TFSF(Grid1DTE &aux_grid);
    void apply_open_right_TFSF(Grid1DTE &aux_grid);

};

class Grid3D {
public:
    //Non-default constructor
    Grid3D(Params in);

    void save_state();
    void update_test(Grid1DTM &aux_grid);
    void xz_cross_section(unsigned long long y, char w);
    void xy_cross_section(unsigned long long z, char w);
    void yz_cross_section(unsigned long long x, char w);
    void test();
    void add_pec_rectangle(Point start, Point finish);
    void add_pec_rectangle_centre(Point centre, unsigned long long xlen, unsigned long long ylen, unsigned long long zlen);
    void add_diel_rectangle(Point start, Point finish, double eps_r, double sigma, double mu_r, double sigma_m);
    void add_diel_rectangle_centre(Point centre, unsigned long long xlen, unsigned long long ylen, unsigned long long zlen,
            double eps_r, double sigma, double mu_r, double sigma_m);
    void add_diel_sphere(Point centre, unsigned long long radius, double eps_r, double sigma, double mu_r, double sigma_m);
    void add_pec_sphere(Point centre, unsigned long long radius);
    void add_pec_cylinder(Point centre, unsigned long long radius, unsigned long long length, char opt);
    void add_diel_cylinder(Point centre, unsigned long long radius, unsigned long long length, char opt,
            double eps_r, double sigma, double mu_r, double sigma_m);

private:
    Flat_vec Ex, Ey, Ez, Hx, Hy, Hz; //Electric and magnetic fields
    Flat_vec c1ex, c2ex, c1ey, c2ey, c1ez, c2ez, c1hx, c2hx, c1hy, c2hy, c1hz, c2hz; //Field coefficients
    Flat_vec Ey_x0, Ez_x0, Ey_xf, Ez_xf, Ex_y0, Ez_y0, Ex_yf, Ez_yf, Ex_z0, Ey_z0, Ex_zf, Ey_zf; //ABC coeff
    Point tfsf_start, tfsf_end;
    unsigned long long curr_time = 0, max_time, size_x, size_y, size_z;
    int source_type, choice = 0;
    double time_delay = 30, dispersion = 10, ppw, wavelength, min_wavelength, dx, dt, courant = 1.0 / sqrt(3.0);
    double cour_prime, A, B, C, D; //Second order absorbing boundary condition coefficients

    void update_Hx();
    void update_Hy();
    void update_Hz();
    void update_Ex();
    void update_Ey();
    void update_Ez();
    void update_magnetic();
    void update_electric();

    //First order ABC (less costly than second order abc)
    void simple_x0();
    void simple_xf();
    void simple_y0();
    void simple_yf();
    void simple_z0();
    void simple_zf();
    //Apply all of the previous functions
    void simple_abc();

    //TFSF formulation
    void correct_Hy_0(Grid1DTM &aux_grid);
    void correct_Hy_f(Grid1DTM &aux_grid);
    void correct_Hx_0(Grid1DTM &aux_grid);
    void correct_Hx_f(Grid1DTM &aux_grid);
    void correct_Ez_0(Grid1DTM &aux_grid);
    void correct_Ez_f(Grid1DTM &aux_grid);
    void correct_Ex_0(Grid1DTM &aux_grid);
    void correct_Ex_f(Grid1DTM &aux_grid);
    //Apply all the previous corrections
    void apply_TFSF(Grid1DTM &aux_grid);

    //Checks if a point is inside the area
    bool bounds_check(Point a);

};




#endif // DEFINITIONS_H_INCLUDED
