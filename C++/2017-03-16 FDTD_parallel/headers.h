#ifndef HEADERS_H_INCLUDED
#define HEADERS_H_INCLUDED

/*Include list*/
#include <iostream>
//#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <complex>

//Multithreading
#include <thread>

//Constants
const double PI = 3.14159265358979; //Pi constant
const double epsilon_0 = 8.854187e-12; //Permittivity of vacuum/free space
const double mu_0 = 4.0 * PI * 1.0e-7; //Permeability of free space
const double c0 = 299792458.0; //Speed of light in vacuum m/s
const double z0 = mu_0 * c0; //Characteristic impedance of free space Ohms
const unsigned int NUM_THREADS = std::thread::hardware_concurrency(); //Detects maximum number of threads

using namespace std;

//Dimensions of a dielectric square in number of nodes
struct Dielectric_sq {
    unsigned long length, width, height;
    double mu_rel, eps_rel;
    Dielectric_sq(); //Default constructor
    Dielectric_sq(unsigned long l, unsigned long w, unsigned long h, double mu, double eps); //Non-Default Constructor
};



struct Grid {
    unsigned long I, J, K, diel_x_start, diel_x_end, diel_y_start, diel_y_end, diel_z_start, diel_z_end;
    unsigned long long vsize;
    double wl, n_core, courant, abs_fact, ds, dt, frequency;
    vector<double> sigma, sigma_ast, mu, epsilon;
    vector<double> E_x, E_y, E_z;
    vector<double> H_x, H_y, H_z;
    vector<double> C_a, C_b, D_a, D_b;
    vector<double> E_y_x0, E_y_x1, E_z_x0, E_z_x1;
    vector<double> E_x_y0, E_x_y1, E_z_y0, E_z_y1;
    vector<double> E_x_z0, E_x_z1, E_y_z0, E_y_z1;


    Dielectric_sq def_obj {1,1,1,1.0,1.0}; //Default object to pass to the default constructor
    Dielectric_sq &obj;
    Grid(); //Default constructor
    Grid(Dielectric_sq &diel, double wll, unsigned long l, unsigned long w, unsigned long h); //Non-default constructor
    unsigned long long get_ind(const unsigned long long i, const unsigned long long j, const unsigned long long k);
    unsigned long long get_ind_2d(const unsigned long long a, const unsigned long long b, const unsigned long long L);
    void init_properties();
    void init_electric();
    void init_magnetic();
    void init_boundary();
    void update_boundary_x(unsigned long J, unsigned long K, double face, vector<double> &E_y_x0, vector<double> &E_y,
                        double abs_fact);
    void update_boundary_y(unsigned long J, unsigned long K, double face, vector<double> &E_y_x0, vector<double> &E_y,
                        double abs_fact);
    void update_boundary_z(unsigned long J, unsigned long K, double face, vector<double> &E_y_x0, vector<double> &E_y,
                        double abs_fact);
    void parallel_update_boundary_all();

    void update_E_x(unsigned long long I_s, unsigned long long I, unsigned long long J, unsigned long long K,
                    vector<double> &H_z, vector<double> &H_y, vector<double> &E_x, vector<double> &C_a, vector<double> &C_b);
    void update_E_y(unsigned long long I_s, unsigned long long I, unsigned long long J, unsigned long long K,
                    vector<double> &H_z, vector<double> &H_y, vector<double> &E_x, vector<double> &C_a, vector<double> &C_b);
    void update_E_z(unsigned long long I_s, unsigned long long I, unsigned long long J, unsigned long long K,
                    vector<double> &H_z, vector<double> &H_y, vector<double> &E_x, vector<double> &C_a, vector<double> &C_b);
    void parallel_update_E_field();

    void update_H_x(unsigned long long I_s, unsigned long long I, unsigned long long J, unsigned long long K,
                    vector<double> &H_z, vector<double> &H_y, vector<double> &E_x, vector<double> &C_a, vector<double> &C_b);
    void update_H_y(unsigned long long I_s, unsigned long long I, unsigned long long J, unsigned long long K,
                    vector<double> &H_z, vector<double> &H_y, vector<double> &E_x, vector<double> &C_a, vector<double> &C_b);
    void update_H_z(unsigned long long I_s, unsigned long long I, unsigned long long J, unsigned long long K,
                    vector<double> &H_z, vector<double> &H_y, vector<double> &E_x, vector<double> &C_a, vector<double> &C_b);
    void parallel_update_H_field();

    void serial_update_boundary_all();
    void serial_update_E_field();
    void serial_update_H_field();






};










#endif // HEADERS_H_INCLUDED

