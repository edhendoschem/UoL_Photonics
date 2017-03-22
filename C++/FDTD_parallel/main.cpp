#include "headers.h"
//Pending: Incorporate source and file writing

using namespace std;

int main()
{

    Dielectric_sq waveguide {80, 20, 20, 1.0, 4.0}; //lxwxh in nodes, mu_rel, eps_rel
    Grid grid {waveguide, 1.55e-6, 100, 40, 40};


    cout<<"E_x.size() = "<<grid.E_y[2222]<<'\n';
    cout<<"ds = "<<grid.ds<<'\n';
    cout<<"dt = "<<grid.dt<<'\n';
    cout<<"wl = "<<grid.wl<<'\n';
    cout<<"courant = "<<grid.courant<<'\n';
    cout<<"n_core = "<<grid.n_core<<'\n';
    cout<<"abs_fact = "<<grid.abs_fact<<'\n';

    double source {0.0};

    for (long n = 0; n < 5000; ++n) {
        cout<<"Time loop = "<<n<<'\n';
        //bad source
        source = sin(2.0 * PI * grid.frequency * grid.dt * (n+1));
        grid.H_y[grid.get_ind(grid.diel_x_start, grid.diel_y_start+10, grid.diel_z_start+10)] =
                                grid.H_y[grid.get_ind(grid.diel_x_start, grid.diel_y_start+10, grid.diel_z_start+10)] + source;


        //serial processing 103 sec
        //grid.serial_update_H_field();
        //grid.serial_update_E_field();
        //grid.serial_update_boundary_all();

        //parallel processing 69 sec
        grid.parallel_update_H_field();
        grid.parallel_update_E_field();
        grid.parallel_update_boundary_all();


    }


    //z = middle of waveguide, E_x
    long long source_z = grid.diel_z_start + (grid.diel_z_end - grid.diel_z_start) / 2.0;

    cout<<"Writing output file\n";

    ofstream file;
    file.open("im.csv");
    file <<"#x, y, z, E_x, E_y, E_z, H_x, H_y, H_z\n";
    for (long i = 0; i < grid.I; i+=1) {
        if (i>0) {
            file<<'\n';
        }
        for (long j = 0; j < grid.J; j+=1) {
                for (long k = 0; k < grid.K; k+=1)
            file<<i<<','<<j<<','<<k<<','<<grid.E_x[grid.get_ind(i, j, k)]<<','<<grid.E_y[grid.get_ind(i, j, k)]<<
            ','<<grid.E_z[grid.get_ind(i, j, k)]<<','<<grid.H_x[grid.get_ind(i, j, k)]<<','<<grid.H_y[grid.get_ind(i, j, k)]<<
            ','<<grid.H_z[grid.get_ind(i, j, k)]<<'\n';
        }
    }

    file.close();


    ofstream file1a;
    file1a.open("imex.csv");
    file1a <<"#x, y, E_x\n";
    for (long i = 0; i < grid.I; i+=1) {
        for (long j = 0; j < grid.J; j+=1) {
            file1a<<i<<','<<j<<','<<grid.E_x[grid.get_ind(i, j, source_z)]<<'\n';
        }
    }

    file1a.close();

    //z = middle of waveguide, E_y
    ofstream file2;
    file2.open("imey.csv");
    file2 <<"#x, y, E_y\n";
    for (long i = 0; i < grid.I; i+=1) {
        for (long j = 0; j < grid.J; j+=1) {

            file2<<i<<','<<j<<','<<grid.E_y[grid.get_ind(i, j, source_z)]<<'\n';
        }
    }

    ofstream file2a;
    file2a.open("imeya.csv");
    file2a <<"#x, z, E_y\n";
    for (long i = 0; i < grid.I; i+=1) {
        for (long k = 0; k < grid.K; k+=1) {

            file2a<<i<<','<<k<<','<<grid.E_y[grid.get_ind(i, source_z, k)]<<'\n';
        }
    }

    file2a.close();

    //z = middle of waveguide, E_z
    ofstream file3;
    file3.open("imez.csv");
    file3 <<"#x, y, E_z\n";
    for (long i = 0; i < grid.I; i+=1) {
        for (long j = 0; j < grid.J; j+=1) {
            file3<<i<<','<<j<<','<<grid.E_z[grid.get_ind(i, j, source_z)]<<'\n';
        }
    }

    file3.close();




    return 0;
}
