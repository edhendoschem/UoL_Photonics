#include "headers.h"
//Pending: Incorporate source and file writing

using namespace std;

int main()
{

    Dielectric_sq waveguide {80, 20, 20, 1.0, 4.0}; //lxwxh in nodes, mu_rel, eps_rel
    Grid grid {waveguide, 1.55e-6, 140, 40, 40};


    cout<<"E_x.size() = "<<grid.E_x.size()<<'\n';
    cout<<"ds = "<<grid.ds<<'\n';
    cout<<"dt = "<<grid.dt<<'\n';
    cout<<"wl = "<<grid.wl<<'\n';
    cout<<"courant = "<<grid.courant<<'\n';
    cout<<"n_core = "<<grid.n_core<<'\n';
    cout<<"abs_fact = "<<grid.abs_fact<<'\n';

    double source {0.0};
    unsigned long long n_max {10000};
    unsigned long long diel_cent_x = grid.diel_x_start + (grid.diel_x_end - grid.diel_x_start) / 2;
    unsigned long long diel_cent_y = grid.diel_y_start + (grid.diel_y_end - grid.diel_y_start) / 2;
    unsigned long long diel_cent_z = grid.diel_z_start + (grid.diel_z_end - grid.diel_z_start) / 2;

    for (unsigned long long n = 0; n < n_max; ++n) {
        cout<<"Time loop = "<<n<<'\n';
        //bad sources, fix
        /*
        //Additive source
        source = sin(2.0 * PI * grid.frequency * grid.dt * (n+1));
        grid.H_y[grid.get_ind(5, diel_cent_y, diel_cent_z)] =
                                grid.H_y[grid.get_ind(5, diel_cent_y, diel_cent_z)] + source;
        grid.H_z[grid.get_ind(5, diel_cent_y, diel_cent_z)] =
                                grid.H_z[grid.get_ind(5, diel_cent_y, diel_cent_z)] + source;

        source = cos(2.0 * PI * grid.frequency * grid.dt * (n+1));
        grid.E_x[grid.get_ind(5, diel_cent_y, diel_cent_z)] =
                                grid.E_x[grid.get_ind(5, diel_cent_y, diel_cent_z)] + source;
        */
        //Hard source
        source = sin(2.0 * PI * grid.frequency * grid.dt * (n+1));
        grid.H_y[grid.get_ind(5, diel_cent_y, diel_cent_z)] = source;
        grid.H_z[grid.get_ind(5, diel_cent_y, diel_cent_z)] = source;

        source = cos(2.0 * PI * grid.frequency * grid.dt * (n+1));
        grid.E_x[grid.get_ind(5, diel_cent_y, diel_cent_z)] = source;

        //serial processing
        //grid.serial_update_H_field();
        //grid.serial_update_E_field();
        //grid.serial_update_boundary_all();

        //parallel processing
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

    ofstream file1b;
    file1b.open("im2.csv");
    file1b <<"#x, y, z, E_x, E_y, E_z, H_x, H_y, H_z\n";
    for (long i = 0; i < grid.I; i+=1) {
        if (i>0) {
            file1b<<'\n';
        }
        for (long j = 0; j < grid.J; j+=1) {
                for (long k = 0; k < grid.K; k+=1)
            file1b<<i * grid.ds * 1e6<<','<<j * grid.ds * 1e6<<','<<k * grid.ds * 1e6<<','<<grid.E_x[grid.get_ind(i, j, k)]<<
            ','<<grid.E_y[grid.get_ind(i, j, k)]<<','<<grid.E_z[grid.get_ind(i, j, k)]<<','<<grid.H_x[grid.get_ind(i, j, k)]<<
            ','<<grid.H_y[grid.get_ind(i, j, k)]<<','<<grid.H_z[grid.get_ind(i, j, k)]<<'\n';
        }
    }

    file1b.close();


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

    ofstream info {"info.txt"};
    info << "Steps info\n";
    info << "ds = "<<grid.ds * 1.0e6<<" um\tdt = "<<grid.dt * 1.0e12<<" ps\n";
    info << "ds = "<<grid.ds<<" m\tdt = "<<grid.dt<<" s\n";
    info << "Total time steps = "<<n_max<<"\ttime elapsed = "<<n_max * grid.dt * 1.0e12<<" ps\n";
    info << "__________________________________________________________\n";
    info << "Dimensions info\n";
    info << "Grid dimensions: x = "<<grid.I<<"\ty = "<<grid.J<<"\tz = "<<grid.K<<'\n';
    info << "Grid dimensions: x = "<<grid.I * grid.ds * 1.0e6<<" um\ty = "<<grid.J * grid.ds * 1.0e6
    <<" um\tz = "<<grid.K * grid.ds * 1.0e6<<" um\n";
    info << "Waveguide start x = "<<grid.diel_x_start<<"\ty = "<<grid.diel_y_start<<"\tz = "<<grid.diel_z_start<<'\n';
    info << "Waveguide end x = "<<grid.diel_x_end<<"\ty = "<<grid.diel_y_end<<"\tz = "<<grid.diel_z_end<<'\n';
    info << "Waveguide start x = "<<grid.diel_x_start * grid.ds * 1.0e6<<" um\ty = "<<grid.diel_y_start * grid.ds * 1.0e6
    <<" um\tz = "<<grid.diel_z_start * grid.ds * 1.0e6<<" um\n";
    info << "Waveguide end x = "<<grid.diel_x_end * grid.ds * 1.0e6<<" um\ty = "<<grid.diel_y_end * grid.ds * 1.0e6
    <<" um\tz = "<<grid.diel_z_end * grid.ds * 1.0e6<<" um\n";

    info.close();


    return 0;
}
