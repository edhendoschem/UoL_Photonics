/*--------------Grid 3D, #7--------------*/
#include "definitions.h"

Grid3D::Grid3D(Params in) : wavelength {in.wavelength}, source_type {in.source_type}, Hx{in.size_x, in.size_y, in.size_z},
Hy {in.size_x, in.size_y, in.size_z}, Hz {in.size_x, in.size_y, in.size_z}, Ex {in.size_x, in.size_y, in.size_z},
Ey {in.size_x, in.size_y, in.size_z}, Ez {in.size_x, in.size_y, in.size_z}, c1hx {in.size_x, in.size_y, in.size_z},
c2hx {in.size_x, in.size_y, in.size_z}, c1hy {in.size_x, in.size_y, in.size_z}, c2hy {in.size_x, in.size_y, in.size_z},
c1hz {in.size_x, in.size_y, in.size_z}, c2hz {in.size_x, in.size_y, in.size_z}, c1ex {in.size_x, in.size_y, in.size_z},
c2ex {in.size_x, in.size_y, in.size_z}, c1ey {in.size_x, in.size_y, in.size_z}, c2ey {in.size_x, in.size_y, in.size_z},
c1ez {in.size_x, in.size_y, in.size_z}, c2ez {in.size_x, in.size_y, in.size_z}, tfsf_start {in.tfsf_start},
tfsf_end {in.tfsf_end}, max_time {in.max_time}, size_x {in.size_x}, size_y {in.size_y}, size_z {in.size_z},
choice {in.choice}, ppw {in.ppw}, Ey_x0 {in.size_y, in.size_z, 6}, Ez_x0 {in.size_y, in.size_z, 6},
Ey_xf {in.size_y, in.size_z, 6}, Ez_xf {in.size_y, in.size_z, 6}, Ex_y0 {in.size_x, in.size_z, 6},
Ez_y0 {in.size_x, in.size_z, 6}, Ex_yf {in.size_x, in.size_z, 6}, Ez_yf {in.size_x, in.size_z, 6},
Ex_z0 {in.size_x, in.size_y, 6}, Ey_z0 {in.size_x, in.size_y, 6}, Ex_zf {in.size_x, in.size_y, 6},
Ey_zf {in.size_x, in.size_y, 6} {

    min_wavelength = wavelength;
    dx = wavelength / ppw;
    dt = dx * courant / c0;
    for (unsigned long long i = 0; i < size_x; ++i) {
        for (unsigned long long j = 0; j< size_y; ++j) {
            for (unsigned long long k = 0; k < size_z; ++k) {
                //Electric field coefficients
                c1ex(i,j,k) = 1.0;
                c2ex(i,j,k) = courant * imp0;
                c1ey(i,j,k) = 1.0;
                c2ey(i,j,k) = courant * imp0;
                c1ez(i,j,k) = 1.0;
                c2ez(i,j,k) = courant * imp0;
                //Magnetic field coefficients
                c1hx(i,j,k) = 1.0;
                c2hx(i,j,k) = courant / imp0;
                c1hy(i,j,k) = 1.0;
                c2hy(i,j,k) = courant / imp0;
                c1hz(i,j,k) = 1.0;
                c2hz(i,j,k) = courant / imp0;

            }
        }
    }

}

void Grid3D::save_state() {
    stringstream name;
    name<<"0-"<<"snapshot3D_"<<curr_time<<".csv";
    ofstream file {name.str()};
    for (unsigned long long x = 0; x < size_x; ++x) {
        file<<'\n';
        for (unsigned long long y = 0; y < size_y; ++y) {
            for (unsigned long long z = 0; z < size_z; ++z) {
                file<<x<<','<<y<<','<<z<<','<<Ex(x,y,z)<<','<<Ey(x,y,z)<<','<<Ez(x,y,z)<<','<<Hx(x,y,z)\
                <<','<<Hy(x,y,z)<<','<<Hz(x,y,z)<<'\n';
            }
        }
    }

    file.close();
    return;
}

void Grid3D::xy_cross_section(unsigned long long z, char w) {
    stringstream name;
    name<<"0-"<<"xy"<<"_snapshot3D_"<<w<<"_"<<curr_time<<".csv";
    ofstream file {name.str()};
    file<<"#z = "<<z<<'\n';
    for (unsigned long long x = 0; x < size_x; ++x) {
        file<<'\n';
        for (unsigned long long y = 0; y < size_y; ++y) {
            file<<x<<','<<y<<','<<Ex(x,y,z)<<','<<Ey(x,y,z)<<','<<Ez(x,y,z)<<','<<Hx(x,y,z)\
            <<','<<Hy(x,y,z)<<','<<Hz(x,y,z)<<'\n';
        }
    }

    file.close();
    return;
}

void Grid3D::yz_cross_section(unsigned long long x, char w) {
    stringstream name;
    name<<"0-"<<"yz"<<"_snapshot3D_"<<w<<"_"<<curr_time<<".csv";
    ofstream file {name.str()};
    file<<"#x = "<<x<<'\n';
    for (unsigned long long y = 0; y < size_y; ++y) {
        file<<'\n';
        for (unsigned long long z = 0; z < size_z; ++z) {
            file<<y<<','<<z<<','<<Ex(x,y,z)<<','<<Ey(x,y,z)<<','<<Ez(x,y,z)<<','<<Hx(x,y,z)\
            <<','<<Hy(x,y,z)<<','<<Hz(x,y,z)<<'\n';
        }
    }

    file.close();
    return;
}

void Grid3D::xz_cross_section(unsigned long long y, char w) {
    stringstream name;
    name<<"0-"<<"xz"<<"_snapshot3D_"<<w<<"_"<<curr_time<<".csv";
    ofstream file {name.str()};
    file<<"#y = "<<y<<'\n';
    for (unsigned long long x = 0; x < size_x; ++x) {
        file<<'\n';
        for (unsigned long long z = 0; z < size_z; ++z) {
            file<<x<<','<<z<<','<<Ex(x,y,z)<<','<<Ey(x,y,z)<<','<<Ez(x,y,z)<<','<<Hx(x,y,z)\
            <<','<<Hy(x,y,z)<<','<<Hz(x,y,z)<<'\n';
        }
    }

    file.close();
    return;
}

void Grid3D::test() {
    stringstream name;
    name<<"0-"<<"INFO"<<".csv";
    ofstream file(name.str(), std::ios_base::app);
    unsigned long long x, y, z;

    x = 42, y = 42, z = 25;
    file<<curr_time<<','<<x<<','<<y<<','<<Ex(x,y,z)<<','<<Ey(x,y,z)<<','<<Ez(x,y,z)<<','<<Hx(x,y,z)\
            <<','<<Hy(x,y,z)<<','<<Hz(x,y,z)<<'\n';
    x = 42, y = 25, z = 42;
    file<<curr_time<<','<<x<<','<<z<<','<<Ex(x,y,z)<<','<<Ey(x,y,z)<<','<<Ez(x,y,z)<<','<<Hx(x,y,z)\
            <<','<<Hy(x,y,z)<<','<<Hz(x,y,z)<<'\n';
    x = 25, y = 42, z = 42;
    file<<curr_time<<','<<y<<','<<z<<','<<Ex(x,y,z)<<','<<Ey(x,y,z)<<','<<Ez(x,y,z)<<','<<Hx(x,y,z)\
            <<','<<Hy(x,y,z)<<','<<Hz(x,y,z)<<'\n';

    file.close();
    return;
}

void Grid3D::update_Hx() {
    unsigned long long i, j, k;

    for (i = 0; i < size_x; ++i) {
        for (j = 0; j < size_y - 1; ++j) {
            for (k = 0; k < size_z - 1; ++k) {
                Hx(i,j,k) = c1hx(i,j,k) * Hx(i,j,k) + c2hx(i,j,k) * ((Ey(i,j,k+1) - Ey(i,j,k)) - (Ez(i,j+1,k) - Ez(i,j,k)));
            }
        }
    }

    return;
}

void Grid3D::update_Hy() {
    unsigned long long i, j, k;

    for (i = 0; i < size_x - 1; ++i) {
        for (j = 0; j < size_y; ++j) {
            for (k = 0; k < size_z - 1; ++k) {
                Hy(i,j,k) = c1hy(i,j,k) * Hy(i,j,k) + c2hy(i,j,k) * ((Ez(i+1,j,k) - Ez(i,j,k)) - (Ex(i,j,k+1) - Ex(i,j,k)));
            }
        }
    }

    return;
}

void Grid3D::update_Hz() {
    unsigned long long i, j, k;

    for (i = 0; i < size_x - 1; ++i) {
        for (j = 0; j < size_y - 1; ++j) {
            for (k = 0; k < size_z; ++k) {
                Hz(i,j,k) = c1hz(i,j,k) * Hz(i,j,k) + c2hz(i,j,k) * ((Ex(i,j+1,k) - Ex(i,j,k)) - (Ey(i+1,j,k) - Ey(i,j,k)));
            }
        }
    }

    return;
}

void Grid3D::update_magnetic() {
    update_Hx();
    update_Hy();
    update_Hz();
    return;
}

void Grid3D::update_Ex() {
    unsigned long long i, j, k;

    for (i = 0; i < size_x - 1; ++i) {
        for (j = 1; j < size_y - 1; ++j) {
            for (k = 1; k < size_z - 1; ++k) {
                Ex(i,j,k) = c1ex(i,j,k) * Ex(i,j,k) + c2ex(i,j,k) * ((Hz(i,j,k) - Hz(i,j-1,k)) - (Hy(i,j,k) - Hy(i,j,k-1)));
            }
        }
    }

    return;
}

void Grid3D::update_Ey() {
    unsigned long long i, j, k;

    for (i = 1; i < size_x - 1; ++i) {
        for (j = 0; j < size_y - 1; ++j) {
            for (k = 1; k < size_z - 1; ++k) {
                Ey(i,j,k) = c1ey(i,j,k) * Ey(i,j,k) + c2ey(i,j,k) * ((Hx(i,j,k) - Hx(i,j,k-1)) - (Hz(i,j,k) - Hz(i-1,j,k)));
            }
        }
    }

    return;
}

void Grid3D::update_Ez() {
    unsigned long long i, j, k;

    for (i = 1; i < size_x - 1; ++i) {
        for (j = 1; j < size_y - 1; ++j) {
            for (k = 0; k < size_z - 1; ++k) {
                Ez(i,j,k) = c1ez(i,j,k) * Ez(i,j,k) + c2ez(i,j,k) * ((Hy(i,j,k) - Hy(i-1,j,k)) - (Hx(i,j,k) - Hx(i,j-1,k)));
            }
        }
    }

    return;
}

void Grid3D::update_electric() {
    update_Ex();
    update_Ey();
    update_Ez();
}


void Grid3D::simple_abc() {
    simple_x0();
    simple_xf();
    simple_y0();
    simple_yf();
    simple_z0();
    simple_zf();
}

void Grid3D::correct_Hx_0(Grid1DTM &aux_grid) {
    unsigned long long i,j, k;
    j = tfsf_start.y;

    for (i = tfsf_start.x; i < tfsf_end.x + 1; ++i) {
        for (k = tfsf_start.z; k < tfsf_end.z; ++k) {
            Hx(i, j-1, k) += c2hx(i,j-1,k) * aux_grid.return_Ez(i - tfsf_start.x+1);
        }
    }

    return;
}

void Grid3D::correct_Hx_f(Grid1DTM &aux_grid) {
    unsigned long long i,j, k;
    j = tfsf_end.y;

    for (i = tfsf_start.x; i < tfsf_end.x + 1; ++i) {
        for (k = tfsf_start.z; k < tfsf_end.z; ++k) {
            Hx(i, j, k) -= c2hx(i,j,k) * aux_grid.return_Ez(i - tfsf_start.x+1);
        }
    }

    return;
}

void Grid3D::correct_Hy_0(Grid1DTM &aux_grid) {
    unsigned long long i,j, k;
    i = tfsf_start.x;

    for (j = tfsf_start.y; j < tfsf_end.y + 1; ++j) {
        for (k = tfsf_start.z; k < tfsf_end.z; ++k) {
            Hy(i-1, j, k) -= c2hy(i-1,j,k) * aux_grid.return_Ez(i - tfsf_start.x+1);
        }
    }

    return;
}

void Grid3D::correct_Hy_f(Grid1DTM &aux_grid) {
    unsigned long long i,j, k;
    i = tfsf_end.x;

    for (j = tfsf_start.y; j < tfsf_end.y + 1; ++j) {
        for (k = tfsf_start.z; k < tfsf_end.z; ++k) {
            Hy(i, j, k) += c2hy(i,j,k) * aux_grid.return_Ez(i - tfsf_start.x+1);
        }
    }

    return;
}

void Grid3D::correct_Ez_0(Grid1DTM &aux_grid) {
    unsigned long long i,j, k;
    i = tfsf_start.x;

    for (j = tfsf_start.y; j < tfsf_end.y + 1; ++j) {
        for (k = tfsf_start.z; k < tfsf_end.z; ++k) {
            Ez(i, j, k) -= c2ez(i,j,k) * aux_grid.return_Hy(i - 1 - tfsf_start.x+1);
        }
    }

    return;
}

void Grid3D::correct_Ez_f(Grid1DTM &aux_grid) {
    unsigned long long i,j, k;
    i = tfsf_end.x;

    for (j = tfsf_start.y; j < tfsf_end.y + 1; ++j) {
        for (k = tfsf_start.z; k < tfsf_end.z; ++k) {
            Ez(i, j, k) += c2ez(i,j,k) * aux_grid.return_Hy(i - tfsf_start.x+1);
        }
    }

    return;
}

void Grid3D::correct_Ex_0(Grid1DTM &aux_grid) {
    unsigned long long i,j, k;
    k = tfsf_start.z;

    for (i = tfsf_start.x; i < tfsf_end.x; ++i) {
        for (j = tfsf_start.y; j < tfsf_end.y + 1; ++j) {
            Ex(i, j, k) += c2ex(i,j,k) * aux_grid.return_Hy(i - tfsf_start.x+1);
        }
    }

    return;
}

void Grid3D::correct_Ex_f(Grid1DTM &aux_grid) {
    unsigned long long i,j, k;
    k = tfsf_end.z;

    for (i = tfsf_start.x; i < tfsf_end.x; ++i) {
        for (j = tfsf_start.y; j < tfsf_end.y + 1; ++j) {
            Ex(i, j, k) -= c2ex(i,j,k) * aux_grid.return_Hy(i - tfsf_start.x+1);
        }
    }

    return;
}

void Grid3D::apply_TFSF(Grid1DTM &aux_grid) {
    correct_Hy_0(aux_grid);
    correct_Hy_f(aux_grid);
    correct_Hx_0(aux_grid);
    correct_Hx_f(aux_grid);
    aux_grid.advance_simulation_hs();
    correct_Ez_0(aux_grid);
    correct_Ez_f(aux_grid);
    correct_Ex_0(aux_grid);
    correct_Ex_f(aux_grid);

}


void Grid3D::update_test(Grid1DTM &aux_grid) {
    update_magnetic();
    //Ez(25, 25, 25) = sin(2.0 * PI * (courant * (curr_time-time_delay) - 12.0) / ppw);
    //Hx(25, 25, 25) = sin(2.0 * PI * (courant * (curr_time-time_delay) - 12.0) / ppw);
    //Hy(25, 25, 25) = sin(2.0 * PI * (courant * (curr_time-time_delay) - 12.0) / ppw);
    //Ez(25, 25, 25) = exp(-pow(((curr_time-time_delay) - 2.2 * dispersion)/ dispersion, 2));
    //Hx(25, 25, 25) = exp(-pow(((curr_time-time_delay) - 2.2 * dispersion)/ dispersion, 2));
    //Hy(25, 25, 25) = exp(-pow(((curr_time-time_delay) - 2.2 * dispersion)/ dispersion, 2));
    apply_TFSF(aux_grid);
    update_electric();
    //apply_abc();
    simple_abc();
    ++curr_time;
}

