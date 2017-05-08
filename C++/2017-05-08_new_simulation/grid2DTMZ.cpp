/*--------------Grid 2D TMz, #5--------------*/
#include "definitions.h"

Params::Params(double wavelengthh, unsigned long long max_timee, double ppww, unsigned long long size_xx,
           unsigned long long size_yy, int source_typee, Point tfsf_startt, Point tfsf_endd, int choicee)
: wavelength {wavelengthh}, max_time {max_timee}, size_x {size_xx}, size_y {size_yy}, source_type {source_typee},
tfsf_start {tfsf_startt}, tfsf_end {tfsf_endd}, choice {choicee}, ppw {ppww} {};

Grid2DTMZ::Grid2DTMZ(Params in)
:wavelength {in.wavelength}, source_type {in.source_type}, Hx{in.size_x, in.size_y}, Hy {in.size_x, in.size_y},
Ez {in.size_x, in.size_y}, c1hx {in.size_x, in.size_y}, c2hx {in.size_x, in.size_y}, c1hy {in.size_x, in.size_y},
c2hy {in.size_x, in.size_y}, c1ez {in.size_x, in.size_y}, c2ez {in.size_x, in.size_y}, tfsf_start {in.tfsf_start},
tfsf_end {in.tfsf_end}, max_time {in.max_time}, Ez_top {in.size_x, 6}, Ez_bottom {in.size_x, 6}, Ez_left {in.size_y, 6},
Ez_right {in.size_y, 6}, size_x {in.size_x}, size_y {in.size_y}, choice {in.choice}, ppw {in.ppw}
{
    min_wavelength = wavelength;
    dx = wavelength / ppw;
    dt = dx * courant / c0;

    //Initialize coeficients for update equations, free space, lossless
    for (unsigned long long i = 0; i < size_x; ++i) {
        for (unsigned long long j = 0; j< size_y; ++j) {
            c1ez(i,j) = 1.0;
            c2ez(i,j) = courant * imp0;
            c1hx(i,j) = 1.0;
            c2hx(i,j) = courant / imp0;
            c1hy(i,j) = 1.0;
            c2hy(i,j) = courant / imp0;
        }

    }
}

void Grid2DTMZ::show_params() {
    cout<<"Courant = "<<courant<<'\n';
    cout<<"Wavelength = "<<wavelength<<" metres"<<'\n';
    cout<<"Min Wavelength = "<<min_wavelength<<" metres"<<'\n';
    cout<<"ppw = "<<ppw<<" points per wavelength"<<'\n';
    cout<<"dx = "<<dx<<" metres"<<'\n';
    cout<<"dt = "<<dt<<" seconds"<<'\n';
    cout<<"time_delay = "<<time_delay<<" dt's\n";
    cout<<"dispersion = "<<dispersion<<" dt's\n";
    cout<<"======================================================================\n";
    return;
}

void Grid2DTMZ::update_magnetic() {
    for (unsigned long long i = 0; i < size_x; ++i) {
        for (unsigned long long j = 0; j < size_y - 1; ++j) {
            Hx(i,j) = c1hx(i,j) * Hx(i,j) - c2hx(i,j) * (Ez(i,j+1) - Ez(i,j));
        }
    }

    for (unsigned long long i = 0; i < size_x - 1; ++i) {
        for (unsigned long long j = 0; j < size_y; ++j) {
            Hy(i,j) = c1hy(i,j) * Hy(i,j) + c2hy(i,j) * (Ez(i+1,j) - Ez(i,j));
        }
    }

    return;
}

bool Grid2DTMZ::bounds_check(Point a) {
    if (a.x < size_x && a.y < size_y) {return true;}
    return false;
}

void Grid2DTMZ::update_electric() {
    for (unsigned long long i = 1; i < size_x - 1; ++i) {
        for (unsigned long long j = 1; j < size_y - 1; ++j) {
            Ez(i,j) = c1ez(i,j) * Ez(i,j) + c2ez(i,j) * ((Hy(i,j) - Hy(i-1,j)) - (Hx(i,j) - Hx(i,j-1)));
        }
    }

    return;
}

void Grid2DTMZ::abc_top() {
    unsigned long long i, j;
    //Notes:
    //Ez_top(i,0) = E0 at time t-1
    //Ez_top(i,1) = E1 at time t-1
    //Ez_top(i,2) = E2 at time t-1
    //Ez_top(i,3) = E0 at time t
    //Ez_top(i,4) = E1 at time t
    //Ez_top(i,5) = E2 at time t

    for(i = 0; i < size_x; ++i) {
        j = size_y - 1; //Top edge
        cour_prime = c2ez(i,j) * c2hy(i,j);
        A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
        B = (1.0/cour_prime) - 2 + cour_prime;
        C = 2.0 * (cour_prime - 1/cour_prime);
        D = -4.0 * (cour_prime + 1/cour_prime);
        Ez(i, j) = A * (B * (Ez(i, j-2) + Ez_top(i, 0)) + C * (Ez_top(i, 3) + Ez_top(i, 5) - Ez(i, j-1) - Ez_top(i, 1))
                + D * Ez_top(i, 4)) - Ez_top(i, 2);

        //Remember old fields, t-1 becomes t and t becomes t+1
        for (j = 0; j < 3; ++j) {
            Ez_top(i, j) = Ez_top(i, j+3);
            Ez_top(i, j+3) = Ez(i, size_y-1 - j);
        }
    }

    return;
}

void Grid2DTMZ::abc_bottom() {
    unsigned long long i, j;
    //Notes:
    //Ez_bottom(i,0) = E0 at time t-1
    //Ez_bottom(i,1) = E1 at time t-1
    //Ez_bottom(i,2) = E2 at time t-1
    //Ez_bottom(i,3) = E0 at time t
    //Ez_bottom(i,4) = E1 at time t
    //Ez_bottom(i,5) = E2 at time t

    for(i = 0; i < size_x; ++i) {
        j = 0; //bottom edge
        cour_prime = c2ez(i,j) * c2hy(i,j);
        A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
        B = (1.0/cour_prime) - 2 + cour_prime;
        C = 2.0 * (cour_prime - 1/cour_prime);
        D = -4.0 * (cour_prime + 1/cour_prime);
        Ez(i, j) = A * (B * (Ez(i, 2) + Ez_bottom(i, 0)) + C * (Ez_bottom(i, 3) + Ez_bottom(i, 5) - Ez(i, 1) - Ez_bottom(i, 1))
                + D * Ez_bottom(i, 4)) - Ez_bottom(i, 2);

        //Remember old fields, t-1 becomes t and t becomes t+1
        for (j = 0; j < 3; ++j) {
            Ez_bottom(i, j) = Ez_bottom(i, j+3);
            Ez_bottom(i, j+3) = Ez(i, j);
        }
    }

    return;
}

void Grid2DTMZ::abc_right() {
    unsigned long long i, j;
    //Notes:
    //Ez_right(i,0) = E0 at time t-1
    //Ez_right(i,1) = E1 at time t-1
    //Ez_right(i,2) = E2 at time t-1
    //Ez_right(i,3) = E0 at time t
    //Ez_right(i,4) = E1 at time t
    //Ez_right(i,5) = E2 at time t

    for(j = 0; j < size_y; ++j) {
        i = size_x - 1; //right edge
        cour_prime = c2ez(i,j) * c2hy(i,j);
        A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
        B = (1.0/cour_prime) - 2 + cour_prime;
        C = 2.0 * (cour_prime - 1/cour_prime);
        D = -4.0 * (cour_prime + 1/cour_prime);
        Ez(i, j) = A * (B * (Ez(i-2, j) + Ez_right(j, 0)) + C * (Ez_right(j, 3) + Ez_right(j, 5) - Ez(i-1, j) - Ez_right(j, 1))
                + D * Ez_right(j, 4)) - Ez_right(j, 2);

        //Remember old fields, t-1 becomes t and t becomes t+1
        for (i = 0; i < 3; ++i) {
            Ez_right(j, i) = Ez_right(j, i+3);
            Ez_right(j, i+3) = Ez(size_x - 1 - i, j);
        }
    }

    return;
}

void Grid2DTMZ::abc_left() {
    unsigned long long i, j;
    //Notes:
    //Ez_left(i,0) = E0 at time t-1
    //Ez_left(i,1) = E1 at time t-1
    //Ez_left(i,2) = E2 at time t-1
    //Ez_left(i,3) = E0 at time t
    //Ez_left(i,4) = E1 at time t
    //Ez_left(i,5) = E2 at time t

    for(j = 0; j < size_y; ++j) {
        i = 0; //left edge
        cour_prime = c2ez(i,j) * c2hy(i,j);
        A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
        B = (1.0/cour_prime) - 2 + cour_prime;
        C = 2.0 * (cour_prime - 1/cour_prime);
        D = -4.0 * (cour_prime + 1/cour_prime);
        Ez(i, j) = A * (B * (Ez(i+2, j) + Ez_left(j, 0)) + C * (Ez_left(j, 3) + Ez_left(j, 5) - Ez(i+1, j) - Ez_left(j, 1))
                + D * Ez_left(j, 4)) - Ez_left(j, 2);

        //Remember old fields, t-1 becomes t and t becomes t+1
        for (i = 0; i < 3; ++i) {
            Ez_left(j, i) = Ez_left(j, i+3);
            Ez_left(j, i+3) = Ez(i, j);
        }
    }

    return;
}

void Grid2DTMZ::apply_abc() {
    abc_top();
    abc_bottom();
    abc_left();
    abc_right();
    return;
}

void Grid2DTMZ::apply_TFSF(Grid1DTM &aux_grid) {
    unsigned long long i, j;
    //Correction along the left edge
    i = tfsf_start.x - 1;

    for (j = tfsf_start.y; j <= tfsf_end.y; ++j) {
        Hy(i, j) -= c2hy(i,j) * aux_grid.return_Ez(i+1 - tfsf_start.x+1);
    }

    //Correct Hy along the right edge
    i = tfsf_end.x;

    for (j = tfsf_start.y; j <= tfsf_end.y; ++j) {
        Hy(i, j) += c2hy(i,j) * aux_grid.return_Ez(i - tfsf_start.x+1);
    }

    //Correct Hx along the bottom
    j = tfsf_start.y - 1;

    for (i = tfsf_start.x; i <= tfsf_end.x; ++i) {
        Hx(i, j) += c2hx(i,j) * aux_grid.return_Ez(i - tfsf_start.x+1);
    }

    //Correct Hx along the top
    j = tfsf_end.y;

    for (i = tfsf_start.x; i <= tfsf_end.x; ++i) {
        Hx(i, j) -= c2hx(i,j) * aux_grid.return_Ez(i - tfsf_start.x+1);
    }


    //Advance one time step the auxiliary grid
    aux_grid.advance_simulation_hs();

    //Correct Ez field along left edge
    i = tfsf_start.x;
    for (j = tfsf_start.y; j <  tfsf_end.y; ++j) {
        Ez(i,j) -= c2ez(i,j) * aux_grid.return_Hy(i-1 - tfsf_start.x+1);
    }

    //Correct Ez field along right edge
    i = tfsf_end.x;
    for (j = tfsf_start.y; j < tfsf_end.y; ++j) {
        Ez(i,j) += c2ez(i,j) * aux_grid.return_Hy(i - tfsf_start.x+1);
    }

    //No need to correct Ez along top and bottom since incident Hx is zero
    return;

}


void Grid2DTMZ::apply_open_right_TFSF(Grid1DTM &aux_grid) {
    unsigned long long i, j;
    //Correction along the left edge
    i = tfsf_start.x - 1;

    for (j = tfsf_start.y; j <= tfsf_end.y; ++j) {
        Hy(i, j) -= c2hy(i,j) * aux_grid.return_Ez(i+1 - tfsf_start.x+1);
    }

    //No Hy correction along the right edge

    //Correct Hx along the bottom
    j = tfsf_start.y - 1;

    for (i = tfsf_start.x; i <= tfsf_end.x; ++i) {
        Hx(i, j) += c2hx(i,j) * aux_grid.return_Ez(i - tfsf_start.x+1);
    }

    //Correct Hx along the top
    j = tfsf_end.y;

    for (i = tfsf_start.x; i <= tfsf_end.x; ++i) {
        Hx(i, j) -= c2hx(i,j) * aux_grid.return_Ez(i - tfsf_start.x+1);
    }


    //Advance one time step the auxiliary grid
    aux_grid.advance_simulation_hs();

    //Correct Ez field along left edge
    i = tfsf_start.x;
    for (j = tfsf_start.y; j <  tfsf_end.y; ++j) {
        Ez(i,j) -= c2ez(i,j) * aux_grid.return_Hy(i-1 - tfsf_start.x+1);
    }

    //No Ez correction along right edge to allow the field to escape

    //No need to correct Ez along top and bottom since incident Hx is zero
    return;

}


void Grid2DTMZ::advance(Grid1DTM &aux_grid) {
    switch(choice) {
    case 0:
        {
        update_magnetic();
        apply_TFSF(aux_grid);
        update_electric();
        apply_abc();
        ++curr_time;
        return;
        }
    case 1:
        {
        update_magnetic();
        apply_open_right_TFSF(aux_grid);
        update_electric();
        apply_abc();
        ++curr_time;
        return;
        }
    default:
        {
        update_magnetic();
        apply_TFSF(aux_grid);
        update_electric();
        apply_abc();
        ++curr_time;
        return;
        }
    }
}

unsigned long long Grid2DTMZ::return_maximum_time() {
    return max_time;
}

unsigned long long Grid2DTMZ::return_current_time() {
    return curr_time;
}

Point Grid2DTMZ::size() {
    Point p {size_x, size_y};
    return p;
}

void Grid2DTMZ::save_state() {
    stringstream name;
    name<<"0-"<<"snapshot2D_TMZ_"<<curr_time<<".csv";
    ofstream file {name.str()};

    for (unsigned long long x = 0; x < size_x; ++x) {
        file<<'\n';
        for (unsigned long long y = 0; y < size_y; ++y) {
            file<<x<<','<<y<<','<<Ez(x, y)<<','<<Hx(x, y)<<','<<Hy(x, y)<<'\n';
        }
    }
    file.close();
    return;
}

void Grid2DTMZ::add_diel_rectangle(Point start, Point finish, double eps_r, double sigma, double mu_r, double sigma_m) {
    double new_min_wl = wavelength / (sqrt(eps_r * mu_r));
    if (bounds_check(start) && bounds_check(finish)) {
        if (new_min_wl < min_wavelength) {
            min_wavelength = new_min_wl;
            double new_dx = (min_wavelength / ppw);
            dx = new_dx;
            ppw = wavelength / dx;
            dt = courant * dx / c0;
            cout<<"Warning following values changed:\n";
            cout<<"min_wavelength = "<<min_wavelength<<'\n';
            cout<<"dx = "<<dx<<'\n';
            cout<<"dt = "<<dt<<'\n';
            cout<<"ppw = "<<ppw<<'\n';
        }

        for (unsigned long long i = start.x; i < finish.x; ++i) {
            for (unsigned long long j = start.y; j < finish.y; ++j) {
                c1ez(i,j) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                c2ez(i,j) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                c1hx(i,j) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                c2hx(i,j) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                c1hy(i,j) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                c2hy(i,j) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
            }
        }

        return;
    } else {
        cout<<"Function add_diel_rectangle() failed, point out of bounds\nNo object added\n";
        return;
    }
}

void Grid2DTMZ::add_diel_circle(Point center, unsigned long long radius, double eps_r, double sigma, double mu_r,
                         double sigma_m) {
    Point test {center.x + radius, center.y + radius};
    Point test2 {center.x - radius, center.y - radius};
    if (bounds_check(test) && bounds_check(test2)) {
        double new_min_wl = wavelength / (sqrt(eps_r * mu_r));

        if (new_min_wl < min_wavelength) {
            min_wavelength = new_min_wl;
            double new_dx = (min_wavelength / ppw);
            dx = new_dx;
            ppw = wavelength / dx;
            dt = courant * dx / c0;
            cout<<"Warning following values changed:\n";
            cout<<"min_wavelength = "<<min_wavelength<<'\n';
            cout<<"dx = "<<dx<<'\n';
            cout<<"dt = "<<dt<<'\n';
            cout<<"ppw = "<<ppw<<'\n';
        }

        for (unsigned long long i = center.x - radius; i < center.x + radius; ++i) {
            for (unsigned long long j = center.y - radius; j < center.y + radius; ++j) {
                double x = i;
                double y = j;
                double val1 = pow(x - center.x, 2.0);
                double val2 = pow(y - center.y, 2.0);
                double rad_sq = pow(radius, 2.0);

                if (val1 + val2 < rad_sq) {
                    c1ez(i,j) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c2ez(i,j) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c1hx(i,j) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c2hx(i,j) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c1hy(i,j) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c2hy(i,j) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                }
            }
        }

        return;

    } else {
        cout<<"Function add_diel_circle() failed, point out of bounds\nNo object added\n";
        return;
    }
}

void Grid2DTMZ::add_pec_circle(Point center, unsigned long long radius) {
    Point test {center.x + radius, center.y + radius};
    Point test2 {center.x - radius, center.y - radius};
    if (bounds_check(test) && bounds_check(test2)) {
        for (unsigned long long i = center.x - radius; i < center.x + radius; ++i) {
            for (unsigned long long j = center.y - radius; j < center.y + radius; ++j) {
                double x = i;
                double y = j;
                double val1 = pow(x - center.x, 2.0);
                double val2 = pow(y - center.y, 2.0);
                double rad_sq = pow(radius, 2.0);

                if (val1 + val2 < rad_sq) {
                    c1ez(i,j) = 0.0;
                    c2ez(i,j) = 0.0;
                }
            }
        }
        return;

    } else {
        cout<<"Function add_pec_circle() failed, point out of bounds\nNo object added\n";
        return;
    }
}

void Grid2DTMZ::add_pec_rectangle(Point start, Point finish) {
    if (bounds_check(start) && bounds_check(finish)) {
        for (unsigned long long i = start.x; i < finish.x; ++i) {
            for (unsigned long long j = start.y; j < finish.y; ++j) {
                c1ez(i,j) = 0.0;
                c2ez(i,j) = 0.0;
            }
        }
        return;

    } else {
        cout<<"Function add_pec_rectangle() failed, point out of bounds\nNo object added\n";
        return;
    }

}

unsigned long long convert_from_metres(double value, double wavelength, double eps_rel_max, double mu_rel_max, double ppw) {

        double new_min_wl = wavelength / (sqrt(eps_rel_max * mu_rel_max));
        double courant = 1.0 / sqrt(2.0);
        double new_dx = (new_min_wl / ppw);
        ppw = wavelength / new_dx;
        unsigned long long result = value / new_dx;
        return result;
}
