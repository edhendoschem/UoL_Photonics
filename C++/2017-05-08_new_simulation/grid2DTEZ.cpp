/*--------------Grid 2D TEz, #6--------------*/
#include "definitions.h"

Grid2DTEZ::Grid2DTEZ(Params in)
:wavelength {in.wavelength}, source_type {in.source_type}, Hz{in.size_x, in.size_y}, Ey {in.size_x, in.size_y},
Ex {in.size_x, in.size_y}, c1hz {in.size_x, in.size_y}, c2hz {in.size_x, in.size_y}, c1ey {in.size_x, in.size_y},
c2ey {in.size_x, in.size_y}, c1ex {in.size_x, in.size_y}, c2ex {in.size_x, in.size_y}, tfsf_start {in.tfsf_start},
tfsf_end {in.tfsf_end}, max_time {in.max_time}, Ex_top {in.size_x, 6}, Ex_bottom {in.size_x, 6}, Ey_left {in.size_y, 6},
Ey_right {in.size_y, 6}, size_x {in.size_x}, size_y {in.size_y}, choice {in.choice}, ppw {in.ppw}
{
    min_wavelength = wavelength;
    dx = wavelength / ppw;
    dt = dx * courant / c0;

    //Initialize coeficients for update equations, free space, lossless
    for (unsigned long long i = 0; i < size_x; ++i) {
        for (unsigned long long j = 0; j< size_y; ++j) {
            c1ey(i,j) = 1.0;
            c2ey(i,j) = courant * imp0;
            c1ex(i,j) = 1.0;
            c2ex(i,j) = courant * imp0;
            c1hz(i,j) = 1.0;
            c2hz(i,j) = courant / imp0;
        }

    }
}

void Grid2DTEZ::show_params() {
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

void Grid2DTEZ::update_magnetic() {
    for (unsigned long long i = 0; i < size_x-1; ++i) {
        for (unsigned long long j = 0; j < size_y - 1; ++j) {
            Hz(i,j) = c1hz(i,j) * Hz(i,j) + c2hz(i,j) * ((Ex(i,j+1) - Ex(i,j))-(Ey(i+1,j) - Ey(i,j)));
        }
    }


    return;
}

bool Grid2DTEZ::bounds_check(Point a) {
    if (a.x < size_x && a.y < size_y) {return true;}
    return false;
}

void Grid2DTEZ::update_electric() {
    for (unsigned long long i = 0; i < size_x - 1; ++i) {
        for (unsigned long long j = 1; j < size_y - 1; ++j) {
            Ex(i,j) = c1ex(i,j) * Ex(i,j) + c2ex(i,j) * ((Hz(i,j) - Hz(i,j-1)));
        }
    }

    for (unsigned long long i = 1; i < size_x - 1; ++i) {
        for (unsigned long long j = 0; j < size_y - 1; ++j) {
            Ey(i,j) = c1ey(i,j) * Ey(i,j) - c2ey(i,j) * ((Hz(i,j) - Hz(i-1,j)));
        }
    }

    return;
}

void Grid2DTEZ::abc_top() {
    unsigned long long i, j;
    //Notes:
    //Ez_top(i,0) = E0 at time t-1
    //Ez_top(i,1) = E1 at time t-1
    //Ez_top(i,2) = E2 at time t-1
    //Ez_top(i,3) = E0 at time t
    //Ez_top(i,4) = E1 at time t
    //Ez_top(i,5) = E2 at time t

    for(i = 0; i < size_x - 1; ++i) {
        j = size_y - 1; //Top edge
        cour_prime = c2hz(i,j) * c2ex(i,j);
        A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
        B = (1.0/cour_prime) - 2 + cour_prime;
        C = 2.0 * (cour_prime - 1/cour_prime);
        D = -4.0 * (cour_prime + 1/cour_prime);
        Ex(i, j) = A * (B * (Ex(i, j-2) + Ex_top(i, 0)) + C * (Ex_top(i, 3) + Ex_top(i, 5) - Ex(i, j-1) - Ex_top(i, 1))
                + D * Ex_top(i, 4)) - Ex_top(i, 2);

        //Remember old fields, t-1 becomes t and t becomes t+1
        for (j = 0; j < 3; ++j) {
            Ex_top(i, j) = Ex_top(i, j+3);
            Ex_top(i, j+3) = Ex(i, size_y-1 - j);
        }
    }

    return;
}

void Grid2DTEZ::abc_bottom() {
    unsigned long long i, j;
    //Notes:
    //Ez_bottom(i,0) = E0 at time t-1
    //Ez_bottom(i,1) = E1 at time t-1
    //Ez_bottom(i,2) = E2 at time t-1
    //Ez_bottom(i,3) = E0 at time t
    //Ez_bottom(i,4) = E1 at time t
    //Ez_bottom(i,5) = E2 at time t

    for(i = 0; i < size_x - 1; ++i) {
        j = 0; //bottom edge
        cour_prime = c2hz(i,j) * c2ex(i,j);
        A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
        B = (1.0/cour_prime) - 2 + cour_prime;
        C = 2.0 * (cour_prime - 1/cour_prime);
        D = -4.0 * (cour_prime + 1/cour_prime);
        Ex(i, j) = A * (B * (Ex(i, 2) + Ex_bottom(i, 0)) + C * (Ex_bottom(i, 3) + Ex_bottom(i, 5) - Ex(i, 1) - Ex_bottom(i, 1))
                + D * Ex_bottom(i, 4)) - Ex_bottom(i, 2);

        //Remember old fields, t-1 becomes t and t becomes t+1
        for (j = 0; j < 3; ++j) {
            Ex_bottom(i, j) = Ex_bottom(i, j+3);
            Ex_bottom(i, j+3) = Ex(i, j);
        }
    }

    return;
}

void Grid2DTEZ::abc_right() {
    unsigned long long i, j;
    //Notes:
    //Ez_right(i,0) = E0 at time t-1
    //Ez_right(i,1) = E1 at time t-1
    //Ez_right(i,2) = E2 at time t-1
    //Ez_right(i,3) = E0 at time t
    //Ez_right(i,4) = E1 at time t
    //Ez_right(i,5) = E2 at time t

    for(j = 0; j < size_y - 1; ++j) {
        i = size_x - 1; //right edge
        cour_prime = c2hz(i,j) * c2ey(i,j);
        A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
        B = (1.0/cour_prime) - 2 + cour_prime;
        C = 2.0 * (cour_prime - 1/cour_prime);
        D = -4.0 * (cour_prime + 1/cour_prime);
        Ey(i, j) = A * (B * (Ey(i-2, j) + Ey_right(j, 0)) + C * (Ey_right(j, 3) + Ey_right(j, 5) - Ey(i-1, j) - Ey_right(j, 1))
                + D * Ey_right(j, 4)) - Ey_right(j, 2);

        //Remember old fields, t-1 becomes t and t becomes t+1
        for (i = 0; i < 3; ++i) {
            Ey_right(j, i) = Ey_right(j, i+3);
            Ey_right(j, i+3) = Ey(size_x - 1 - i, j);
        }
    }

    return;
}

void Grid2DTEZ::abc_left() {
    unsigned long long i, j;
    //Notes:
    //Ez_left(i,0) = E0 at time t-1
    //Ez_left(i,1) = E1 at time t-1
    //Ez_left(i,2) = E2 at time t-1
    //Ez_left(i,3) = E0 at time t
    //Ez_left(i,4) = E1 at time t
    //Ez_left(i,5) = E2 at time t

    for(j = 0; j < size_y - 1; ++j) {
        i = 0; //left edge
        cour_prime = c2hz(i,j) * c2ey(i,j);
        A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
        B = (1.0/cour_prime) - 2 + cour_prime;
        C = 2.0 * (cour_prime - 1/cour_prime);
        D = -4.0 * (cour_prime + 1/cour_prime);
        Ey(i, j) = A * (B * (Ey(i+2, j) + Ey_left(j, 0)) + C * (Ey_left(j, 3) + Ey_left(j, 5) - Ey(i+1, j) - Ey_left(j, 1))
                + D * Ey_left(j, 4)) - Ey_left(j, 2);

        //Remember old fields, t-1 becomes t and t becomes t+1
        for (i = 0; i < 3; ++i) {
            Ey_left(j, i) = Ey_left(j, i+3);
            Ey_left(j, i+3) = Ey(i, j);
        }
    }

    return;
}

void Grid2DTEZ::apply_abc() {
    abc_top();
    abc_bottom();
    abc_left();
    abc_right();
    return;
}

void Grid2DTEZ::apply_TFSF(Grid1DTE &aux_grid) {
    unsigned long long i, j;
    //Correction along the left edge
    i = tfsf_start.x - 1;

    for (j = tfsf_start.y; j < tfsf_end.y; ++j) {
        Hz(i, j) += c2hz(i,j) * aux_grid.return_Ey(i+1 - tfsf_start.x+1);
    }

    //Correct Hz along the right edge
    i = tfsf_end.x;

    for (j = tfsf_start.y; j < tfsf_end.y; ++j) {
        Hz(i, j) -= c2hz(i,j) * aux_grid.return_Ey(i - tfsf_start.x+1);
    }

    //Advance 1 time step in the auxiliary grid
    aux_grid.advance_simulation_hs();

    //Correct Ex along the bottom
    j = tfsf_start.y;

    for (i = tfsf_start.x; i < tfsf_end.x; ++i) {
        Ex(i, j) -= c2ex(i,j) * aux_grid.return_Hz(i - tfsf_start.x+1);
    }

    //Correct Ex along the top
    j = tfsf_end.y;

    for (i = tfsf_start.x; i < tfsf_end.x; ++i) {
        Ex(i, j) += c2ex(i,j) * aux_grid.return_Hz(i - tfsf_start.x+1);
    }

    //Correct Ey field along left edge
    i = tfsf_start.x;
    for (j = tfsf_start.y; j <  tfsf_end.y; ++j) {
        Ey(i,j) += c2ey(i,j) * aux_grid.return_Hz(i-1 - tfsf_start.x+1);
    }

    //Correct Ey field along right edge
    i = tfsf_end.x;
    for (j = tfsf_start.y; j < tfsf_end.y; ++j) {
        Ey(i,j) -= c2ey(i,j) * aux_grid.return_Hz(i - tfsf_start.x+1);
    }

    return;
}


void Grid2DTEZ::apply_open_right_TFSF(Grid1DTE &aux_grid) {
    unsigned long long i, j;
    //Correction along the left edge
    i = tfsf_start.x - 1;

    for (j = tfsf_start.y; j < tfsf_end.y; ++j) {
        Hz(i, j) += c2hz(i,j) * aux_grid.return_Ey(i+1 - tfsf_start.x+1);
    }

    //No Correction Hz along the right edge

    //Advance 1 time step in the auxiliary grid
    aux_grid.advance_simulation_hs();

    //Correct Ex along the bottom
    j = tfsf_start.y;

    for (i = tfsf_start.x; i < tfsf_end.x; ++i) {
        Ex(i, j) -= c2ex(i,j) * aux_grid.return_Hz(i - tfsf_start.x+1);
    }

    //Correct Ex along the top
    j = tfsf_end.y;

    for (i = tfsf_start.x; i < tfsf_end.x; ++i) {
        Ex(i, j) += c2ex(i,j) * aux_grid.return_Hz(i - tfsf_start.x+1);
    }

    //Correct Ey field along left edge
    i = tfsf_start.x;
    for (j = tfsf_start.y; j <  tfsf_end.y; ++j) {
        Ey(i,j) += c2ey(i,j) * aux_grid.return_Hz(i-1 - tfsf_start.x+1);
    }

    //No Correction Ey field along right edge


    return;

}


void Grid2DTEZ::advance(Grid1DTE &aux_grid) {
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

unsigned long long Grid2DTEZ::return_maximum_time() {
    return max_time;
}

unsigned long long Grid2DTEZ::return_current_time() {
    return curr_time;
}

void Grid2DTEZ::save_state() {
    stringstream name;
    name<<"0-"<<"snapshot2D_TEZ_"<<curr_time<<".csv";
    ofstream file {name.str()};

    for (unsigned long long x = 0; x < size_x; ++x) {
        file<<'\n';
        for (unsigned long long y = 0; y < size_y; ++y) {
            file<<x<<','<<y<<','<<Hz(x, y)<<','<<Ex(x, y)<<','<<Ey(x, y)<<'\n';
        }
    }
    file.close();
    return;
}

void Grid2DTEZ::add_diel_rectangle(Point start, Point finish, double eps_r, double sigma, double mu_r, double sigma_m) {
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
                c1ey(i,j) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                c2ey(i,j) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                c1ex(i,j) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                c2ex(i,j) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                c1hz(i,j) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                c2hz(i,j) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
            }
        }

        return;
    } else {
        cout<<"Function add_diel_rectangle() failed, point out of bounds\nNo object added\n";
        return;
    }
}

void Grid2DTEZ::add_diel_circle(Point center, unsigned long long radius, double eps_r, double sigma, double mu_r,
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
                c1ey(i,j) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                c2ey(i,j) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                c1ex(i,j) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                c2ex(i,j) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                c1hz(i,j) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                c2hz(i,j) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                }
            }
        }

        return;

    } else {
        cout<<"Function add_diel_circle() failed, point out of bounds\nNo object added\n";
        return;
    }
}

void Grid2DTEZ::add_pec_circle(Point center, unsigned long long radius) {
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
                    c1ex(i,j) = 0.0;
                    c2ex(i,j) = 0.0;
                    c1ey(i,j) = 0.0;
                    c2ey(i,j) = 0.0;
                }
            }
        }
        return;

    } else {
        cout<<"Function add_pec_circle() failed, point out of bounds\nNo object added\n";
        return;
    }
}

void Grid2DTEZ::add_pec_rectangle(Point start, Point finish) {
    if (bounds_check(start) && bounds_check(finish)) {
        for (unsigned long long i = start.x; i < finish.x; ++i) {
            for (unsigned long long j = start.y; j < finish.y; ++j) {
                c1ex(i,j) = 0.0;
                c2ex(i,j) = 0.0;
                c1ey(i,j) = 0.0;
                c2ey(i,j) = 0.0;
            }
        }
        return;

    } else {
        cout<<"Function add_pec_rectangle() failed, point out of bounds\nNo object added\n";
        return;
    }

}


