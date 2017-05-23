/*--------------Grid 3D, #10--------------*/
#include "definitions.h"

bool Grid3D::bounds_check(Point a) {
    if (a.x < size_x && a.y < size_y && a.z < size_z) {return true;}
    return false;
}

void Grid3D::add_pec_rectangle(Point start, Point finish) {
    unsigned long long i, j, k;
    if (bounds_check(start) && bounds_check(finish)) {
        for (i = start.x; i < finish.x; ++i) {
            for (j = start.y;j < finish.y; ++j) {
                for (k = start.z; k < finish.z; ++k) {
                    //Ex nodes
                    c1ex(i,j,k) = 0.0;
                    c1ex(i,j+1,k) = 0.0;
                    c1ex(i,j,k+1) = 0.0;
                    c1ex(i,j+1,k+1) = 0.0;
                    c2ex(i,j,k) = 0.0;
                    c2ex(i,j+1,k) = 0.0;
                    c2ex(i,j,k+1) = 0.0;
                    c2ex(i,j+1,k+1) = 0.0;
                    //Ey nodes
                    c1ey(i,j,k) = 0.0;
                    c1ey(i+1,j,k) = 0.0;
                    c1ey(i,j,k+1) = 0.0;
                    c1ey(i+1,j,k+1) = 0.0;
                    c2ey(i,j,k) = 0.0;
                    c2ey(i+1,j,k) = 0.0;
                    c2ey(i,j,k+1) = 0.0;
                    c2ey(i+1,j,k+1) = 0.0;
                    //Ez nodes
                    c1ez(i,j,k) = 0.0;
                    c1ez(i+1,j,k) = 0.0;
                    c1ez(i,j+1,k) = 0.0;
                    c1ez(i+1,j+1,k) = 0.0;
                    c2ez(i,j,k) = 0.0;

                    c2ez(i+1,j,k) = 0.0;
                    c2ez(i,j+1,k) = 0.0;
                    c2ez(i+1,j+1,k) = 0.0;
                }
            }
        }

        return;
    } else {
        cout<<"Function add_pec_rectangle() failed, Point out of bounds\nNo object added\n";
        return;
    }
}

void Grid3D::add_pec_rectangle_centre(Point centre, unsigned long long xlen, unsigned long long ylen,
            unsigned long long zlen) {
unsigned long long i, j, k;
    xlen = xlen/2;
    ylen = ylen/2;
    zlen = zlen/2;

    Point start {centre.x - xlen, centre.y - ylen, centre.z - zlen},
    finish {centre.x + xlen, centre.y + ylen, centre.z + zlen};


    if (bounds_check(start) && bounds_check(finish)) {
        for (i = centre.x - xlen; i < centre.x + xlen; ++i) {
            for (j = centre.y - zlen;j < centre.y + ylen; ++j) {
                for (k = centre.z - zlen; k < centre.z + zlen; ++k) {
                    //Ex nodes
                    c1ex(i,j,k) = 0.0;
                    c1ex(i,j+1,k) = 0.0;
                    c1ex(i,j,k+1) = 0.0;
                    c1ex(i,j+1,k+1) = 0.0;
                    c2ex(i,j,k) = 0.0;
                    c2ex(i,j+1,k) = 0.0;
                    c2ex(i,j,k+1) = 0.0;
                    c2ex(i,j+1,k+1) = 0.0;
                    //Ey nodes
                    c1ey(i,j,k) = 0.0;
                    c1ey(i+1,j,k) = 0.0;
                    c1ey(i,j,k+1) = 0.0;
                    c1ey(i+1,j,k+1) = 0.0;
                    c2ey(i,j,k) = 0.0;
                    c2ey(i+1,j,k) = 0.0;
                    c2ey(i,j,k+1) = 0.0;
                    c2ey(i+1,j,k+1) = 0.0;
                    //Ez nodes
                    c1ez(i,j,k) = 0.0;
                    c1ez(i+1,j,k) = 0.0;
                    c1ez(i,j+1,k) = 0.0;
                    c1ez(i+1,j+1,k) = 0.0;
                    c2ez(i,j,k) = 0.0;
                    c2ez(i+1,j,k) = 0.0;
                    c2ez(i,j+1,k) = 0.0;
                    c2ez(i+1,j+1,k) = 0.0;
                }
            }
        }

        return;
    } else {
        cout<<"Function add_pec_rectangle_centre() failed, Point out of bounds\nNo object added\n";
        return;
    }
}

void Grid3D::add_diel_rectangle(Point start, Point finish, double eps_r, double sigma, double mu_r, double sigma_m) {
    unsigned long long i, j, k;
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

        for (i = start.x; i < finish.x; ++i) {
            for (j = start.y; j < finish.y; ++j) {
                for (k = start.z; k < finish.z; ++k) {
                    c1ez(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c2ez(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c1ex(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c2ex(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c1ey(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c2ey(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c1hx(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c2hx(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c1hy(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c2hy(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c1hz(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c2hz(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));

                }
            }
        }

        return;
    } else {
        cout<<"Function add_diel_rectangle() failed, point out of bounds\nNo object added\n";
        return;
    }
}

void Grid3D::add_diel_rectangle_centre(Point centre, unsigned long long xlen, unsigned long long ylen, unsigned long long zlen,
            double eps_r, double sigma, double mu_r, double sigma_m) {
    unsigned long long i, j, k;
    Point start {centre.x - xlen, centre.y - ylen, centre.z - zlen},
    finish {centre.x + xlen, centre.y + ylen, centre.z + zlen};
    xlen = xlen / 2;
    ylen = ylen / 2;
    zlen = zlen / 2;

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

        for (i = start.x; i < finish.x; ++i) {
            for (j = start.y; j < finish.y; ++j) {
                for (k = start.z; k < finish.z; ++k) {
                    c1ez(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c2ez(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c1ex(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c2ex(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c1ey(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c2ey(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                    c1hx(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c2hx(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c1hy(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c2hy(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c1hz(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    c2hz(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                }
            }
        }

        return;
    } else {
        cout<<"Function add_diel_rectangle_centre() failed, point out of bounds\nNo object added\n";
        return;
    }
}

void Grid3D::add_diel_sphere(Point centre, unsigned long long radius, double eps_r, double sigma, double mu_r,
                             double sigma_m) {
    unsigned long long i, j, k;
    Point start {centre.x - radius, centre.y - radius, centre.z - radius},
    finish {centre.x + radius, centre.y + radius, centre.z + radius};
    double x_sq, y_sq, z_sq, rad_sq, ii, jj, kk, centx, centy, centz;
    rad_sq = pow(radius, 2.0);
    centx = centre.x;
    centy = centre.y;
    centz = centre.z;

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

        for (i = start.x; i < finish.x; ++i) {
            ii = i;
            x_sq = pow((ii - centx), 2.0);
            for (j = start.y; j < finish.y; ++j) {
                jj = j;
                y_sq = pow((jj - centy), 2.0);
                for (k = start.z; k < finish.z; ++k) {
                    kk = k;
                    z_sq = pow((kk - centz), 2.0);

                    if (x_sq + y_sq + z_sq <= rad_sq ) {
                        c1ez(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c2ez(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c1ex(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c2ex(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c1ey(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c2ey(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c1hx(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c2hx(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c1hy(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c2hy(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c1hz(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c2hz(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    }
                }
            }
        }

        return;

    } else {
        cout<<"Function add_diel_sphere() failed, point out of bounds\nNo object added\n";
        return;
    }
}

void Grid3D::add_pec_sphere(Point centre, unsigned long long radius) {
    unsigned long long i, j, k;
    Point start {centre.x - radius, centre.y - radius, centre.z - radius},
    finish {centre.x + radius, centre.y + radius, centre.z + radius};
    double x_sq, y_sq, z_sq, rad_sq, ii, jj, kk, centx, centy, centz;
    rad_sq = pow(radius, 2.0);
    centx = centre.x;
    centy = centre.y;
    centz = centre.z;

    if (bounds_check(start) && bounds_check(finish)) {

        for (i = start.x; i < finish.x; ++i) {
            ii = i;
            x_sq = pow((ii - centx), 2.0);
            for (j = start.y; j < finish.y; ++j) {
                jj = j;
                y_sq = pow((jj - centy), 2.0);
                for (k = start.z; k < finish.z; ++k) {
                    kk = k;
                    z_sq = pow((kk - centz), 2.0);

                    if (x_sq + y_sq + z_sq <= rad_sq ) {
                        //Ex nodes
                        c1ex(i,j,k) = 0.0;
                        c1ex(i,j+1,k) = 0.0;
                        c1ex(i,j,k+1) = 0.0;
                        c1ex(i,j+1,k+1) = 0.0;
                        c2ex(i,j,k) = 0.0;
                        c2ex(i,j+1,k) = 0.0;
                        c2ex(i,j,k+1) = 0.0;
                        c2ex(i,j+1,k+1) = 0.0;
                        //Ey nodes
                        c1ey(i,j,k) = 0.0;
                        c1ey(i+1,j,k) = 0.0;
                        c1ey(i,j,k+1) = 0.0;
                        c1ey(i+1,j,k+1) = 0.0;
                        c2ey(i,j,k) = 0.0;
                        c2ey(i+1,j,k) = 0.0;
                        c2ey(i,j,k+1) = 0.0;
                        c2ey(i+1,j,k+1) = 0.0;
                        //Ez nodes
                        c1ez(i,j,k) = 0.0;
                        c1ez(i+1,j,k) = 0.0;
                        c1ez(i,j+1,k) = 0.0;
                        c1ez(i+1,j+1,k) = 0.0;
                        c2ez(i,j,k) = 0.0;
                        c2ez(i+1,j,k) = 0.0;
                        c2ez(i,j+1,k) = 0.0;
                        c2ez(i+1,j+1,k) = 0.0;
                    }
                }
            }
        }

        return;

    } else {
        cout<<"Function add_pec_sphere() failed, point out of bounds\nNo object added\n";
        return;
    }
}


void Grid3D::add_pec_cylinder(Point centre, unsigned long long radius, unsigned long long length, char opt) {
    switch (opt) {
    case 'x':
    {
    unsigned long long i, j, k;
    Point start {centre.x, centre.y - radius, centre.z - radius},
    finish {centre.x + length, centre.y + radius, centre.z + radius};
    double y_sq, z_sq, rad_sq, jj, kk, centy, centz;
    rad_sq = pow(radius, 2.0);
    centy = centre.y;
    centz = centre.z;

    if (bounds_check(start) && bounds_check(finish)) {

        for (i = start.x; i < finish.x; ++i) {
            for (j = start.y; j < finish.y; ++j) {
                jj = j;
                y_sq = pow((jj - centy), 2.0);
                for (k = start.z; k < finish.z; ++k) {
                    kk = k;
                    z_sq = pow((kk - centz), 2.0);

                    if (y_sq + z_sq <= rad_sq ) {
                        //Ex nodes
                        c1ex(i,j,k) = 0.0;
                        c1ex(i,j+1,k) = 0.0;
                        c1ex(i,j,k+1) = 0.0;
                        c1ex(i,j+1,k+1) = 0.0;
                        c2ex(i,j,k) = 0.0;
                        c2ex(i,j+1,k) = 0.0;
                        c2ex(i,j,k+1) = 0.0;
                        c2ex(i,j+1,k+1) = 0.0;
                        //Ey nodes
                        c1ey(i,j,k) = 0.0;
                        c1ey(i+1,j,k) = 0.0;
                        c1ey(i,j,k+1) = 0.0;
                        c1ey(i+1,j,k+1) = 0.0;
                        c2ey(i,j,k) = 0.0;
                        c2ey(i+1,j,k) = 0.0;
                        c2ey(i,j,k+1) = 0.0;
                        c2ey(i+1,j,k+1) = 0.0;
                        //Ez nodes
                        c1ez(i,j,k) = 0.0;
                        c1ez(i+1,j,k) = 0.0;
                        c1ez(i,j+1,k) = 0.0;
                        c1ez(i+1,j+1,k) = 0.0;
                        c2ez(i,j,k) = 0.0;
                        c2ez(i+1,j,k) = 0.0;
                        c2ez(i,j+1,k) = 0.0;
                        c2ez(i+1,j+1,k) = 0.0;
                    }
                }
            }
        }

        return;

    } else {
        cout<<"Function add_pec_cylinder() failed, point out of bounds\nNo object added\n";
        return;
    }
    } //End case 'x'

    case 'y':
    {
    unsigned long long i, j, k;
    Point start {centre.x - radius, centre.y, centre.z - radius},
    finish {centre.x + radius, centre.y + length, centre.z + radius};
    double x_sq, z_sq, rad_sq, ii, kk, centx, centz;
    rad_sq = pow(radius, 2.0);
    centx = centre.x;
    centz = centre.z;

    if (bounds_check(start) && bounds_check(finish)) {
        for (i = start.x; i < finish.x; ++i) {
            ii = i;
            x_sq = pow((ii - centx), 2.0);
            for (j = start.y; j < finish.y; ++j) {
                for (k = start.z; k < finish.z; ++k) {
                    kk = k;
                    z_sq = pow((kk - centz), 2.0);

                    if (x_sq + z_sq <= rad_sq ) {
                        //Ex nodes
                        c1ex(i,j,k) = 0.0;
                        c1ex(i,j+1,k) = 0.0;
                        c1ex(i,j,k+1) = 0.0;
                        c1ex(i,j+1,k+1) = 0.0;
                        c2ex(i,j,k) = 0.0;
                        c2ex(i,j+1,k) = 0.0;
                        c2ex(i,j,k+1) = 0.0;
                        c2ex(i,j+1,k+1) = 0.0;
                        //Ey nodes
                        c1ey(i,j,k) = 0.0;
                        c1ey(i+1,j,k) = 0.0;
                        c1ey(i,j,k+1) = 0.0;
                        c1ey(i+1,j,k+1) = 0.0;
                        c2ey(i,j,k) = 0.0;
                        c2ey(i+1,j,k) = 0.0;
                        c2ey(i,j,k+1) = 0.0;
                        c2ey(i+1,j,k+1) = 0.0;
                        //Ez nodes
                        c1ez(i,j,k) = 0.0;
                        c1ez(i+1,j,k) = 0.0;
                        c1ez(i,j+1,k) = 0.0;
                        c1ez(i+1,j+1,k) = 0.0;
                        c2ez(i,j,k) = 0.0;
                        c2ez(i+1,j,k) = 0.0;
                        c2ez(i,j+1,k) = 0.0;
                        c2ez(i+1,j+1,k) = 0.0;
                    }
                }
            }
        }

        return;

    } else {
        cout<<"Function add_pec_cylinder() failed, point out of bounds\nNo object added\n";
        return;
    }
    } //End case 'y'

        case 'z':
    {
    unsigned long long i, j, k;
    Point start {centre.x - radius, centre.y - radius, centre.z},
    finish {centre.x + radius, centre.y + radius, centre.z + length};
    double x_sq, y_sq, rad_sq, ii, jj, centx, centy;
    rad_sq = pow(radius, 2.0);
    centx = centre.x;
    centy = centre.y;

    if (bounds_check(start) && bounds_check(finish)) {
        for (i = start.x; i < finish.x; ++i) {
            ii = i;
            x_sq = pow((ii - centx), 2.0);
            for (j = start.y; j < finish.y; ++j) {
                jj = j;
                y_sq = pow((jj - centy), 2.0);
                for (k = start.z; k < finish.z; ++k) {

                    if (x_sq + y_sq <= rad_sq ) {
                        //Ex nodes
                        c1ex(i,j,k) = 0.0;
                        c1ex(i,j+1,k) = 0.0;
                        c1ex(i,j,k+1) = 0.0;
                        c1ex(i,j+1,k+1) = 0.0;
                        c2ex(i,j,k) = 0.0;
                        c2ex(i,j+1,k) = 0.0;
                        c2ex(i,j,k+1) = 0.0;
                        c2ex(i,j+1,k+1) = 0.0;
                        //Ey nodes
                        c1ey(i,j,k) = 0.0;
                        c1ey(i+1,j,k) = 0.0;
                        c1ey(i,j,k+1) = 0.0;
                        c1ey(i+1,j,k+1) = 0.0;
                        c2ey(i,j,k) = 0.0;
                        c2ey(i+1,j,k) = 0.0;
                        c2ey(i,j,k+1) = 0.0;
                        c2ey(i+1,j,k+1) = 0.0;
                        //Ez nodes
                        c1ez(i,j,k) = 0.0;
                        c1ez(i+1,j,k) = 0.0;
                        c1ez(i,j+1,k) = 0.0;
                        c1ez(i+1,j+1,k) = 0.0;
                        c2ez(i,j,k) = 0.0;
                        c2ez(i+1,j,k) = 0.0;
                        c2ez(i,j+1,k) = 0.0;
                        c2ez(i+1,j+1,k) = 0.0;
                    }
                }
            }
        }

        return;

    } else {
        cout<<"Function add_pec_cylinder() failed, point out of bounds\nNo object added\n";
        return;
    }
    } //End case 'z'

    default:
        {
            cout<<"Function add_pec_cylinder() failed, no such case\n";
            return;
        }
    }

}

void Grid3D::add_diel_cylinder(Point centre, unsigned long long radius, unsigned long long length, char opt,
double eps_r, double sigma, double mu_r, double sigma_m) {
    switch (opt) {
    case 'x':
    {
    unsigned long long i, j, k;
    double new_min_wl = wavelength / (sqrt(eps_r * mu_r));
    Point start {centre.x, centre.y - radius, centre.z - radius},
    finish {centre.x + length, centre.y + radius, centre.z + radius};
    double y_sq, z_sq, rad_sq, jj, kk, centy, centz;
    rad_sq = pow(radius, 2.0);
    centy = centre.y;
    centz = centre.z;

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

        for (i = start.x; i < finish.x; ++i) {
            for (j = start.y; j < finish.y; ++j) {
                jj = j;
                y_sq = pow((jj - centy), 2.0);
                for (k = start.z; k < finish.z; ++k) {
                    kk = k;
                    z_sq = pow((kk - centz), 2.0);

                    if (y_sq + z_sq <= rad_sq ) {
                        c1ez(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c2ez(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c1ex(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c2ex(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c1ey(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c2ey(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c1hx(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c2hx(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c1hy(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c2hy(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c1hz(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c2hz(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    }
                }
            }
        }

        return;

    } else {
        cout<<"Function add_diel_cylinder() failed, point out of bounds\nNo object added\n";
        return;
    }
    } //End case 'x'

    case 'y':
    {
    unsigned long long i, j, k;
    double new_min_wl = wavelength / (sqrt(eps_r * mu_r));
    Point start {centre.x - radius, centre.y, centre.z - radius},
    finish {centre.x + radius, centre.y + length, centre.z + radius};
    double x_sq, z_sq, rad_sq, ii, kk, centx, centz;
    rad_sq = pow(radius, 2.0);
    centx = centre.x;
    centz = centre.z;

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

        for (i = start.x; i < finish.x; ++i) {
            ii = i;
            x_sq = pow((ii - centx), 2.0);
            for (j = start.y; j < finish.y; ++j) {
                for (k = start.z; k < finish.z; ++k) {
                    kk = k;
                    z_sq = pow((kk - centz), 2.0);

                    if (x_sq + z_sq <= rad_sq ) {
                        c1ez(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c2ez(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c1ex(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c2ex(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c1ey(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c2ey(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c1hx(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c2hx(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c1hy(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c2hy(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c1hz(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c2hz(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    }
                }
            }
        }

        return;

    } else {
        cout<<"Function add_diel_cylinder() failed, point out of bounds\nNo object added\n";
        return;
    }
    } //End case 'y'

        case 'z':
    {
    unsigned long long i, j, k;
    double new_min_wl = wavelength / (sqrt(eps_r * mu_r));
    Point start {centre.x - radius, centre.y - radius, centre.z},
    finish {centre.x + radius, centre.y + radius, centre.z + length};
    double x_sq, y_sq, rad_sq, ii, jj, centx, centy;
    rad_sq = pow(radius, 2.0);
    centx = centre.x;
    centy = centre.y;

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

        for (i = start.x; i < finish.x; ++i) {
            ii = i;
            x_sq = pow((ii - centx), 2.0);
            for (j = start.y; j < finish.y; ++j) {
                jj = j;
                y_sq = pow((jj - centy), 2.0);
                for (k = start.z; k < finish.z; ++k) {

                    if (x_sq + y_sq <= rad_sq ) {
                        c1ez(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c2ez(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c1ex(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c2ex(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c1ey(i,j,k) = (1.0 - sigma * dt / (2.0 * eps_r * eps_0))/(1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c2ey(i,j,k) = (imp0 * courant / eps_r) / (1.0 + sigma * dt / (2.0 * eps_r * eps_0));
                        c1hx(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c2hx(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c1hy(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c2hy(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c1hz(i,j,k) = (1.0 - sigma_m * dt / (2.0 * mu_r * mu_0))/(1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                        c2hz(i,j,k) = (courant / (mu_r * imp0)) / (1.0 + sigma * dt / (2.0 * mu_r * mu_0));
                    }
                }
            }
        }

        return;

    } else {
        cout<<"Function add_diel_cylinder() failed, point out of bounds\nNo object added\n";
        return;
    }
    } //End case 'z'

    default:
        {
            cout<<"Function add_diel_cylinder() failed, no such case\n";
            return;
        }
    }

}
