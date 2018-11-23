/*--------------First order Absorbing Boundary Condition 3D, #9--------------*/
#include "definitions.h"

//First order ABC
void Grid3D::simple_x0() {
    unsigned long long i, j, k;
    i = 0;
    double coef = (courant - 1) / (courant + 1);

    for (j = 0; j < size_y - 1; ++j) {
        for (k = 0; k < size_z; ++k) {
            Ey(i,j,k) = Ey_x0(j,k,0) + coef * (Ey(i+1,j,k) - Ey(i,j,k));
            Ey_x0(j,k,0) = Ey(i+1,j,k);
        }
    }

    for (j = 0; j < size_y; ++j) {
        for (k = 0; k < size_z - 1; ++k) {
            Ez(i,j,k) = Ez_x0(j,k,0) + coef * (Ez(i+1,j,k) - Ez(i,j,k));
            Ez_x0(j,k,0) = Ez(i+1,j,k);
        }
    }

    return;
}

void Grid3D::simple_xf() {
    unsigned long long i, j, k;
    i = size_x - 1;
    double coef = (courant - 1) / (courant + 1);

    for (j = 0; j < size_y - 1; ++j) {
        for (k = 0; k < size_z; ++k) {
            Ey(i,j,k) = Ey_xf(j,k,0) + coef * (Ey(i-1,j,k) - Ey(i,j,k));
            Ey_xf(j,k,0) = Ey(i-1,j,k);
        }
    }

    for (j = 0; j < size_y; ++j) {
        for (k = 0; k < size_z - 1; ++k) {
            Ez(i,j,k) = Ez_xf(j,k,0) + coef * (Ez(i-1,j,k) - Ez(i,j,k));
            Ez_xf(j,k,0) = Ez(i-1,j,k);
        }
    }

    return;
}

void Grid3D::simple_y0() {
    unsigned long long i, j, k;
    j = 0;
    double coef = (courant - 1) / (courant + 1);

    for (i = 0; i < size_x - 1; ++i) {
        for (k = 0; k < size_z; ++k) {
            Ex(i,j,k) = Ex_y0(i,k,0) + coef * (Ex(i,j+1,k) - Ex(i,j,k));
            Ex_y0(i,k,0) = Ex(i,j+1,k);
        }
    }

    for (i = 0; i < size_x; ++i) {
        for (k = 0; k < size_z - 1; ++k) {
            Ez(i,j,k) = Ez_y0(i,k,0) + coef * (Ez(i,j+1,k) - Ez(i,j,k));
            Ez_y0(i,k,0) = Ez(i,j+1,k);
        }
    }

    return;
}

void Grid3D::simple_yf() {
    unsigned long long i, j, k;
    j = size_y - 1;
    double coef = (courant - 1) / (courant + 1);

    for (i = 0; i < size_x - 1; ++i) {
        for (k = 0; k < size_z; ++k) {
            Ex(i,j,k) = Ex_yf(i,k,0) + coef * (Ex(i,j-1,k) - Ex(i,j,k));
            Ex_yf(i,k,0) = Ex(i,j-1,k);
        }
    }

    for (i = 0; i < size_x; ++i) {
        for (k = 0; k < size_z - 1; ++k) {
            Ez(i,j,k) = Ez_yf(i,k,0) + coef * (Ez(i,j-1,k) - Ez(i,j,k));
            Ez_yf(i,k,0) = Ez(i,j-1,k);
        }
    }

    return;
}

void Grid3D::simple_z0() {
    unsigned long long i, j, k;
    k = 0;
    double coef = (courant - 1) / (courant + 1);

    for (i = 0; i < size_x - 1; ++i) {
        for (j = 0; j < size_y; ++j) {
            Ex(i,j,k) = Ex_z0(i,j,0) + coef * (Ex(i,j,k+1) - Ex(i,j,k));
            Ex_z0(i,j,0) = Ex(i,j,k+1);
        }
    }

    for (i = 0; i < size_x; ++i) {
        for (j = 0; j < size_y - 1; ++j) {
            Ey(i,j,k) = Ey_z0(i,j,0) + coef * (Ey(i,j,k+1) - Ey(i,j,k));
            Ey_z0(i,j,0) = Ey(i,j,k+1);
        }
    }

    return;
}


void Grid3D::simple_zf() {
    unsigned long long i, j, k;
    k = size_z - 1;
    double coef = (courant - 1) / (courant + 1);

    for (i = 0; i < size_x - 1; ++i) {
        for (j = 0; j < size_y; ++j) {
            Ex(i,j,k) = Ex_zf(i,j,0) + coef * (Ex(i,j,k-1) - Ex(i,j,k));
            Ex_zf(i,j,0) = Ex(i,j,k-1);
        }
    }

    for (i = 0; i < size_x; ++i) {
        for (j = 0; j < size_y - 1; ++j) {
            Ey(i,j,k) = Ey_zf(i,j,0) + coef * (Ey(i,j,k-1) - Ey(i,j,k));
            Ey_zf(i,j,0) = Ey(i,j,k-1);
        }
    }

    return;
}

