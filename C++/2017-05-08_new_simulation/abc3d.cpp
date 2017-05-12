/*--------------Second order Absorbing Boundary Condition 3D, #8--------------*/
#include "definitions.h"


/////y face
void Grid3D::abc_Ex_y0() {
    unsigned long long i, j, k;
    j = 0; //y0 face

    for(i = 0; i < size_x; ++i) {
        for (k = 0; k < size_z; ++k) {
            cour_prime = c2ex(i,j,k) * c2hy(i,j,k);
            A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
            B = (1.0/cour_prime) - 2 + cour_prime;
            C = 2.0 * (cour_prime - 1/cour_prime);
            D = -4.0 * (cour_prime + 1/cour_prime);
            Ex(i, j, k) = A * (B * (Ex(i, 2, k) + Ex_y0(i, k, 0)) + C * (Ex_y0(i, k, 3) + Ex_y0(i, k, 5) - Ex(i, 1, k)
                - Ex_y0(i, k, 1)) + D * Ex_y0(i, k, 4)) - Ex_y0(i, k, 2);


            //Remember old fields, t-1 becomes t and t becomes t+1
            for (j = 0; j < 3; ++j) {
                Ex_y0(i, k, j) = Ex_y0(i, k, j+3);
                Ex_y0(i, k, j+3) = Ex(i, j, k);
            }
        }
    }

    return;
}

void Grid3D::abc_Ez_y0() {
    unsigned long long i, j, k;
    j = 0; //y0 face

    for(i = 0; i < size_x; ++i) {
        for (k = 0; k < size_z; ++k) {
            cour_prime = c2ez(i,j,k) * c2hx(i,j,k);
            A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
            B = (1.0/cour_prime) - 2 + cour_prime;
            C = 2.0 * (cour_prime - 1/cour_prime);
            D = -4.0 * (cour_prime + 1/cour_prime);
            Ez(i, j, k) = A * (B * (Ez(i, 2, k) + Ez_y0(i, k, 0)) + C * (Ez_y0(i, k, 3) + Ez_y0(i, k, 5) - Ez(i, 1, k)
                - Ez_y0(i, k, 1)) + D * Ez_y0(i, k, 4)) - Ez_y0(i, k, 2);


            //Remember old fields, t-1 becomes t and t becomes t+1
            for (j = 0; j < 3; ++j) {
                Ez_y0(i, k, j) = Ez_y0(i, k, j+3);
                Ez_y0(i, k, j+3) = Ez(i, j, k);
            }
        }
    }

    return;
}

void Grid3D::abc_Ex_yf() {
    unsigned long long i, j, k;
    j = size_y - 1; //yf face

    for(i = 0; i < size_x; ++i) {
        for (k = 0; k < size_z; ++k) {
            cour_prime = c2ex(i,j,k) * c2hy(i,j,k);
            A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
            B = (1.0/cour_prime) - 2 + cour_prime;
            C = 2.0 * (cour_prime - 1/cour_prime);
            D = -4.0 * (cour_prime + 1/cour_prime);
            Ex(i, j, k) = A * (B * (Ex(i, j-2, k) + Ex_yf(i, k, 0)) + C * (Ex_yf(i, k, 3) + Ex_yf(i, k, 5) - Ex(i, j-1, k)
                - Ex_yf(i, k, 1)) + D * Ex_yf(i, k, 4)) - Ex_yf(i, k, 2);


            //Remember old fields, t-1 becomes t and t becomes t+1
            for (j = 0; j < 3; ++j) {
                Ex_yf(i, k, j) = Ex_yf(i, k, j+3);
                Ex_yf(i, k, j+3) = Ex(i, size_y-1-j, k);
            }
        }
    }

    return;
}

void Grid3D::abc_Ez_yf() {
    unsigned long long i, j, k;
    j = size_y - 1; //yf face

    for(i = 0; i < size_x; ++i) {
        for (k = 0; k < size_z; ++k) {
            cour_prime = c2ez(i,j,k) * c2hx(i,j,k);
            A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
            B = (1.0/cour_prime) - 2 + cour_prime;
            C = 2.0 * (cour_prime - 1/cour_prime);
            D = -4.0 * (cour_prime + 1/cour_prime);
            Ez(i, j, k) = A * (B * (Ez(i, j-2, k) + Ez_yf(i, k, 0)) + C * (Ez_yf(i, k, 3) + Ez_yf(i, k, 5) - Ez(i, j-1, k)
                - Ez_yf(i, k, 1)) + D * Ez_yf(i, k, 4)) - Ez_yf(i, k, 2);


            //Remember old fields, t-1 becomes t and t becomes t+1
            for (j = 0; j < 3; ++j) {
                Ez_yf(i, k, j) = Ez_yf(i, k, j+3);
                Ez_yf(i, k, j+3) = Ez(i, size_y-1-j, k);
            }
        }
    }

    return;
}
//////////x face
void Grid3D::abc_Ey_x0() {
    unsigned long long i, j, k;
    i = 0; //x0 face

    for(j = 0; j < size_y; ++j) {
        for (k = 0; k < size_z; ++k) {
            cour_prime = c2ey(i,j,k) * c2hz(i,j,k);
            A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
            B = (1.0/cour_prime) - 2 + cour_prime;
            C = 2.0 * (cour_prime - 1/cour_prime);
            D = -4.0 * (cour_prime + 1/cour_prime);
            Ey(i, j, k) = A * (B * (Ey(2, j, k) + Ey_x0(j, k, 0)) + C * (Ey_x0(j, k, 3) + Ey_x0(j, k, 5) - Ey(1, j, k)
                - Ey_x0(j, k, 1)) + D * Ey_x0(j, k, 4)) - Ey_x0(j, k, 2);


            //Remember old fields, t-1 becomes t and t becomes t+1
            for (i = 0; i < 3; ++i) {
                Ey_x0(j, k, i) = Ey_x0(j, k, i+3);
                Ey_x0(j, k, i+3) = Ey(i, j, k);
            }
        }
    }

    return;
}

void Grid3D::abc_Ez_x0() {
    unsigned long long i, j, k;
    i = 0; //x0 face

    for(j = 0; j < size_y; ++j) {
        for (k = 0; k < size_z; ++k) {
            cour_prime = c2ez(i,j,k) * c2hx(i,j,k);
            A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
            B = (1.0/cour_prime) - 2 + cour_prime;
            C = 2.0 * (cour_prime - 1/cour_prime);
            D = -4.0 * (cour_prime + 1/cour_prime);
            Ez(i, j, k) = A * (B * (Ez(2, j, k) + Ez_x0(j, k, 0)) + C * (Ez_x0(j, k, 3) + Ez_x0(j, k, 5) - Ez(1, j, k)
                - Ez_x0(j, k, 1)) + D * Ez_x0(j, k, 4)) - Ez_x0(j, k, 2);


            //Remember old fields, t-1 becomes t and t becomes t+1
            for (i = 0; i < 3; ++i) {
                Ez_x0(j, k, i) = Ez_x0(j, k, i+3);
                Ez_x0(j, k, i+3) = Ez(i, j, k);
            }
        }
    }

    return;
}

void Grid3D::abc_Ey_xf() {
    unsigned long long i, j, k;
    i = size_x - 1; //xf face

    for(j = 0; j < size_y; ++j) {
        for (k = 0; k < size_z; ++k) {
            cour_prime = c2ey(i,j,k) * c2hz(i,j,k);
            A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
            B = (1.0/cour_prime) - 2 + cour_prime;
            C = 2.0 * (cour_prime - 1/cour_prime);
            D = -4.0 * (cour_prime + 1/cour_prime);
            Ey(i, j, k) = A * (B * (Ey(i-2, j, k) + Ey_xf(j, k, 0)) + C * (Ey_xf(j, k, 3) + Ey_xf(j, k, 5) - Ey(i-1, j, k)
                - Ey_xf(j, k, 1)) + D * Ey_xf(j, k, 4)) - Ey_xf(j, k, 2);


            //Remember old fields, t-1 becomes t and t becomes t+1
            for (i = 0; i < 3; ++i) {
                Ey_xf(j, k, i) = Ey_xf(j, k, i+3);
                Ey_xf(j, k, i+3) = Ey(size_x-1-i, j, k);
            }
        }
    }

    return;
}

void Grid3D::abc_Ez_xf() {
    unsigned long long i, j, k;
    i = size_x - 1; //x0 face

    for(j = 0; j < size_y; ++j) {
        for (k = 0; k < size_z; ++k) {
            cour_prime = c2ez(i,j,k) * c2hx(i,j,k);
            A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
            B = (1.0/cour_prime) - 2 + cour_prime;
            C = 2.0 * (cour_prime - 1/cour_prime);
            D = -4.0 * (cour_prime + 1/cour_prime);
            Ez(i, j, k) = A * (B * (Ez(i-2, j, k) + Ez_xf(j, k, 0)) + C * (Ez_xf(j, k, 3) + Ez_xf(j, k, 5) - Ez(i-1, j, k)
                - Ez_xf(j, k, 1)) + D * Ez_xf(j, k, 4)) - Ez_xf(j, k, 2);


            //Remember old fields, t-1 becomes t and t becomes t+1
            for (i = 0; i < 3; ++i) {
                Ez_xf(j, k, i) = Ez_xf(j, k, i+3);
                Ez_xf(j, k, i+3) = Ez(size_x-1-i, j, k);
            }
        }
    }

    return;
}

////z face
void Grid3D::abc_Ex_z0() {
    unsigned long long i, j, k;
    k = 0; //z0 face

    for(i = 0; i < size_x; ++i) {
        for (j = 0; j < size_y; ++j) {
            cour_prime = c2ex(i,j,k) * c2hy(i,j,k);
            A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
            B = (1.0/cour_prime) - 2 + cour_prime;
            C = 2.0 * (cour_prime - 1/cour_prime);
            D = -4.0 * (cour_prime + 1/cour_prime);
            Ex(i, j, k) = A * (B * (Ex(i, j, 2) + Ex_z0(i, j, 0)) + C * (Ex_z0(i, j, 3) + Ex_z0(i, j, 5) - Ex(i, j, 1)
                - Ex_z0(i, j, 1)) + D * Ex_z0(i, j, 4)) - Ex_z0(i, j, 2);


            //Remember old fields, t-1 becomes t and t becomes t+1
            for (k = 0; k < 3; ++k) {
                Ex_z0(i, j, k) = Ex_z0(i, j, k+3);
                Ex_z0(i, j, k+3) = Ex(i, j, k);
            }
        }
    }

    return;
}

void Grid3D::abc_Ey_z0() {
    unsigned long long i, j, k;
    k = 0; //z0 face

    for(i = 0; i < size_x; ++i) {
        for (j = 0; j < size_y; ++j) {
            cour_prime = c2ey(i,j,k) * c2hz(i,j,k);
            A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
            B = (1.0/cour_prime) - 2 + cour_prime;
            C = 2.0 * (cour_prime - 1/cour_prime);
            D = -4.0 * (cour_prime + 1/cour_prime);
            Ey(i, j, k) = A * (B * (Ey(i, j, 2) + Ey_z0(i, j, 0)) + C * (Ey_z0(i, j, 3) + Ey_z0(i, j, 5) - Ey(i, j, 1)
                - Ey_z0(i, j, 1)) + D * Ey_z0(i, j, 4)) - Ey_z0(i, j, 2);


            //Remember old fields, t-1 becomes t and t becomes t+1
            for (k = 0; k < 3; ++k) {
                Ey_z0(i, j, k) = Ey_z0(i, j, k+3);
                Ey_z0(i, j, k+3) = Ey(i, j, k);
            }
        }
    }

    return;
}

void Grid3D::abc_Ex_zf() {
    unsigned long long i, j, k;
    k = size_z - 1; //zf face

    for(i = 0; i < size_x; ++i) {
        for (j = 0; j < size_y; ++j) {
            cour_prime = c2ex(i,j,k) * c2hy(i,j,k);
            A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
            B = (1.0/cour_prime) - 2 + cour_prime;
            C = 2.0 * (cour_prime - 1/cour_prime);
            D = -4.0 * (cour_prime + 1/cour_prime);
            Ex(i, j, k) = A * (B * (Ex(i, j, k-2) + Ex_zf(i, j, 0)) + C * (Ex_zf(i, j, 3) + Ex_zf(i, j, 5) - Ex(i, j, k-1)
                - Ex_zf(i, j, 1)) + D * Ex_zf(i, j, 4)) - Ex_zf(i, j, 2);


            //Remember old fields, t-1 becomes t and t becomes t+1
            for (k = 0; k < 3; ++k) {
                Ex_zf(i, j, k) = Ex_zf(i, j, k+3);
                Ex_zf(i, j, k+3) = Ex(i, j, size_z-1-k);
            }
        }
    }

    return;
}

void Grid3D::abc_Ey_zf() {
    unsigned long long i, j, k;
    k = size_z - 1; //z0 face

    for(i = 0; i < size_x; ++i) {
        for (j = 0; j < size_y; ++j) {
            cour_prime = c2ey(i,j,k) * c2hz(i,j,k);
            A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
            B = (1.0/cour_prime) - 2 + cour_prime;
            C = 2.0 * (cour_prime - 1/cour_prime);
            D = -4.0 * (cour_prime + 1/cour_prime);
            Ey(i, j, k) = A * (B * (Ey(i, j, k-2) + Ey_zf(i, j, 0)) + C * (Ey_zf(i, j, 3) + Ey_zf(i, j, 5) - Ey(i, j, k-1)
                - Ey_zf(i, j, 1)) + D * Ey_zf(i, j, 4)) - Ey_zf(i, j, 2);


            //Remember old fields, t-1 becomes t and t becomes t+1
            for (k = 0; k < 3; ++k) {
                Ey_zf(i, j, k) = Ey_zf(i, j, k+3);
                Ey_zf(i, j, k+3) = Ey(i, j, size_z-1-k);
            }
        }
    }

    return;
}
