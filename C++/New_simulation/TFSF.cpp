/*--------------Grid 3D, #7.b--------------*/
#include "definitions.h"

//Rectangular TFSF
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
//=======================================================================
//TFSF as a cylinder on the x direction
Point Grid3D::find_TFSF_centre() {
    //Find the centre of TFSF area, centered around the starting face of the TFSF area
    unsigned long long deltay, deltaz;
    Point st {tfsf_start.x, tfsf_start.y, tfsf_start.z},
            fin {tfsf_end.x, tfsf_end.y, tfsf_end.z};
    deltay = fin.y - st.y;
    deltaz = fin.z - st.z;
    Point cent {st.x, st.y + deltay/2, st.z + deltaz/2};
    return cent;
}

double Grid3D::find_TFSF_radius() {
    //This function finds the radius of the TFSF area, defined as the distance between the centre of the Square area
    //to a corner vertex
    Point centre =find_TFSF_centre();
    double radius = sqrt(pow(tfsf_end.y - centre.y, 2.0) + pow(tfsf_end.z - centre.z, 2.0));

    return radius;
}

/*
void Grid3D::correct_Hx_0_cilynder(Grid1DTM &aux_grid) {
    unsigned long long i,j, k;
    double jj, kk, radius;
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
*/
