#include "headers.h"

using namespace std;

//Dielectric_sq functions
Dielectric_sq::Dielectric_sq()
: length{1}, width{1}, height{1}, mu_rel{1.0}, eps_rel{1.0} {
    cout<<"Warning Dielectric_sq created without dimensions, please change dimensions\n";
}

Dielectric_sq::Dielectric_sq(unsigned long l, unsigned long w, unsigned long h, double mu, double eps)
: length{l}, width{w}, height{h}, mu_rel{mu}, eps_rel{eps} {
    cout<<"Dielectric_sq created\n";
}

//Grid functions
Grid::Grid() : obj{def_obj}, I{5}, J{5}, K{5} {
    cout<<"Warning no Dielectric / grid values specified, using default\n";
}

Grid::Grid(Dielectric_sq &diel, double wll, unsigned long l, unsigned long w, unsigned long h)
: obj{diel}, wl{wll}, I{l}, J{w}, K{h} {
    this->vsize = l * w * h;
    //Place the dielectric object in the middle of the grid
    unsigned long long ll = diel.length;
    unsigned long long hh = diel.height;
    unsigned long long ww = diel.width;
    this->diel_x_start = (l/2) - (ll/2);
    this->diel_x_end = (l/2) + (ll/2);
    this->diel_y_start = (w/2) - (ww/2);
    this->diel_y_end = (w/2) + (ww/2);
    this->diel_z_start = (h/2) - (hh/2);
    this->diel_z_end = (h/2) + (hh/2);

    double wavelength = wll; //1550 nm
    double mu_rel = diel.mu_rel; //Relative permeability
    double epsilon_rel = diel.eps_rel; //Relative permittivity
    this->n_core = sqrt(mu_rel * epsilon_rel); //Refractive index of the core
    double frequency = c0 / wavelength;
    this->frequency = frequency;
    double adjusted_wl = (c0 / this->n_core) / frequency;
    this->ds = adjusted_wl / 10.0; //10 steps per wavelength
    this->dt = this->ds / (2.0 * c0);
    this->courant = c0 * this->dt / this->ds; //min courant for convergence is ~0.577 1/sqrt(3)
    this->abs_fact = ((this->courant)-1.0) / ((this->courant)+1.0);
    init_properties(); //Initializes sigma, sigma ast, mu, epsilon, C_a, C_b, D_a, D_b vectors
    init_electric(); //Initializes E_x, E_y, E_z vectors
    init_magnetic(); //Initializes H_x, H_y, H_z vectors
    init_boundary(); //Initializes E_y_x0, E_y_x1, E_z_x0, E_z_x1, E_x_y0, E_x_y1, E_z_y0, E_z_y1,
                          //E_x_z0, E_x_z1, E_y_z0, E_y_z1;


    cout<<"Grid object created\n";
}

//Function to get index from flattened 3D vetors
unsigned long long Grid::get_ind(const unsigned long long i, const unsigned long long j, const unsigned long long k) {
    unsigned long long I = this->I;
    unsigned long long J = this->J;
    unsigned long long K = this->K;
    unsigned long long ind = k + K * (j + J * i);
    if (ind > (I * J * K - 1) || ind < 0) {
        cout<<"Warning, index out of range\n";
    }
    return ind;
}

//Function to get index from flattened 2D vectors
unsigned long long Grid::get_ind_2d(const unsigned long long a, const unsigned long long b, const unsigned long long L) {
    unsigned long long ind = b + L * (a);
    return ind;
}

void Grid::init_properties() {
    for (unsigned long long i = 0; i<this->I; ++i) {
        for (unsigned long long j = 0; j<this->J; ++j) {
            for (unsigned long long k = 0; k<this->K; ++k) {
                if (diel_x_start< i < this->diel_x_end && diel_y_start< i < this->diel_y_end
                     && diel_z_start< i < this->diel_z_end) {
                    this->mu.push_back(mu_0 * this->obj.mu_rel);
                    this->epsilon.push_back(epsilon_0 * this->obj.eps_rel);
                } else {
                    this->mu.push_back(mu_0 * 1.0);
                    this->epsilon.push_back(epsilon_0 * 1.0);
                }

                this->sigma.push_back(0.0);
                this->sigma_ast.push_back(0.0);
                this->C_a.push_back((1.0 - this->sigma[get_ind(i,j,k)] * dt/ (2.0 * this->epsilon[get_ind(i,j,k)])) /
                                    (1.0 + this->sigma[get_ind(i,j,k)] * dt/ (2.0 * this->epsilon[get_ind(i,j,k)])));

                this->C_b.push_back((this->dt / (this->epsilon[get_ind(i,j,k)])) /
                                    (1.0 + this->sigma[get_ind(i,j,k)] * this->dt/ (2.0 * this->epsilon[get_ind(i,j,k)])));

                this->D_a.push_back((1.0 - this->sigma_ast[get_ind(i,j,k)] * this->dt/ (2.0 * this->mu[get_ind(i,j,k)])) /
                                    (1.0 + this->sigma_ast[get_ind(i,j,k)] * this->dt/ (2.0 * this->mu[get_ind(i,j,k)])));

                this->D_b.push_back((this->dt / (this->mu[get_ind(i,j,k)])) /
                                    (1.0 + this->sigma_ast[get_ind(i,j,k)] * this->dt/ (2.0 * this->mu[get_ind(i,j,k)])));

            }
        }
    }
}

void Grid::init_electric() {
    for (unsigned long long i = 0; i<this->I; ++i) {
        for (unsigned long long j = 0; j<this->J; ++j) {
            for (unsigned long long k = 0; k<this->K; ++k) {
                this->E_x.push_back(0.0);
                this->E_y.push_back(0.0);
                this->E_z.push_back(0.0);

            }
        }
    }
}

void Grid::init_magnetic() {
    for (unsigned long long i = 0; i<this->I; ++i) {
        for (unsigned long long j = 0; j<this->J; ++j) {
            for (unsigned long long k = 0; k<this->K; ++k) {
                this->H_x.push_back(0.0);
                this->H_y.push_back(0.0);
                this->H_z.push_back(0.0);

            }
        }
    }
}


void Grid::init_boundary() {
    for (unsigned long long i = 0; i<this->I; ++i) {
        for (unsigned long long j = 0; j<this->J; ++j) {
            E_x_z0.push_back(0.0);
            E_x_z1.push_back(0.0);
            E_y_z0.push_back(0.0);
            E_y_z1.push_back(0.0);

        }
    }

    for (unsigned long long i = 0; i<this->I; ++i) {
        for (unsigned long long k = 0; k<this->K; ++k) {
            E_x_y0.push_back(0.0);
            E_x_y1.push_back(0.0);
            E_z_y0.push_back(0.0);
            E_z_y1.push_back(0.0);
        }
    }

    for (unsigned long long j = 0; j<this->J; ++j) {
        for (unsigned long long k = 0; k<this->K; ++k) {
            E_y_x0.push_back(0.0);
            E_y_x1.push_back(0.0);
            E_z_x0.push_back(0.0);
            E_z_x1.push_back(0.0);
        }
    }
}


void Grid::update_boundary_x(unsigned long J, unsigned long K, double face, vector<double> &E_y_x0, vector<double> &E_y,
                           double abs_fact) {
    if (face == 0) {
        for (unsigned long j = 0; j<J; ++j) {
            for (unsigned long k = 0; k<K; ++k) {
                E_y[get_ind(face,j,k)] = E_y_x0[get_ind_2d(j,k,K)] +
                abs_fact * (E_y[get_ind(face + 1, j, k)] - E_y[get_ind(face, j, k)]);

                E_y_x0[get_ind_2d(j,k,K)] = E_y[get_ind(face + 1, j, k)];
            }
        }
    } else {
        for (unsigned long j = 0; j<J; ++j) {
            for (unsigned long k = 0; k<K; ++k) {
                E_y[get_ind(face,j,k)] = E_y_x0[get_ind_2d(j,k,K)] +
                abs_fact * (E_y[get_ind(face - 1, j, k)] - E_y[get_ind(face, j, k)]);

                E_y_x0[get_ind_2d(j,k,K)] = E_y[get_ind(face - 1, j, k)];
            }
        }
    }
}

void Grid::update_boundary_y(unsigned long J, unsigned long K, double face, vector<double> &E_y_x0, vector<double> &E_y,
                           double abs_fact) {
    if (face == 0) {
        for (unsigned long j = 0; j<J; ++j) {
            for (unsigned long k = 0; k<K; ++k) {
                E_y[get_ind(j,face,k)] = E_y_x0[get_ind_2d(j,k,K)] +
                abs_fact * (E_y[get_ind(j, face + 1, k)] - E_y[get_ind(j, face, k)]);

                E_y_x0[get_ind_2d(j,k,K)] = E_y[get_ind(j, face + 1, k)];
            }
        }
    } else {
        for (unsigned long j = 0; j<J; ++j) {
            for (unsigned long k = 0; k<K; ++k) {
                E_y[get_ind(j, face, k)] = E_y_x0[get_ind_2d(j,k,K)] +
                abs_fact * (E_y[get_ind(j, face - 1, k)] - E_y[get_ind(j, face, k)]);

                E_y_x0[get_ind_2d(j,k,K)] = E_y[get_ind(j, face - 1, k)];
            }
        }
    }
}

void Grid::update_boundary_z(unsigned long J, unsigned long K, double face, vector<double> &E_y_x0, vector<double> &E_y,
                           double abs_fact) {
    if (face == 0) {
        for (unsigned long j = 0; j<J; ++j) {
            for (unsigned long k = 0; k<K; ++k) {
                E_y[get_ind(j,k, face)] = E_y_x0[get_ind_2d(j,k,K)] +
                abs_fact * (E_y[get_ind(j, k, face + 1)] - E_y[get_ind(j, k, face)]);

                E_y_x0[get_ind_2d(j,k,K)] = E_y[get_ind(j, k, face + 1)];
            }
        }
    } else {
        for (unsigned long j = 0; j<J; ++j) {
            for (unsigned long k = 0; k<K; ++k) {
                E_y[get_ind(j,k, face)] = E_y_x0[get_ind_2d(j,k,K)] +
                abs_fact * (E_y[get_ind(j, k, face - 1)] - E_y[get_ind(j, k, face)]);

                E_y_x0[get_ind_2d(j,k,K)] = E_y[get_ind(j, k, face - 1)];
            }
        }
    }
}

void Grid::parallel_update_boundary_all() {
    vector<thread> threads;
    //x face
    threads.push_back(thread(&Grid::update_boundary_x, this, this->J-1, this->K, 0, ref(this->E_y_x0), ref(this->E_y), this->abs_fact));
    threads.push_back(thread(&Grid::update_boundary_x, this, this->J, this->K-1, 0, ref(this->E_z_x0), ref(this->E_z), this->abs_fact));
    threads.push_back(thread(&Grid::update_boundary_x, this, this->J-1, this->K, this->I-1, ref(this->E_y_x1), ref(this->E_y),
                             this->abs_fact));
    threads.push_back(thread(&Grid::update_boundary_x, this, this->J, this->K-1, this->I-1, ref(this->E_z_x1), ref(this->E_z),
                             this->abs_fact));
    for (thread &th : threads) {
        th.join();
    }
    threads.clear();

    //y face
    threads.push_back(thread(&Grid::update_boundary_y, this, this->I-1, this->K, 0, ref(this->E_x_y0), ref(this->E_x),
                             this->abs_fact));
    threads.push_back(thread(&Grid::update_boundary_y, this, this->I, this->K-1, 0, ref(this->E_z_y0), ref(this->E_z),
                      this->abs_fact));
    threads.push_back(thread(&Grid::update_boundary_y, this, this->I-1, this->K, this->J-1, ref(this->E_x_y1), ref(this->E_x),
                      this->abs_fact));
    threads.push_back(thread(&Grid::update_boundary_y, this, this->I, this->K-1, this->J-1, ref(this->E_z_y1), ref(this->E_z),
                      this->abs_fact));
    for (thread &th : threads) {
        th.join();
    }
    threads.clear();

    //z face
    threads.push_back(thread(&Grid::update_boundary_z, this, this->I-1, this->J, 0, ref(this->E_x_z0), ref(this->E_x),
                      this->abs_fact));
    threads.push_back(thread(&Grid::update_boundary_z, this, this->I, this->J-1, 0, ref(this->E_y_z0), ref(this->E_y),
                      this->abs_fact));
    threads.push_back(thread(&Grid::update_boundary_z, this, this->I-1, this->J, this->K-1, ref(this->E_x_z1), ref(this->E_x),
                      this->abs_fact));
    threads.push_back(thread(&Grid::update_boundary_z, this, this->I, this->J-1, this->K-1, ref(this->E_y_z1), ref(this->E_y),
                      this->abs_fact));
    for (thread &th : threads) {
        th.join();
    }
}

void Grid::update_E_x(unsigned long long I_s, unsigned long long I, unsigned long long J, unsigned long long K,
                      vector<double> &H_z, vector<double> &H_y, vector<double> &E_x, vector<double> &C_a, vector<double> &C_b) {
    for (unsigned long i = I_s; i < I-1; ++i) {
        for (unsigned long j = 1; j < J-1; ++j) {
            for (unsigned long k = 1; k < K-1; ++k) {
                E_x[get_ind(i, j, k)] = C_a[get_ind(i, j, k)] * E_x[get_ind(i, j, k)] +
                C_b[get_ind(i, j, k)] * ((H_z[get_ind(i, j, k)] - H_z[get_ind(i, j-1, k)])/ds -
                (H_y[get_ind(i, j, k)] - H_y[get_ind(i, j, k-1)])/ds);
            }
        }
    }
}

void Grid::update_E_y(unsigned long long I_s, unsigned long long I, unsigned long long J, unsigned long long K,
                      vector<double> &H_z, vector<double> &H_y, vector<double> &E_x, vector<double> &C_a, vector<double> &C_b) {
    for (unsigned long i = I_s + 1; i < I-1; ++i) {
        for (unsigned long j = 0; j < J-1; ++j) {
            for (unsigned long k = 1; k < K-1; ++k) {
                E_x[get_ind(i, j, k)] = C_a[get_ind(i, j, k)] * E_x[get_ind(i, j, k)] +
                C_b[get_ind(i, j, k)] * ((H_z[get_ind(i, j, k)] - H_z[get_ind(i, j, k-1)])/ds -
                (H_y[get_ind(i, j, k)] - H_y[get_ind(i-1, j, k)])/ds);
            }
        }
    }
}

void Grid::update_E_z(unsigned long long I_s, unsigned long long I, unsigned long long J, unsigned long long K,
                      vector<double> &H_z, vector<double> &H_y, vector<double> &E_x, vector<double> &C_a, vector<double> &C_b) {
    for (unsigned long i = I_s + 1; i < I-1; ++i) {
        for (unsigned long j = 1; j < J-1; ++j) {
            for (unsigned long k = 0; k < K-1; ++k) {
                E_x[get_ind(i, j, k)] = C_a[get_ind(i, j, k)] * E_x[get_ind(i, j, k)] +
                C_b[get_ind(i, j, k)] * ((H_z[get_ind(i, j, k)] - H_z[get_ind(i-1, j, k)])/ds -
                (H_y[get_ind(i, j, k)] - H_y[get_ind(i, j-1, k)])/ds);
            }
        }
    }
}

void Grid::parallel_update_E_field() {
    vector<thread> threads;

    threads.push_back(thread(&Grid::update_E_x, this, 0, this->I, this->J, this->K, ref(this->H_z), ref(this->H_y),
                      ref(this->E_x), ref(this->C_a), ref(this->C_b)));
    threads.push_back(thread(&Grid::update_E_y, this, 0, this->I, this->J, this->K, ref(this->H_x), ref(this->H_z),
                      ref(this->E_y), ref(this->C_a), ref(this->C_b)));
    threads.push_back(thread(&Grid::update_E_z, this, 0, this->I, this->J, this->K, ref(this->H_y), ref(this->H_x),
                      ref(this->E_z), ref(this->C_a), ref(this->C_b)));
    for (thread &th : threads) {
        th.join();
    }
}

void Grid::update_H_x(unsigned long long I_s, unsigned long long I, unsigned long long J, unsigned long long K,
                      vector<double> &H_z, vector<double> &H_y, vector<double> &E_x, vector<double> &C_a, vector<double> &C_b) {
    for (unsigned long i = I_s; i < I; ++i) {
        for (unsigned long j = 0; j < J-1; ++j) {
            for (unsigned long k = 0; k < K-1; ++k) {
                E_x[get_ind(i, j, k)] = C_a[get_ind(i, j, k)] * E_x[get_ind(i, j, k)] +
                C_b[get_ind(i, j, k)] * ((H_z[get_ind(i, j, k+1)] - H_z[get_ind(i, j, k)])/ds -
                (H_y[get_ind(i, j+1, k)] - H_y[get_ind(i, j, k)])/ds);
            }
        }
    }
}

void Grid::update_H_y(unsigned long long I_s, unsigned long long I, unsigned long long J, unsigned long long K,
                      vector<double> &H_z, vector<double> &H_y, vector<double> &E_x, vector<double> &C_a, vector<double> &C_b) {
    for (unsigned long i = I_s; i < I-1; ++i) {
        for (unsigned long j = 0; j < J; ++j) {
            for (unsigned long k = 0; k < K-1; ++k) {
                E_x[get_ind(i, j, k)] = C_a[get_ind(i, j, k)] * E_x[get_ind(i, j, k)] +
                C_b[get_ind(i, j, k)] * ((H_z[get_ind(i+1, j, k)] - H_z[get_ind(i, j, k)])/ds -
                (H_y[get_ind(i, j, k+1)] - H_y[get_ind(i, j, k)])/ds);
            }
        }
    }
}

void Grid::update_H_z(unsigned long long I_s, unsigned long long I, unsigned long long J, unsigned long long K,
                      vector<double> &H_z, vector<double> &H_y, vector<double> &E_x, vector<double> &C_a, vector<double> &C_b) {
    for (unsigned long i = I_s; i < I-1; ++i) {
        for (unsigned long j = 0; j < J-1; ++j) {
            for (unsigned long k = 0; k < K; ++k) {
                E_x[get_ind(i, j, k)] = C_a[get_ind(i, j, k)] * E_x[get_ind(i, j, k)] +
                C_b[get_ind(i, j, k)] * ((H_z[get_ind(i, j+1, k)] - H_z[get_ind(i, j, k)])/ds -
                (H_y[get_ind(i+1, j, k)] - H_y[get_ind(i, j, k)])/ds);
            }
        }
    }
}

void Grid::parallel_update_H_field() {
    vector<thread> threads;

    threads.push_back(thread(&Grid::update_H_x, this, 0, this->I, this->J, this->K, ref(this->E_y), ref(this->E_z),
                      ref(this->H_x), ref(this->D_a), ref(this->D_b)));
    threads.push_back(thread(&Grid::update_H_y, this, 0, this->I, this->J, this->K, ref(this->E_z), ref(this->E_x),
                      ref(this->H_y), ref(this->D_a), ref(this->D_b)));
    threads.push_back(thread(&Grid::update_H_z, this, 0, this->I, this->J, this->K, ref(this->E_x), ref(this->E_y),
                      ref(this->H_z), ref(this->D_a), ref(this->D_b)));
    for (thread &th : threads) {
        th.join();
    }
}


void Grid::serial_update_boundary_all() {
    //x face
    update_boundary_x(this->J-1, this->K, 0, this->E_y_x0, this->E_y, this->abs_fact);
    update_boundary_x(this->J, this->K-1, 0, this->E_z_x0, this->E_z, this->abs_fact);
    update_boundary_x(this->J-1, this->K, this->I-1, this->E_y_x1, this->E_y, this->abs_fact);
    update_boundary_x(this->J, this->K-1, this->I-1, this->E_z_x1, this->E_z, this->abs_fact);

    //y face
    update_boundary_y(this->I-1, this->K, 0, this->E_x_y0, this->E_x, this->abs_fact);
    update_boundary_y(this->I, this->K-1, 0, this->E_z_y0, this->E_z, this->abs_fact);
    update_boundary_y(this->I-1, this->K, this->J-1, this->E_x_y1, this->E_x, this->abs_fact);
    update_boundary_y(this->I, this->K-1, this->J-1, this->E_z_y1, this->E_z, this->abs_fact);

    //z face
    update_boundary_z(this->I-1, this->J, 0, this->E_x_z0, this->E_x, this->abs_fact);
    update_boundary_z(this->I, this->J-1, 0, this->E_y_z0, this->E_y, this->abs_fact);
    update_boundary_z(this->I-1, this->J, this->K-1, this->E_x_z1, this->E_x, this->abs_fact);
    update_boundary_z(this->I, this->J-1, this->K-1, this->E_y_z1, this->E_y, this->abs_fact);

}

void Grid::serial_update_E_field() {

    update_E_x(0, this->I, this->J, this->K, this->H_z, this->H_y, this->E_x, this->C_a, this->C_b);
    update_E_y(0, this->I, this->J, this->K, this->H_x, this->H_z, this->E_y, this->C_a, this->C_b);
    update_E_z(0, this->I, this->J, this->K, this->H_y, this->H_x, this->E_z, this->C_a, this->C_b);

}

void Grid::serial_update_H_field() {
    update_H_x(0, this->I, this->J, this->K, this->E_y, this->E_z, this->H_x, this->D_a, this->D_b);
    update_H_y(0, this->I, this->J, this->K, this->E_z, this->E_x, this->H_y, this->D_a, this->D_b);
    update_H_z(0, this->I, this->J, this->K, this->E_x, this->E_y, this->H_z, this->D_a, this->D_b);

}
