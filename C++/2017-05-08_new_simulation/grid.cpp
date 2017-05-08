/*--------------Grid 1D TM, #3--------------*/
#include "definitions.h"

ostream &operator << (ostream &os, Point &a) {
    cout<<"x = "<<a.x<<'\t'<<"y = "<<a.y<<'\n';
    return os;
}

Grid1DTM::Grid1DTM(unsigned long long max_timee, unsigned long long lengthh, double wavelengthh, double courantt, double ppww,
unsigned long long source_nodee, int typee)
: max_time {max_timee}, length {lengthh}, wavelength {wavelengthh}, courant {courantt}, ppw {ppww}, Ez {lengthh}, Hy {lengthh},
c1ez {lengthh}, c2ez {lengthh}, c1hy {lengthh}, c2hy {lengthh}, source_node {source_nodee}, type {typee}, Ez_left{6},
Ez_right{6} {
    curr_time = 0;
    dx = wavelength / ppw;
    dt = courant * dx / c0;
    min_wavelength = wavelength;
    time_delay = (ppw - source_node) / courant;

    for (unsigned long long i = 0; i<c1ez.size(); ++i) {
        c1ez(i) = 1.0;
        c2ez(i) = courant * imp0;
    }

    for (unsigned long long i = 0; i<c1hy.size(); ++i) {
        c1hy(i) = 1.0;
        c2hy(i) = courant / imp0;
    }

}


void Grid1DTM::show_params() const {
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

void Grid1DTM::change_courant(double new_val) {
    cout<<"Warning: Changing the Courant number from default values may lead to numerical instability\n";
    cout<<"Previous value courant = "<<courant<<'\n';
    courant = new_val;
    cout<<"Changed value to courant = "<<courant<<"\n Note: dx and dt adjusted to new Courant number"<<'\n';
    cout<<"======================================================================\n";
    dx = wavelength / ppw;
    dt = courant * dx / c0;
    return;
}

void Grid1DTM::change_dx(double new_val) {
    cout<<"Warning: Changing dx may lead to numerical instability\n";
    cout<<"Previous value dx = "<<dx<<'\n';
    dx = new_val;
    cout<<"Changed value to dx = "<<dx<<"\nNote: dt and ppw adjusted to new dx value to preserve Courant number"<<'\n';
    cout<<"======================================================================\n";
    ppw = wavelength / dx;
    dt = courant * dx / c0;
    return;
}

void Grid1DTM::change_dt(double new_val) {
    cout<<"Warning: Changing dt may lead to numerical instability\n";
    cout<<"Previous value dt = "<<dt<<'\n';
    dt = new_val;
    cout<<"Changed value to dt = "<<dt<<"\nNote: dx and ppw adjusted to new dt value to preserve Courant number"<<'\n';
    cout<<"======================================================================\n";
    dx = dt * c0 / courant;
    ppw = wavelength / dx;
    return;
}

void Grid1DTM::change_ppw(double new_val) {
    cout<<"Warning: A ppw value below the default may lead to numerical instability\n";
    cout<<"Previous value ppw = "<<ppw<<'\n';
    ppw = new_val;
    cout<<"Changed value to ppw = "<<ppw<<"\nNote: dx and dt adjusted to new ppw value to preserve Courant number"<<'\n';
    cout<<"======================================================================\n";
    dx = wavelength / ppw;
    dt = courant * dx / c0;
    return;
}

double Grid1DTM::source() {
    switch(type) {
    case 0: //Harmonic source
        return sin(2.0 * PI * (courant * (curr_time-time_delay) - source_node) / ppw);
    case 1: //Gaussian source
        return exp(-pow(((curr_time-time_delay) - 2.2 * dispersion)/ dispersion, 2)); //exp(-pow((curr_time - time_delay)/ dispersion, 2));
    case 2: //Ricker wavelet travelling in the positive x direction
        { //Scope delimiters to prevent compiler errors
        double a = (courant * (curr_time - time_delay) - source_node) / ppw;
        double b = pow(a - 2.0, 2); //(a - Md)^2
        double c = 1.0 - 2.0 * PI * PI * b;
        double d = exp(-1.0 * PI * PI * b);
        double e = c * d;
        return e;
        }
    case 3: //Gaussian pulse
        { //Scope delimiters to prevent compiler errors
        double f = sin(2.0 * PI * (courant * (curr_time-time_delay) - source_node) / ppw); //Harmonic pulse
        double g = exp(-pow(((curr_time-time_delay) - 2.2 * dispersion)/ dispersion, 2)); //Gaussian envelope
        return f * g;
        }
    //default: No source
    default:
        return 0.0;

    }
}

void Grid1DTM::update_magnetic() {
    for (unsigned long long i = 0; i < Hy.size()-1; ++i) {
        Hy(i) = c1hy(i) * Hy(i) + c2hy(i) * (Ez(i + 1) - Ez(i));
    }

    Hy(source_node) -= source() / imp0;
}

void Grid1DTM::update_electric() {
    for (unsigned long long i = 1; i < Ez.size()-1; ++i) {
        Ez(i) = c1ez(i) * Ez(i) + c2ez(i) * (Hy(i) - Hy(i - 1));
    }

    Ez(source_node) += source();
}

void Grid1DTM::update_magnetic_hs() {
    for (unsigned long long i = 0; i < Hy.size()-1; ++i) {
        Hy(i) = c1hy(i) * Hy(i) + c2hy(i) * (Ez(i + 1) - Ez(i));
    }

    Hy(source_node) = source() / imp0;
}

void Grid1DTM::update_electric_hs() {
    for (unsigned long long i = 1; i < Ez.size()-1; ++i) {
        Ez(i) = c1ez(i) * Ez(i) + c2ez(i) * (Hy(i) - Hy(i - 1));
    }

    Ez(source_node) = source();
}

void Grid1DTM::place_source(unsigned long long node) {
    source_node = node;
    return;
}

void Grid1DTM::abc_left() {
    //Left edge
    //Notes:
    //Ez_left(0) = E0 at time t-1
    //Ez_left(1) = E1 at time t-1
    //Ez_left(2) = E2 at time t-1
    //Ez_left(3) = E0 at time t
    //Ez_left(4) = E1 at time t
    //Ez_left(5) = E2 at time t
    unsigned long long i = 0;
    double cour_prime = c2ez(i) * c2hy(i);
    A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
    B = (1.0/cour_prime) - 2 + cour_prime;
    C = 2.0 * (cour_prime - 1/cour_prime);
    D = -4.0 * (cour_prime + 1/cour_prime);
    Ez(i) = A * (B * (Ez(i+2) + Ez_left(0)) + C * (Ez_left(3) + Ez_left(5) - Ez(i+1) - Ez_left(1)) + D * Ez_left(4))
    - Ez_left(2);

    //Remember old fields
    for (i = 0; i < 3; ++i) {
        Ez_left(i) = Ez_left(i+3);
        Ez_left(i+3) = Ez(i);
    }

    return;
}

void Grid1DTM::abc_right() {
    //Right edge
    //Notes:
    //Ez_right(0) = E0 at time t-1
    //Ez_right(1) = E1 at time t-1
    //Ez_right(2) = E2 at time t-1
    //Ez_light(3) = E0 at time t
    //Ez_light(4) = E1 at time t
    //Ez_light(5) = E2 at time t
    unsigned long long i = length - 1;
    double cour_prime = c2ez(i) * c2hy(i);
    A = -1.0 / ((1/cour_prime)+ 2 + cour_prime);
    B = (1.0/cour_prime) - 2 + cour_prime;
    C = 2.0 * (cour_prime - 1/cour_prime);
    D = -4.0 * (cour_prime + 1/cour_prime);
    Ez(i) = A * (B * (Ez(i-2) + Ez_right(0)) + C * (Ez_right(3) + Ez_right(5) - Ez(i-1) - Ez_right(1)) + D * Ez_right(4))
    - Ez_right(2);

    //Remember old fields
    for (i = 0; i < 3; ++i) {
        Ez_right(i) = Ez_right(i+3);
        Ez_right(i+3) = Ez(length - 1 - i);
    }

    return;
}

void Grid1DTM::update_abc() {
    abc_left();
    abc_right();
    return;
}

void Grid1DTM::advance_simulation() {
    if (curr_time == max_time) {return;}
    update_magnetic();
    update_electric();
    update_abc();
    ++curr_time;
}

void Grid1DTM::advance_simulation_hs() {
    if (curr_time == max_time) {return;}
    update_magnetic_hs();
    update_electric_hs();
    update_abc();
    ++curr_time;
}


void Grid1DTM::save_state() {
    stringstream name;
    name<<"0-"<<"snapshot_"<<curr_time<<".csv";
    ofstream file {name.str()};
    unsigned long long min_size = Ez.size() < Hy.size() ? Ez.size() : Hy.size(); //sets the smallest size
    for (unsigned long long x = 0; x<min_size; ++x) {
        file<<x<<','<<Ez(x)<<','<<Hy(x)<<'\n';
    }
    file.close();
    return;
}

unsigned long long Grid1DTM::return_maximum_time() const{
    return max_time;
}

unsigned long long Grid1DTM::return_current_time() const{
    return curr_time;
}

void Grid1DTM::enable_lossy_termination() {

    double depthInLayer, lossFactor, MAX_LOSS = 0.35;
    for (unsigned long long i = 0; i < c1hy.size() - 1; ++i) {
        if (i < c1hy.size() - 1 - 20) {
            c1ez(i) = 1.0; //Equivalent to ceze
            c2ez(i) = courant * imp0; //equivalent to cezh
            c1hy(i) = 1.0; //equivalent to chyh
            c2hy(i) = courant / imp0; //equivalent to chye
        } else {
            depthInLayer = i - (c1hy.size() - 1 - 20) + 0.5;
            lossFactor = MAX_LOSS * pow(depthInLayer / 20, 2);
            c1ez(i) = (1.0 - lossFactor) / (1.0 + lossFactor);
            c2ez(i) = courant * imp0 / (1.0 + lossFactor);
            depthInLayer += 0.5;
            lossFactor = MAX_LOSS * pow(depthInLayer / 20, 2);
            c1hy(i) = (1.0 - lossFactor) / (1.0 + lossFactor); //equivalent to
            c2hy(i) = courant / imp0 / (1.0 + lossFactor);
        }
    }
}

void Grid1DTM::add_object(unsigned long long start, unsigned long long finish, double eps_rel, double mu_rel, double sigma,
double sigma_mag) {
    double a1 = sigma * dt / (2.0 * eps_rel);
    double a2 = sigma_mag * dt / (2.0 * mu_rel);
    double b1 = imp0 * courant / eps_rel;
    double b2 = courant / (imp0 * mu_rel);

    double new_min_wavelength = wavelength / (sqrt(eps_rel * mu_rel));

    if (new_min_wavelength < min_wavelength) {
        min_wavelength = new_min_wavelength;
        double new_dx = (min_wavelength / ppw);
        change_dx(new_dx);
    }

    unsigned long long maximum_length = c1ez.size() < c1hy.size() ? c1ez.size() : c1hy.size();
    if (start > maximum_length || finish > maximum_length) {
        cout<<"start or end point cannot exceed vector length"<<'\n';
        return;
    }

    for (unsigned long long i = start; i < finish; ++i) {
            c1ez(i) = (1.0 - a1) / (1.0 + a1); //Equivalent to ceze
            c2ez(i) = b1 / (1.0 + a1); //equivalent to cezh
            c1hy(i) = (1.0 - a2) / (1.0 + a2); //equivalent to chyh
            c2hy(i) = b2 / (1.0 + a2); //equivalent to chye
    }

    return;
}


void Grid1DTM::set_time_delay(double val) {
    cout<<"Changing time delay: Previous value of time_delay = "<<time_delay<<'\n';
    time_delay = val;
    cout<<"Changed value to: time_delay = "<<time_delay<<'\n';
    cout<<"======================================================================\n";
    return;
}

void Grid1DTM::start_source_at_0() {
    cout<<"Changing time delay: Previous value of time_delay = "<<time_delay<<'\n';
    time_delay = (ppw - source_node) / courant;
    cout<<"Changed value to: time_delay = "<<time_delay<<'\n';
    cout<<"======================================================================\n";
    return;
}

void Grid1DTM::set_dispersion(double val) {
    cout<<"Changing time : Previous value of dispersion = "<<dispersion<<'\n';
    dispersion = val;
    cout<<"Changed value to: dispersion = "<<dispersion<<'\n';
    cout<<"======================================================================\n";
    return;
}

void Grid1DTM::save_information() {
    stringstream name;
    name<<"0-Information_1DTM.csv";
    ofstream file {name.str()};
    file<<"Courant = "<<courant<<'\n';
    file<<"Wavelength = "<<wavelength<<" metres"<<'\n';
    file<<"Min Wavelength = "<<min_wavelength<<" metres"<<'\n';
    file<<"ppw = "<<ppw<<" points per wavelength"<<'\n';
    file<<"dx = "<<dx<<" metres"<<'\n';
    file<<"dt = "<<dt<<" seconds"<<'\n';
    file<<"time_delay = "<<time_delay<<" dt's\n";
    file<<"dispersion = "<<dispersion<<" dt's\n";
    file.close();
    return;
}

double Grid1DTM::return_Ez(unsigned long long i) {
    return Ez(i);
}

double Grid1DTM::return_Hy(unsigned long long i) {
    return Hy(i);
}

