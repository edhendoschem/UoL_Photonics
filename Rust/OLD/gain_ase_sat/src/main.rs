
extern crate gnuplot;
use gnuplot::*;
use std::{f64, i64};
use std::f64::consts;

const steps: i32 = 1000; //Number of steps to propagate the signals
const nloops: i32 = 20; //Number of times to run forward and backward propagation
const h: f64 = 6.626070e-34; //Planck's constant, J.s
const c: f64 =  299792458.0; //Speed of light, m/s

struct starting_params{
    signal_freq: f64, // Signal frequency, 1/s
    signal_wavelength: f64, // Signal Wavelength, m
    pump_freq: f64, // Pump frequency, 1/s
    pump_wavelength: f64, // Pump wavelength, m
    channel_spacing: f64, // Channel spacing 125 GHz = 1 nm,  1/s
    sigma_s_a: f64, // Signal absorption cross section, m^2
    sigma_s_e: f64, // Signal emission cross section, m^2
    sigma_p_a: f64, // Pump absorption cross section, m^2
    sigma_p_e: f64, // Pump emission cross section, m^2
    tau_21: f64, // Lifetime from 13/2 to 15/2, s
    gamma_s: f64, // Overlap factor for signal and ASE, Adimensional
    gamma_p: f64, // Overlap factor for pump, Adimensional
    width: f64, // Waveguide width, m
    depth: f64, // Waveguide depth, m
    length: f64, // Waveguide length, m
    A: f64, // Area, m^2
    N: f64, // Erbium ion density, ions/m^3
    Pp0: f64, // Pump power at length 0, Watts
    Ps0: f64, // Signal power at length 0, Watts
}


//This function takes wavelength in m and outputs frequency in 1/s
fn find_freq(wavelength: &f64) -> f64 {
    c / wavelength
}

//Inputs width and depth in m and outputs area in m^2
fn find_area(width: &f64, depth: &f64) -> f64 {
    width * depth
}

//Calculates gain in dB based on output power vs input power in Watts (or other
//units)
fn find_gain(output: &f64, input: &f64) -> f64 {
    10.0 * (output / input).log(10.0)
}

//Calculates Normalized Excited Ion density (between 0 to 1)
fn N2_norm(A: &f64, B: &f64, C: &f64, D: &f64, E: &f64,
    F: &f64, Ps: &f64, PA: &f64, Pp: &f64) -> f64 {
        (A * Ps + B * PA + C * Pp) / (D * Ps + E * PA + F * Pp + 1.0)
    }

//Calculates dPump power/dz, G is sigma_p_a + sigma_p_e
fn dPp(N: &f64, G: &f64, sigma_p_a: &f64,
    gamma_p: &f64, N2_norm: &f64, Pp: &f64) -> f64 {
    N * (N2_norm * G - sigma_p_a) * gamma_p * Pp
}

//Calculates dSignal power/dz, H is sigma_s_a + sigma_s_e
fn dPs(N: &f64, H: &f64, sigma_s_a: &f64,
    gamma_s: &f64, N2_norm: &f64, Ps: &f64) -> f64 {
    N * (N2_norm * H - sigma_s_a) * gamma_s * Ps
}

//Calculates dForward ASE/dz, I is sigma_s_e * h * signal_freq * channel_spacing
fn dPAf(N: &f64, H: &f64, sigma_s_a: &f64,
    gamma_s: &f64, I: &f64, N2_norm: &f64, PAf: &f64) -> f64 {
    N * (N2_norm * H - sigma_s_a) * gamma_s * PAf + N * N2_norm * I
}

//Calculates dBackward ASE/dz
fn dPAb(N: &f64, H: &f64, sigma_s_a: &f64,
    gamma_s: &f64, I: &f64, N2_norm: &f64, PAb: &f64) -> f64 {
    -1.0 * (N * (N2_norm * H - sigma_s_a) * gamma_s * PAb + N * N2_norm * I)
}

fn main() {
    //Gain for 10 mW pump no ASE = 10.8260
    let mut params: starting_params = starting_params {
        signal_freq: 0.0, //This value is later calculated Signal frequency, 1/s
        signal_wavelength: 1536e-9, // Signal Wavelength, m
        pump_freq: 0.0, // This values is calculated later, Pump frequency, 1/s
        pump_wavelength: 981e-9, // Pump wavelength, m
        channel_spacing: (125.0e9 * 70.0), // Channel spacing, 1/s
        sigma_s_a: 5.8e-25, // Signal absorption cross section, m^2
        sigma_s_e: 7.27e-25, // Signal emission cross section, m^2
        sigma_p_a: 3.12e-25, // Pump absorption cross section, m^2
        sigma_p_e: 0.0, // Pump emission cross section no emission assumed, m^2
        tau_21: 12.94e-3, // Lifetime from 13/2 to 15/2, s
        gamma_s: 0.693, // Overlap factor for signal and ASE, Adimensional (from Marcatilis)
        gamma_p: 0.863, // Overlap factor for pump, Adimensional (from MArcatilis)
        width: 5.0e-6, // Waveguide width, m
        depth: 5.0e-6, // Waveguide depth, m
        length: 10.0e-2, // Waveguide length, m
        A: 0.0, //This will be calculated later Area, m^2
        N: 4.45e26, // Erbium ion density, ions/m^3
        Pp0: 10.0e-3, // 10 mW Pump power at length 0, Watts
        Ps0: 1.0e-6, // -40 dBm (1e-7 W)Signal power at length 0, Watts
    };

    //Calculating and storing the signal frequency
    params.signal_freq = find_freq(&params.signal_wavelength);

    //Calculating and storing the pump frequency
    params.pump_freq = find_freq(&params.pump_wavelength);

    //Calculating and storing the Area
    params.A = find_area(&params.width, &params.depth);

    //Definition of equation constants
    let A: f64 = params.tau_21 * params.sigma_s_a * params.gamma_s /
                (params.A * h * params.signal_freq);

    //Same parameters as signal therefore same as signal
    let B: f64 = A.clone();

    let C: f64 = params.tau_21 * params.sigma_p_a * params.gamma_p /
                (params.A * h * params.pump_freq);

    let D: f64 = params.tau_21 * (params.sigma_s_a + params.sigma_s_e)
                * params.gamma_s / (params.A * h * params.signal_freq);

    //Same parameters as signal therefore same as signal
    let E: f64 = D.clone();

    let F: f64 = params.tau_21 * (params.sigma_p_a + params.sigma_p_e)
                * params.gamma_p / (params.A * h * params.pump_freq);

    let G: f64 = params.sigma_p_e + params.sigma_p_a;

    let H: f64 = params.sigma_s_e + params.sigma_s_a;

    let I: f64 = params.sigma_s_e * params.gamma_s * h * params.signal_freq
                * params.channel_spacing;

    let N: f64 = params.N; //Erbium Ion density

    let gamma_s: f64 = params.gamma_s;
    let gamma_p: f64 = params.gamma_p;
    let sigma_s_a: f64 = params.sigma_s_a;
    let sigma_s_e: f64 = params.sigma_s_e;
    let sigma_p_a: f64 = params.sigma_p_a;
    let sigma_p_e: f64 = params.sigma_p_e;


    //Defining initial conditions
    //let mut Pp0: f64 = 10e-3;
    let mut Pp0: f64 = params.Pp0; //Initial pump power
    let mut Ps0: f64 = params.Ps0; //Initial signal power
    let mut PAf0: f64 = 0.0; //Initial Forward ASE at L = 0
    let mut PAb0: f64 = 0.0; //Initial Backwards ASE at L = 0, incorrect guess,
    //at L = 0, backwards ASE is max
    let delta_z: f64 = params.length / (steps as f64); //Number of divisions
    let mut N20 = N2_norm(&A, &B, &C, &D, &E, &F, &Ps0, &0.0, &Pp0); //Starting N2
    let mut z_tot: f64 = 0.0; //current position
    let mut N21: f64 = 0.0; // Next N2 value
    let input: f64 = Ps0.clone(); //Reference value used to calculate gain
    let mut dir: f64 = 1.0; //Direction of propagation, Forward +1.0, backwards -1.0

println!("A = {}, B = {}, C = {}, D = {}, E = {}, F = {}, Ps0 = {}, PA = {}, Pp0 = {}",
    A, B, C, D, E, F, Ps0, 0, Pp0);


///////////////////////////////////////////////////////////////////////////////
    //The order is N2, Ps, Pp, PAf, PAb
    //A vector containing all the positions of z
    let mut z_position: [f64; (steps) as usize] = [0.0; (steps) as usize];
    //A vector containing the value of Ps with respect of z
    let mut Ps_position: [f64; (steps) as usize] = [0.0; (steps) as usize];
    Ps_position[0] = Ps0;
    //A vector containing the value of Pp with respect of z
    let mut Pp_position: [f64; (steps) as usize] = [0.0; (steps) as usize];
    Pp_position[0] = Pp0;
    //A vector containing the value of PAf with respect of z
    let mut PAf_position: [f64; (steps) as usize] = [0.0; (steps) as usize];
    PAf_position[0] = PAf0;
    //A vector containing the value of PAb with respect of z
    let mut PAb_position: [f64; (steps) as usize] = [0.0; (steps) as usize];
    PAb_position[0] = PAb0;
    //A vector containing the value of N2 with respect of z
    let mut N2_position: [f64; (steps) as usize] = [0.0; (steps) as usize];
    N2_position[0] = N20;

    let mut PAb1: f64 = 0.0;
    let mut PAf1: f64 = 0.0;
    let mut Pp1: f64 = 0.0;
    let mut Ps1: f64 = 0.0;


///////////////////////////////////////////////////////////////////////////////
//Final Runge Kutta loop

for n in 0..nloops {
    //Initial Approximation
        for i in 0..(steps) {
            //4th Order Runge Kutta for Signal
            //dPs(N: &f64, H: &f64, sigma_s_a: &f64, gamma_s: &f64, N2_norm: &f64 Ps: &f64)
            //z = 0
            let k1_s: f64 = dir * dPs(&N, &H, &sigma_s_a, &gamma_s, &N20, &Ps0);

            let k2_s: f64 = dir * dPs(&N, &H, &sigma_s_a, &gamma_s, &N20,
                            &(Ps0 + k1_s * delta_z / 2.0));
            let k3_s: f64 = dir * dPs(&N, &H, &sigma_s_a, &gamma_s, &N20,
                            &(Ps0 + k2_s * delta_z / 2.0));
            let k4_s: f64 = dir * dPs(&N, &H, &sigma_s_a, &gamma_s, &N20,
                            &(Ps0 + k3_s * delta_z));
            //Approximation of the value at z+delta_z
            Ps1 = Ps0 + (1.0/6.0) * (k1_s + 2.0 * k2_s + 2.0 * k3_s + k4_s) * delta_z;
            if Ps1 < 0.0 {Ps1 = 0.0}

            //4th Order Runge Kutta for Pump
            //dPp(N: &f64, G: &f64, sigma_p_a: &f64, gamma_p: &f64, N2_norm: &f64, Pp: &f64)
            //z = 0

            let k1_p: f64 = dir * dPp(&N, &G, &sigma_p_a, &gamma_p, &N20, &Pp0);

            let k2_p: f64 = dir * dPp(&N, &G, &sigma_p_a, &gamma_p, &N20,
                            &(Pp0 + k1_p * delta_z / 2.0));
            let k3_p: f64 = dir * dPp(&N, &G, &sigma_p_a, &gamma_p, &N20,
                            &(Pp0 + k2_p * delta_z / 2.0));
            let k4_p: f64 = dir * dPp(&N, &G, &sigma_p_a, &gamma_p, &N20,
                            &(Pp0 + k3_p * delta_z));
            //Approximation of the value at z+delta_z
            Pp1 = Pp0 + (1.0/6.0) * (k1_p + 2.0 * k2_p + 2.0 * k3_p + k4_p) * delta_z;

            if Pp1 < 0.0 {Pp1 = 0.0}
            if Pp1 > params.Pp0 {Pp1 = params.Pp0}
            //4th Order Runge Kutta for Forward ASE
            //dPAf(N: &f64, H: &f64, sigma_s_a: &f64, gamma_s: &f64, I: &f64, N2_norm: &f64, PAf: &f64)
            //z = 0
            let k1_Af: f64 = dir * dPAf(&N, &H, &sigma_s_a, &gamma_s, &I, &N20, &PAf0);

            let k2_Af: f64 = dir * dPAf(&N, &H, &sigma_s_a, &gamma_s, &I, &N20,
                            &(PAf0 + k1_Af * delta_z / 2.0));
            let k3_Af: f64 = dir * dPAf(&N, &H, &sigma_s_a, &gamma_s, &I, &N20,
                            &(PAf0 + k2_Af * delta_z / 2.0));
            let k4_Af: f64 = dir * dPAf(&N, &H, &sigma_s_a, &gamma_s, &I, &N20,
                            &(PAf0 + k3_Af * delta_z));
            //Approximation of the value at z+delta_z
            PAf1 = PAf0 + (1.0/6.0) * (k1_Af + 2.0 * k2_Af + 2.0 * k3_Af + k4_Af) * delta_z;

            if PAf1 < 0.0 {PAf1 = 0.0}
            //4th Order Runge Kutta for Backwards ASE
            //dPAb(N: &f64, H: &f64, sigma_s_a: &f64, gamma_s: &f64, I: &f64, N2_norm: &f64, PAb: &f64)
            //z = 0
            let k1_Ab: f64 = dir * dPAb(&N, &H, &sigma_s_a, &gamma_s, &I, &N20, &PAb0);

            let k2_Ab: f64 = dir * dPAb(&N, &H, &sigma_s_a, &gamma_s, &I, &N20,
                                &(PAb0 + k1_Ab * delta_z / 2.0));
            let k3_Ab: f64 = dir * dPAb(&N, &H, &sigma_s_a, &gamma_s, &I, &N20,
                                &(PAb0 + k2_Ab * delta_z / 2.0));
            let k4_Ab: f64 = dir * dPAb(&N, &H, &sigma_s_a, &gamma_s, &I, &N20,
                                &(PAb0 + k3_Ab * delta_z));
            //Approximation of the value at z+delta_z
            PAb1 = PAb0 + (1.0/6.0) * (k1_Ab + 2.0 * k2_Ab + 2.0 * k3_Ab + k4_Ab) * delta_z;

            if PAb1 < 0.0 {PAb1 = 0.0}
            //N2_norm(A: &f64, B: &f64, C: &f64, D: &f64, E: &f64, F: &f64, Ps: &f64, PA: &f64, Pp: &f64)
            N21 = N2_norm(&A, &B, &C, &D, &E, &F, &Ps1, &(PAf1+PAb1), &Pp1);
            N20 = N21;
            Pp0 = Pp1;
            Ps0 = Ps1;
            PAf0 = PAf1;
            PAb0 = PAb1;
            z_tot += delta_z;
            z_position[i as usize] = z_tot;
            Pp_position[i as usize] = Pp1;
            Ps_position[i as usize] = Ps1;
            PAf_position[i as usize] = PAf1;
            PAb_position[i as usize] = PAb1;
            N2_position[i as usize] = N21;
            /*
            println!("Ps0 = {}, Pp0 = {}, PAf0 = {}, PAb0 = {}", Ps0, Pp0, PAf0, PAb0);
            println!("Gain {}", find_gain(&Ps0, &input));
            println!("___________________________________");
            */
        }
    ///////////////////////////////////////////////////////////////////////////////
        //Second Approximation

        N2_position.reverse();
        N20 = N2_position[0];
        Ps_position.reverse();
        Ps0 = Ps_position[0];
        Pp_position.reverse();
        Pp0 = Pp_position[0];
        PAf_position.reverse();
        PAf0 = PAf_position[0];
        PAb_position.reverse();
        PAb_position[0] = 0.0; //Restoring boundary conditions
        PAb0 = PAb_position[0];
        z_tot = 0.0;
        dir = -1.0;

        for i in 0..(steps) {
            //4th Order Runge Kutta for Signal
            //dPs(N: &f64, H: &f64, sigma_s_a: &f64, gamma_s: &f64, N2_norm: &f64 Ps: &f64)
            //z = 0
            let k1_s: f64 = dir * dPs(&N, &H, &sigma_s_a, &gamma_s, &N20, &Ps0);

            let k2_s: f64 = dir * dPs(&N, &H, &sigma_s_a, &gamma_s, &N20,
                            &(Ps0 + k1_s * delta_z / 2.0));
            let k3_s: f64 = dir * dPs(&N, &H, &sigma_s_a, &gamma_s, &N20,
                            &(Ps0 + k2_s * delta_z / 2.0));
            let k4_s: f64 = dir * dPs(&N, &H, &sigma_s_a, &gamma_s, &N20,
                            &(Ps0 + k3_s * delta_z));
            //Approximation of the value at z+delta_z
            Ps1 = Ps0 + (1.0/6.0) * (k1_s + 2.0 * k2_s + 2.0 * k3_s + k4_s) * delta_z;

            if Ps1 < 0.0 {Ps1 = 0.0}
            //4th Order Runge Kutta for Pump
            //dPp(N: &f64, G: &f64, sigma_p_a: &f64, gamma_p: &f64, N2_norm: &f64, Pp: &f64)
            //z = 0
            let k1_p: f64 = dir * dPp(&N, &G, &sigma_p_a, &gamma_p, &N20, &Pp0);

            let k2_p: f64 = dir * dPp(&N, &G, &sigma_p_a, &gamma_p, &N20,
                            &(Pp0 + k1_p * delta_z / 2.0));
            let k3_p: f64 = dir * dPp(&N, &G, &sigma_p_a, &gamma_p, &N20,
                            &(Pp0 + k2_p * delta_z / 2.0));
            let k4_p: f64 = dir * dPp(&N, &G, &sigma_p_a, &gamma_p, &N20,
                            &(Pp0 + k3_p * delta_z));
            //Approximation of the value at z+delta_z
            Pp1 = Pp0 + (1.0/6.0) * (k1_p + 2.0 * k2_p + 2.0 * k3_p + k4_p) * delta_z;

            if Pp1 < 0.0 {Pp1 = 0.0}
            if Pp1 > params.Pp0 {Pp1 = params.Pp0}
            //4th Order Runge Kutta for Forward ASE
            //dPAf(N: &f64, H: &f64, sigma_s_a: &f64, gamma_s: &f64, I: &f64, N2_norm: &f64, PAf: &f64)
            //z = 0
            let k1_Af: f64 = dir * dPAf(&N, &H, &sigma_s_a, &gamma_s, &I, &N20, &PAf0);

            let k2_Af: f64 = dir * dPAf(&N, &H, &sigma_s_a, &gamma_s, &I, &N20,
                            &(PAf0 + k1_Af * delta_z / 2.0));
            let k3_Af: f64 = dir * dPAf(&N, &H, &sigma_s_a, &gamma_s, &I, &N20,
                            &(PAf0 + k2_Af * delta_z / 2.0));
            let k4_Af: f64 = dir * dPAf(&N, &H, &sigma_s_a, &gamma_s, &I, &N20,
                            &(PAf0 + k3_Af * delta_z));
            //Approximation of the value at z+delta_z
            let mut PAf1: f64 = PAf0 +
                        (1.0/6.0) * (k1_Af + 2.0 * k2_Af + 2.0 * k3_Af + k4_Af) * delta_z;

            if PAf1 < 0.0 {PAf1 = 0.0}
            //4th Order Runge Kutta for Backwards ASE
            //dPAb(N: &f64, H: &f64, sigma_s_a: &f64, gamma_s: &f64, I: &f64, N2_norm: &f64, PAb: &f64)
            //z = 0
            let k1_Ab: f64 = dir * dPAb(&N, &H, &sigma_s_a, &gamma_s, &I, &N20, &PAb0);

            let k2_Ab: f64 = dir * dPAb(&N, &H, &sigma_s_a, &gamma_s, &I, &N20,
                                &(PAb0 + k1_Ab * delta_z / 2.0));
            let k3_Ab: f64 = dir * dPAb(&N, &H, &sigma_s_a, &gamma_s, &I, &N20,
                                &(PAb0 + k2_Ab * delta_z / 2.0));
            let k4_Ab: f64 = dir * dPAb(&N, &H, &sigma_s_a, &gamma_s, &I, &N20,
                                &(PAb0 + k3_Ab * delta_z));
            //Approximation of the value at z+delta_z
            PAb1 = PAb0 + (1.0/6.0) * (k1_Ab + 2.0 * k2_Ab + 2.0 * k3_Ab + k4_Ab) * delta_z;

            if PAb1 < 0.0 {PAb1 = 0.0}
            //N2_norm(A: &f64, B: &f64, C: &f64, D: &f64, E: &f64, F: &f64, Ps: &f64, PA: &f64, Pp: &f64)
            N21 = N2_norm(&A, &B, &C, &D, &E, &F, &Ps1, &(PAf1+PAb1), &Pp1);
            N20 = N21;
            Pp0 = Pp1;
            Ps0 = Ps1;
            PAf0 = PAf1;
            PAb0 = PAb1;
            z_tot += delta_z;
            z_position[i as usize] = z_tot;
            Pp_position[i as usize] = Pp1;
            Ps_position[i as usize] = Ps1;
            PAf_position[i as usize] = PAf1;
            PAb_position[i as usize] = PAb1;
            N2_position[i as usize] = N21;
            /*
            println!("Ps0 = {}, Pp0 = {}, PAf0 = {}, PAb0 = {}", Ps0, Pp0, PAf0, PAb0);
            println!("Gain {}", find_gain(&Ps0, &input));
            println!("___________________________________");
            */
        }

//Restoring boundary conditions

        Ps_position.reverse();
        Ps_position[0] = params.Ps0;
        Ps0 = Ps_position[0];
        Pp_position.reverse();
        Pp0 = params.Pp0; //Restoring Original pump power = 0
        Pp_position[0] = Pp0;
        PAf_position.reverse();//We use this value as our starting point
        PAf0 = 0.0;
        PAf_position[0] = 0.0;
        PAb_position.reverse(); //We use this as our starting point
        PAb0 = PAb_position[0];
        N2_position.reverse();
        N20 = N2_norm(&A, &B, &C, &D, &E, &F, &Ps0, &(PAf0+PAb0), &Pp0);
        N2_position[0] = N20;
        z_tot = 0.0;
        dir = 1.0;
    println!("Loop number = {}", n);
    println!("Gain {}", find_gain(&Ps_position[Ps_position.len() - 1], &input));
    println!("___________________________________");
}


//////////////////////////////////////////////////////////////////////////////
//Plotting the graphs
/*
    println!("z_position = {:?}, length = {}", z_position, z_position.len());
    println!("N2_position = {:?}, length = {}", N2_position, N2_position.len());
    println!("Ps_position = {:?}, length = {}", Ps_position, Ps_position.len());
    println!("Pp_position = {:?}, length = {}", Pp_position, Pp_position.len());
    println!("PAf_position = {:?}, length = {}", PAf_position, PAf_position.len());
    println!("PAb_position = {:?}, length = {}", PAb_position, PAb_position.len());
*/
    //This loop helps adjust values for displaying in the graphs
    for i in 0..steps {
        N2_position[i as usize] = N2_position[i as usize] * 1.0;
        Ps_position[i as usize] = Ps_position[i as usize] * 1000.0;
        Pp_position[i as usize] = Pp_position[i as usize] * 1000.0;
        PAf_position[i as usize] = PAf_position[i as usize] * 1000.0;
        PAb_position[i as usize] = PAb_position[i as usize] * 1000.0;
        z_position[i as usize] = z_position[i as usize] * 100.0;

    }


    let mut fg1 = Figure::new();
    fg1.axes2d()
    .lines(z_position.iter(), PAf_position.iter(), &[Caption("Forward ASE"), Color("blue")])
    .lines(z_position.iter(), PAb_position.iter(), &[Caption("Backwards ASE"), Color("green")])
    .set_x_label("Z position in cm", &[])
    .set_y_label("Power in mW",&[]);
    fg1.show();

    let mut fg2 = Figure::new();
    fg2.axes2d()
    .lines(z_position.iter(), Pp_position.iter(), &[Caption("Pump Power"), Color("black")])
    .set_x_label("Z position in cm", &[])
    .set_y_label("Power in mW",&[]);
    fg2.show();

    let mut fg3 = Figure::new();
    fg3.axes2d()
    .lines(z_position.iter(), Ps_position.iter(), &[Caption("Signal Power"), Color("red")])
    .set_x_label("Z position in cm", &[])
    .set_y_label("Power in mW",&[]);
    fg3.show();

    let mut fg4 = Figure::new();
    fg4.axes2d()
    .lines(z_position.iter(), N2_position.iter(), &[Caption("Excited ions"), Color("yellow")])
    .set_x_label("Z position in cm", &[])
    .set_y_label("Fractional Population Inversion",&[]);
    fg4.show();

}
