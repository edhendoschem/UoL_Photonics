extern crate rand;
extern crate gnuplot;

use rand::distributions::{IndependentSample, Range};
use gnuplot::*;
use std::{f64, i64};
use std::f64::consts;
use std::error::Error;
use std::str::FromStr;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

///////////////////////CONSTANTS FOR MODE CALCULATION///////////////////////////
const STEPSMODE: usize = 1000000;
//Calculation
const PI: f64 = consts::PI;
const mu_0: f64 = (4.0 * PI) * 1e-7;                                            // Permeability of free space in H/m
const epsilon_0: f64 = 8.854187817 * 1e-12;                                     // Permittivity of free space in F/m
const P: usize = 2;
const Q: usize = 2;
///////////////////////////CONSTANTS FOR GAIN CALCULATION///////////////////////
//const STEPS: usize = 1000; // Number of subdivision in the fiber
const MATCHING: u32 = 4; //Number of cycles to try and match N21 in each Runge Kutta
const CYCLES: u32 = 10; //Number of back and forth propagations
const h: f64 = 6.626070e-34; // Planck's constant, J.s
const c: f64 =  299792458.0; // Speed of light, m/s
const SIGNALS: usize = 35; // Number of signals to consider
const FPUMPS: usize = 1; // Number of forward Pump signals to consider
const BPUMPS: usize = 1; // Number of forward Pump signals to consider
const FASE: usize = 100; // Number of forward ASE signals to consider
const BASE: usize = 100; // Number of forward ASE signals to consider
const LIFETIMES: usize = 1; // Number of lifetimes to consider

///////////////////////STRUCTURES FOR MODE CALCULATION//////////////////////////
struct mode_params {
    n_0: f64,                                                                   // Refractive index of the cladding
    n_1: f64,                                                                   // Refractive index of the core
    wavelength: f64,                                                            // Wavelength in m
    frequency: f64,                                                             // Frequency in Hz
    w_ang: f64,                                                                // Angular frequency in rad/s
    k: f64,                                                                     // Wavenumber in Vacuum
    d: f64,                                                                     // depth
    a: f64,                                                                     // width
}

struct mode_vectors {
    k_x: [f64;P],
    k_y: [f64;Q],
    gamma_x: [f64; P],
    gamma_y: [f64; Q],
    beta: [[f64;P];Q],
    sols_Ex: [[f64; 13];P*Q], //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
    sols_Ey: [[f64; 13];P*Q], //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
}

////////////////////////STRUCTS FOR GAIN CALCULATION////////////////////////////

struct starting_params {
    // Signals frequency, 1/s
    s_freq: Vec<f64>,
    // Signals starting wavelength
    s_start: f64,
    // Signals Wavelength, m
    s_wavelength: Vec<f64>,
    // Signals absorption cross section [i][0] and emission cross section [i][1], m^2
    s_sig: Vec<Vec<f64>>,
    // Signals gamma, adimensional
    s_gam: Vec<f64>,
    //Signal attenuation coefficient to account for fiber losses
    s_alp: Vec<f64>,
    // Forward Pump frequency, 1/s
    fp_freq: Vec<f64>,
    // Forward Pump wavelength, m
    fp_wavelength: Vec<f64>,
    // Forward Pump absorption cross section [i][0] and emission cross section [i][1], m^2
    fp_sig: Vec<Vec<f64>>,
    // Forward Pump gamma, adimensional
    fp_gam: Vec<f64>,
    //Forward pump attenuation coefficient to account for fiber losses
    fp_alp: Vec<f64>,
    // Backwards pump frequency, 1/s
    bp_freq: Vec<f64>,
    // Backwards Pump wavelength, m
    bp_wavelength: Vec<f64>,
    // Backwards Pump absorption cross section [i][0] and emission cross section [i][1], m^2
    bp_sig: Vec<Vec<f64>>,
    // Backwards Pump gamma, adimensional
    bp_gam: Vec<f64>,
    //Backwards pump attenuation coefficient to account for fiber losses, Adimensional
    bp_alp: Vec<f64>,
    // Waveguide width, m
    width: f64,
    // Waveguide depth, m
    depth: f64,
    // Waveguide length, m
    length: f64,
    // Area, m^2
    area: f64,
    // Erbium ion density, ions/m^3
    N: f64,
    // Lifetimes in s
    tau: Vec<f64>,
    // Forward Pump power at length 0, Watts
    fp0: Vec<f64>,
    // Backward Pump power at length L, Watts
    bp0: Vec<f64>,
    // Signal power at length 0, Watts
    s0: Vec<f64>,
}


struct ASE_signals {
    // Starting wavelength of backwards ASE in meters
    fASEstart_wavelength: f64,
    // Channel spacing of both forward and backwards ASE, 125 GHz = 1 nm,  1/s
    channel_spacing: f64,
    // Attenuation constant of the ASE, depends of fibre, adimensional
    fASE_alp: Vec<f64>,
    // Forward ASE cross section, [i][0] absorption, [i][1] emission,  m^2
    fASE_sig: Vec<Vec<f64>>,
    // Wavelenght of forward ASE signals, m
    fASE_wavelength: Vec<f64>,
    // Frequency of forward ASE signals, 1/s
    fASE_freq: Vec<f64>,
    // Overlap factor at signal wavelength, Adimensional
    fASE_gam: Vec<f64>,
    // Forward ASE starting power z = 0
    fASE0: Vec<f64>,
    // Starting wavelength of backwards ASE in meters
    bASEstart_wavelength: f64,
    // Attenuation constant of the ASE, depends of fibre, adimensional
    bASE_alp: Vec<f64>,
    // Backwards ASE cross section, [i][0] absorption, [i][1] emission,  m^2
    bASE_sig: Vec<Vec<f64>>,
    // Wavelenght of backwards ASE signals, m
    bASE_wavelength: Vec<f64>,
    // Frequency of backwards ASE signals, 1/s
    bASE_freq: Vec<f64>,
    // Overlap factor at signal wavelength, Adimensional
    bASE_gam: Vec<f64>,
    // Backwards ASE starting power z = l
    bASE0: Vec<f64>,
}

#[derive(Clone)]
struct constants {
    //A constant in the N2 equation for each signal
    A: Vec<f64>,
    //D constant in the N2 equation for each signal
    D: Vec<f64>,
    //H constant in the dPs equation for each signal
    H: Vec<f64>,
    // Forward Pump constant in the N2 equation
    fC: Vec<f64>,
    // Forward Pump constant in the N2 equation
    fF: Vec<f64>,
    // Forward Pump constant in the dPp equations
    fG: Vec<f64>,
    // Backwards Pump constant in the N2 equation
    bC: Vec<f64>,
    // Backwards Pump constant in the N2 equation
    bF: Vec<f64>,
    // Backwards Pump constant in the dPp equations
    bG: Vec<f64>,
    // B constant for N2 functions = tau_21 * sigma_a * gamma_ASE/ A * h * frequency of signal
    fB: Vec<f64>,
    // E constant for N2 function = tau_21 * (sigma_a + sigma_e) * gamma_ASE/ A * h * frequency of signal
    fE: Vec<f64>,
    // I constant for dPASE function = h * freq * channel spacing * overlap
    fI: Vec<f64>,
    // J constant for dPASE function =  sigma_e + sigma_a for each signal
    fJ: Vec<f64>,
    // B constant for N2 functions = tau_21 * sigma_a * gamma_ASE/ A * h * frequency of signal
    bB: Vec<f64>,
    // E constant for N2 function = tau_21 * (sigma_a + sigma_e) * gamma_ASE/ A * h * frequency of signal
    bE: Vec<f64>,
    // I constant for dPASE function = h * freq * channel spacing * overlap
    bI: Vec<f64>,
    // J constant for dPASE function =  sigma_e + sigma_a for each signal
    bJ: Vec<f64>,
}

#[derive(Clone)]
struct vectors {
    //Contains the evolution of each individual signal per subdivision
    signals:Vec<Vec<f64>>,
    //Contains the evolution of each individual forward pump signal per subdivision
    fpumps: Vec<Vec<f64>>,
    //Contains the evolution of each individual backwards pump signal per subdivision
    bpumps: Vec<Vec<f64>>,
    //Contains the evolution of each individual forward ASE signal per subdivision
    fASE: Vec<Vec<f64>>,
    //Contains the evolution of each individual backwards ASE signal per subdivision
    bASE: Vec<Vec<f64>>,
    //Contains the sums of each signal,
    //[i][0] = N2, [i][1] = signal, [i][2] = fpumps, [i][3] = bpumps, [i][4] = fASE, [i][5] = bASE
    sums: Vec<Vec<f64>>,
    //Z position in the waveguide
    position: Vec<f64>,
    //Evolution of the signals from z = l to z = 0
    rsignals:Vec<Vec<f64>>,
    //Evolution of forward pump from z = l to z = 0
    rfpumps: Vec<Vec<f64>>,
    //Evolution of backwards pump from z = l to z = 0
    rbpumps: Vec<Vec<f64>>,
    //Evolution of forwards ASE from z = l to z = 0
    rfASE: Vec<Vec<f64>>,
    //Evolution of backwards ASE from z = l to z = 0
    rbASE: Vec<Vec<f64>>,
    //Evolution of sums from z = l to z = 0
    rsums: Vec<Vec<f64>>,
}

//Runge Kutta coefficients
struct rk {
    c1: f64,
    c2: f64,
    c3: f64,
    c4: f64,
}

//////////////////////////FUNCTIONS FOR MODE CALCULATION////////////////////////

fn find_ang_freq(frequency: f64) -> f64 {
    2.0 * PI * frequency
}

//Unused function
/*
fn is_unique(k_x1: f64, k_y1: f64, gamma_x1: f64, gamma_y1: f64,
    k_x2: f64, k_y2: f64, gamma_x2: f64, gamma_y2: f64) -> bool {
    let mut comp1: bool = false;
    let mut comp2: bool = false;
    let mut comp3: bool = false;
    let mut comp4: bool = false;

    if k_x1 > 0.0 && k_x2 > 0.0 {
        comp1 = ((k_x1 - k_x2)* 100.0/k_x1).abs() > 1.0
    } else {comp1 = true}

    if k_y1 > 0.0 && k_y2 > 0.0 {
        comp2 = ((k_y1 - k_y2)* 100.0/k_y1).abs() > 1.0
    } else {comp2 = true}

    if gamma_x1 > 0.0 && gamma_x2 > 0.0 {
        comp3 = ((gamma_x1 - gamma_x2)* 100.0/gamma_x1).abs() > 1.0
    } else {comp3 = true}

    if gamma_y1 > 0.0 && gamma_y2 > 0.0 {
        comp4 = ((gamma_y1 - gamma_y2)* 100.0/gamma_y1).abs() > 1.0
    } else {comp4 = true}

    comp1 == comp2 && comp3 == comp4

}
*/

////////////////////FUNCTIONS FOR GAIN CALCULATION//////////////////////////////

//This function takes wavelength in m and outputs frequency in 1/s
fn find_freq(wavelength: f64) -> f64 {
        c / wavelength
}

//Inputs width and depth in m and outputs area in m^2
fn find_area(width: f64, depth: f64) -> f64 {
    width * depth
}

//Calculates gain in dB based on output power vs input power in Watts (or other
//units)
fn find_gain(output: f64, input: f64) -> f64 {
    10.0 * (output / input).log(10.0)
}

//Calculates absorption cross section for Ge-Al doped silica
fn find_absorption(wavelen: f64) -> f64 {
    let wavelength: f64 = wavelen * 1.0e9;
    if wavelength >= 975.0 && wavelength <= 982.0 {
        1.70e-25
    } else if wavelength >= 1449.0 && wavelength < 1530.0 {
        (1.63944165217437E-33) * (wavelength).powi(5) +
        (-1.22017921379809E-29) * (wavelength).powi(4) +
        (3.63241277092243E-26) * (wavelength).powi(3) +
        (-5.40655205496438E-23) * (wavelength).powi(2) +
        (4.02345299618307E-20) * (wavelength) +
        (-1.19762418396094E-17)

    } else if wavelength >= 1530.0 && wavelength < 1600.0 {
        (-5.11559925330538E-33) * (wavelength).powi(5) +
        (4.01173413861064E-29) * (wavelength).powi(4) +
        (-1.25835271657889E-25) * (wavelength).powi(3) +
        (1.97341030852107E-22) * (wavelength).powi(2) +
        (-1.54730983982823E-19) * (wavelength) +
        (4.85257456503714E-17)
    } else {-1.0}

}

//Calculates emission cross section for Ge-Al doped silica
fn find_emission(value: f64) -> f64 {
    let wavelength: f64 = value * 1.0e9;
    if wavelength >= 975.0 && wavelength < 981.0 {
        0.0
    } else if wavelength >= 1449.0 && wavelength < 1530.0 {
        (4.00655760004503E-34) * (wavelength).powi(5) +
        (-2.93490351214106E-30) * (wavelength).powi(4) +
        (8.59912422523874E-27) * (wavelength).powi(3) +
        (-1.25968467454396E-23) * (wavelength).powi(2) +
        (9.22605495880702E-21) * (wavelength) +
        (-2.7027520389058E-18)

    } else if wavelength >= 1530.0 && wavelength < 1600.0 {
        (-7.55679820575314E-33) * (wavelength).powi(5) +
        (5.92783629831709E-29) * (wavelength).powi(4) +
        (-1.85987343916262E-25) * (wavelength).powi(3) +
        (2.91749017649647E-22) * (wavelength).powi(2) +
        (-2.28809598569961E-19) * (wavelength) +
        (7.1774126020988E-17)
    } else {-1.0}
}

//Calculates inverted populations
fn N2 (N: f64, APS: f64, BPA: f64, CPP: f64, DPS: f64, EPA: f64, FPP: f64)
-> f64 {
    N * (APS + BPA + CPP) / (DPS + EPA + FPP + 1.0)
}

fn dPs (N: f64, H: f64, sigma_signal: f64, gamma_signal: f64, alpha_signal: f64,
N2: f64, Ps: f64) -> f64 {
     (N2 * H - N * sigma_signal) * gamma_signal * Ps - alpha_signal * Ps
}

fn dPp (N: f64, G: f64, sigma_pump: f64, gamma_pump: f64, alpha_pump: f64,
N2: f64, Pp: f64) -> f64 {
     (N2 * G - N * sigma_pump) * gamma_pump * Pp - alpha_pump * Pp
}

fn dPASE (N: f64, J: f64, I: f64, sigma_ASE: f64, gamma_ASE: f64, alpha_ASE: f64,
N2: f64, PASE: f64) -> f64 {
     (N2 * J - N * sigma_ASE) * gamma_ASE * PASE + N2 * I - alpha_ASE * PASE
}

//////////////////////////MODS FOR MODE CALCULATION/////////////////////////////
mod E_x {

    pub fn gamma_x(n_0: f64, n_1: f64, k: f64, k_x: f64) -> f64 {
        ((n_1.powi(2) - n_0.powi(2)) * k.powi(2) - k_x.powi(2)).sqrt()
    }

    pub fn gamma_y(n_0: f64, n_1: f64, k: f64, k_y: f64) -> f64 {
        ((n_1.powi(2)-n_0.powi(2)) * k.powi(2) - k_y.powi(2)).sqrt()
    }

    pub fn find_k_x(n_0: f64, n_1: f64, p: f64, a: f64, PI: f64, gamma_x: f64, k_x: f64)
    -> f64 {
        (k_x * a - ((p - 1.0) * 0.5 * PI +
        ((n_1.powi(2) * gamma_x) / (n_0.powi(2) * k_x)).atan())).abs()
    }

    pub fn find_k_y(n_0: f64, n_1: f64, q: f64, d: f64, PI: f64, gamma_y: f64, k_y: f64)
    -> f64 {
        (k_y * d - ((q - 1.0) * 0.5 * PI + ((gamma_y) / (k_y)).atan())).abs()
    }

    pub fn find_beta(n_0: f64, n_1: f64, k: f64, k_x: f64, k_y: f64) -> f64 {
        (k.powi(2) * n_1.powi(2) - (k_x.powi(2) + k_y.powi(2))).sqrt()
    }

    pub fn find_n_eff(beta: f64, k: f64) -> f64 {
        beta / k
    }

    pub fn phi(p: f64, PI: f64) -> f64 {
        0.5 * (p - 1.0) * PI
    }

    pub fn psi(q: f64, PI: f64) -> f64 {
        0.5 * (q - 1.0) * PI
    }

    pub fn Hy1(A: f64, k_x: f64, k_y: f64, phi: f64, psi: f64, x: f64, y: f64)
    -> f64 {
        A * (k_x * x - phi).cos() * (k_y * y - psi).cos()
    }

    pub fn Hy2(A: f64, k_x: f64, k_y: f64, gamma_x: f64, phi: f64, psi: f64,
        x: f64, y: f64, a: f64) -> f64 {
        A * (k_x * a - phi).cos() * (-1.0 * gamma_x * (x - a)).exp()
        * (k_y * y - psi).cos()
    }

    pub fn Hy3(A: f64, k_x: f64, k_y: f64, gamma_y: f64, phi: f64, psi: f64,
        x: f64, y: f64, d: f64) -> f64 {
        A * (k_x * x - phi).cos() * (-1.0 * gamma_y * (y - d)).exp()
        * (k_y * d - psi).cos()
    }
}

mod E_y {

    pub fn gamma_x(n_0: f64, n_1: f64, k: f64, k_x: f64) -> f64 {
        ((n_1.powi(2) - n_0.powi(2)) * k.powi(2) - k_x.powi(2)).sqrt()
    }

    pub fn gamma_y(n_0: f64, n_1: f64, k: f64, k_y: f64) -> f64 {
        ((n_1.powi(2)-n_0.powi(2)) * k.powi(2) - k_y.powi(2)).sqrt()
    }

    pub fn find_k_x(n_0: f64, n_1: f64, p: f64, a: f64, PI: f64, gamma_x: f64, k_x: f64)
    -> f64 {
        (k_x * a - ((p - 1.0) * 0.5 * PI +
        (gamma_x / k_x).atan())).abs()
    }

    pub fn find_k_y(n_0: f64, n_1: f64, q: f64, d: f64, PI: f64, gamma_y: f64, k_y: f64)
    -> f64 {
        (k_y * d - ((q - 1.0) * 0.5 * PI +
        ((n_1.powi(2) * (gamma_y)) / (n_0.powi(2) * (k_y))).atan())).abs()
    }

    pub fn find_beta(n_0: f64, n_1: f64, k: f64, k_x: f64, k_y: f64) -> f64 {
        (k.powi(2) * n_1.powi(2) - (k_x.powi(2) + k_y.powi(2))).sqrt()
    }

    pub fn find_n_eff(beta: f64, k: f64) -> f64 {
        beta / k
    }

    pub fn phi(p: f64, PI: f64) -> f64 {
        0.5 * (p - 1.0) * PI
    }

    pub fn psi(q: f64, PI: f64) -> f64 {
        0.5 * (q - 1.0) * PI
    }

    pub fn Hx1(A: f64, k_x: f64, k_y: f64, phi: f64, psi: f64, x: f64, y: f64)
    -> f64 {
        A * (k_x * x - phi).cos() * (k_y * y - psi).cos()
    }

    pub fn Hx2(A: f64, k_x: f64, k_y: f64, gamma_x: f64, phi: f64, psi: f64,
        x: f64, y: f64, a: f64) -> f64 {
        A * (k_x * a - phi).cos() * (-1.0 * gamma_x * (x - a)).exp()
        * (k_y * y - psi).cos()
    }

    pub fn Hx3(A: f64, k_x: f64, k_y: f64, gamma_y: f64, phi: f64, psi: f64,
        x: f64, y: f64, d: f64) -> f64 {
        A * (k_x * x - phi).cos() * (-1.0 * gamma_y * (y - d)).exp()
        * (k_y * d - psi).cos()
    }
}


fn main() {
    println!("Starting modes calculation");
    println!("________________________________________________");
    /////////////////////////////////Reading config.txt/////////////////////////
    let path = Path::new("config.txt");
    let display = path.display();
    println!("display = {}", display);
    let mut file = match File::open(path) {
        Err(why) => panic!("couldn't open {:?}, {}", path, why.description()),
        Ok(file) => file,
    };

    let mut contents: String = String::new();
    match file.read_to_string(&mut contents) {
        Ok(_) => println!("File read successfully!"),
        Err(why) => panic!("Couldn't read file {:?}, {}", display, why.description()),
    }

    let val: Vec<&str> = contents.split("\n").collect();

    let mut val2: Vec<&str> = vec![];
    let mut val3: Vec<f64> = vec![];
    let mut val4: Vec<bool> = vec![];
    let mut inter: Vec<&str> = vec![];
    for elements in val.into_iter() {
        inter = elements.split("=").collect();
        for a in inter.iter() {
            val2.push(a.trim());
        }
    }

    for elements in val2.iter() {
        match FromStr::from_str(elements) {
            Ok(value) => val4.push(value),
            Err(_) => continue,
            }
    }

    for elements in val2.into_iter() {
        match elements.parse::<f64>() {
            Ok(var) => val3.push(var),
            Err(_) => continue,
        }
    }

    let STEPS: usize = val3[16] as usize;
    let mut values_to_iterate: Vec<f64> = vec![];
    values_to_iterate.push(val3[1] * 1e-9);                                            //Forward pump wavelength
    values_to_iterate.push(val3[2] * 1e-9);                                            //Backwards pump wavelength
    for g in 0..SIGNALS {
        let value: f64 = val3[0] * 1e-9 + (g as f64) * 1e-9;
        values_to_iterate.push(value);
    }
    let size_of_values: Box<i32> = Box::new(values_to_iterate.len() as i32);

    if val4[1] == true
//start of mode block
{
    let mut result_file = match File::create("mode_results.csv") {
        Ok(content) => content,
        Err(why)    => panic!("Couldn't create file, {}", why.description()),
    };
    let mut completion: i32 = 0;
    if val4[0] == false {
    for g in values_to_iterate.into_iter() {
    println!("{} signal out of {:?} completed", completion + 1,size_of_values);
    println!("Calculating modes for {} nm", (g * 1e9).round());
    println!("__________________________________");
    /////////////////////////////////PARAMETER SET UP///////////////////////////
    let mut params: mode_params = mode_params {
        n_0: val3[14],                                                          // Refractive index of the cladding
        n_1: val3[15],                                                          // Refractive index of the core
        wavelength: g,                                                          // Wavelength in m
        frequency: 0.0,                                                         // Frequency in Hz
        w_ang: 0.0,                                                             // Angular frequency in rad/s
        k: 0.0,                                                                 // Wavenumber in Vacuum
        a: val3[4] * 1.0e-6 * 0.5,                                              // width
        d: val3[3] * 1.0e-6 * 0.5,                                              // depth
    };

    let mut vectors: mode_vectors = mode_vectors {
        k_x: [0.0; P],
        k_y: [0.0; Q],
        gamma_x: [0.0; P],
        gamma_y: [0.0; Q],
        beta:[[0.0; P]; Q],
        sols_Ex: [[0.0; 13]; P*Q], //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
        sols_Ey: [[0.0; 13]; P*Q], //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
    };

    params.frequency = find_freq(params.wavelength);
    params.w_ang = find_ang_freq(params.frequency);
    params.k = params.w_ang / c;

    let mut a: f64 = params.a;
    let mut d: f64 = params.d;
    let mut k: f64 = params.k;
    let mut n_0: f64 = params.n_0;
    let mut n_1: f64 = params.n_1;
    let mut stepsize: f64 = k * n_1 / (STEPSMODE as f64);
    let mut p: f64 = 0.0;
    let mut q: f64 = 0.0;
    let mut beta: f64 = 0.0;

    let mut gamma_x: f64 = 0.0;
    let mut k_x: f64 = 0.0;
    let mut gamma_y: f64 = 0.0;
    let mut k_y: f64= 0.0;
    let mut n_eff: f64 = 0.0;
    let mut comparator: f64 = 10.0;
    let mut index: usize = 0;
    let mut newlow: f64 = 1000000.0;

    ////////////////////////////////CALCULATION/////////////////////////////////
    //////////////////////////////////E_X mode/////////////////////////////////
    for w in 0..P {
        newlow = 1000000.0;
        for i in 0..STEPSMODE {
            /*gamma_x, k_x calculation*/
            k_x = (i as f64) * stepsize;
            gamma_x = E_x::gamma_x(n_0, n_1, k, k_x);
            p = w as f64 + 1.0;
            comparator = E_x::find_k_x(n_0, n_1, p, a, PI, gamma_x, k_x);
            if comparator < 1e-3 && comparator < newlow &&
            k_x > 10.0 {
                vectors.k_x[w] = k_x;
                vectors.gamma_x[w] = gamma_x;
                newlow = comparator;
            }
        }
    }

    for w in 0..Q {
        newlow = 1000000.0;
        for i in 0..STEPSMODE {
            /*gamma_y, k_y calculation*/
            k_y = (i as f64) * stepsize;
            gamma_y = E_x::gamma_y(n_0, n_1, k, k_y);
            q = w as f64 + 1.0;
            comparator = E_x::find_k_y(n_0, n_1, q, d, PI, gamma_y, k_y);
            if comparator < 1e-3 && comparator < newlow &&
            k_y > 10.0 {
                vectors.k_y[w] = k_y;
                vectors.gamma_y[w] = gamma_y;
                newlow = comparator;
            }
        }
    }

    'first:for w in 0..P {
        'second:for e in 0..Q{
            k_x = vectors.k_x[w];
            k_y = vectors.k_y[e];
            gamma_x = vectors.gamma_x[w];
            gamma_y = vectors.gamma_y[e];
            beta = E_x::find_beta(n_0, n_1, k, k_x, k_y);
            n_eff = E_x::find_n_eff(beta, k);
            if (n_eff > n_0 && n_eff < n_1) && (k_x > 1.0 && k_y > 1.0) {
                vectors.beta[w][e] = beta;
                vectors.sols_Ex[index] = [(w+1) as f64, (e+1) as f64, a, d, n_0,
                    n_1, k_x, k_y, gamma_x, gamma_y, beta, n_eff, 0.0];
                index += 1;
                //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
            }
        }
    }

    index = 0;
{
    let mut counter: i32 = 0;

    for i in vectors.sols_Ex.iter() {
        if i[0] != 0.0 {
            counter += 1
        }

    }
    if counter > 0 {
        for t in 0..vectors.sols_Ex.len() {
            if vectors.sols_Ex[t][6] > 0.5 {
                let c1: f64 = (params.w_ang * mu_0 * (vectors.sols_Ex[t][10].powi(2))
                - (vectors.sols_Ex[t][6].powi(2)))/(vectors.sols_Ex[t][10].powi(3));
                let c2: f64 = (params.w_ang * mu_0 * (vectors.sols_Ex[t][10].powi(2))
                + (vectors.sols_Ex[t][8].powi(2)))/(vectors.sols_Ex[t][10].powi(3));
                let c3: f64 = c1.clone();
                let mut sum_x: f64 = 0.0;
                let mut sum_y: f64 = 0.0;
                let mut x: f64 = 0.0;
                let mut y: f64 = 0.0;
            //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
                let a: f64 = vectors.sols_Ex[t][2];
                let d: f64 = vectors.sols_Ex[t][3];
                let k_x: f64 = vectors.sols_Ex[t][6];
                let k_y: f64 = vectors.sols_Ex[t][7];
                let gamma_x: f64 = vectors.sols_Ex[t][8];
                let gamma_y: f64 = vectors.sols_Ex[t][9];
                let beta: f64 = vectors.sols_Ex[t][11];
                let phi: f64 = E_x::phi(vectors.sols_Ex[t][0], PI);
                let psi: f64 = E_x::psi(vectors.sols_Ex[t][1], PI);

                let range_x  = Range::new(-1.0 * a, a);
                let range_y  = Range::new(-1.0 * d, d);
                let mut rng = rand::thread_rng();
                let N: i32 = 1000;
                let cycles: i32 = 1;
                //Region 1, -a<=x<=a, -d<=y<=d
                for z in 0..cycles {
                    for i in 0..N {
                        y = range_y.ind_sample(&mut rng);
                        for j in 0..N {
                            x = range_x.ind_sample(&mut rng);
                            sum_x += (E_x::Hy1(1.0, k_x, k_y, phi, psi, x, y)).powi(2)

                        }
                        sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                        sum_x = 0.0;
                    }
                }
                let P1: f64 = c1 * (2.0 * d) * sum_y/((N as f64) * (cycles as f64));
                sum_x = 0.0;
                sum_y = 0.0;
                let range_x  = Range::new(a, 3.0 * a);
                let range_y  = Range::new(-1.0 * d, d);
                let mut rng = rand::thread_rng();
                //Region 2, a<=x<=3a, -d<=y<=d
                for z in 0..cycles {
                    for i in 0..N {
                        y = range_y.ind_sample(&mut rng);
                        for j in 0..N {
                            x = range_x.ind_sample(&mut rng);
                            sum_x += (E_x::Hy2(1.0, k_x, k_y, gamma_x, phi, psi, x, y, a)).powi(2)

                        }
                        sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                        sum_x = 0.0;
                    }
                }

                let P2: f64 = c2 * (2.0 * d) * sum_y / ((N as f64) * (cycles as f64));
                let range_x  = Range::new(-1.0 * a, a);
                let range_y  = Range::new(d, 3.0 * d);
                let mut rng = rand::thread_rng();
                sum_x = 0.0;
                sum_y = 0.0;
                //Region 3, -a<=x<=a, d<=y<=3d
                for z in 0..cycles {
                    for i in 0..N {
                        y = range_y.ind_sample(&mut rng);
                        for j in 0..N {
                            x = range_x.ind_sample(&mut rng);
                            sum_x += (E_x::Hy3(1.0, k_x, k_y, gamma_y, phi, psi, x, y, d)).powi(2)

                        }
                        sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                        sum_x = 0.0;
                    }
                }
                let P3: f64 = c3 * (2.0 * d) * sum_y / ((N as f64) * (cycles as f64));
                let containment = P1 / (P1+P2+P3);
                vectors.sols_Ex[t][12] = containment;
            }
        }
    /*
    println!("vectors.sols_Ex[0] = {:?}", vectors.sols_Ex);
    println!("0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x");
    */

    }

    write!(result_file, "E_x mode,{},nm\n", (g * 1e9).round());
    write!(result_file,
        "P, Q, a/2, d/2, n_0, n_1, kx, ky, gamma_x, gamma_y, beta, n_eff, cont_x\n");
    for elements in vectors.sols_Ex.iter() {
        if elements[0] > 0.0 {
            write!(result_file, "{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
                elements[0], elements[1], elements[2], elements[3], elements[4],
                elements[5], elements[6], elements[7], elements[8], elements[9],
                elements[10], elements[11], elements[12]
            );

        }
    }

//end of block
//0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
}

////////////////////////////////RESET////////////////////////////////////////
vectors.k_x = [0.0; P];
vectors.k_y = [0.0; Q];
vectors.gamma_x = [0.0; P];
vectors.gamma_y = [0.0; Q];
vectors.beta = [[0.0; P]; Q];

a = params.a;
d = params.d;
k = params.k;
n_0 = params.n_0;
n_1 = params.n_1;
stepsize = k * n_1 / (STEPSMODE as f64);
p = 0.0;
q = 0.0;
beta = 0.0;
gamma_x = 0.0;
k_x = 0.0;
gamma_y = 0.0;
k_y = 0.0;
n_eff = 0.0;
comparator = 10.0;
index = 0;
newlow = 1000000.0;

////////////////////////////////CALCULATION/////////////////////////////////
//////////////////////////////////E_y mode/////////////////////////////////
for w in 0..P {
    newlow = 1000000.0;
    for i in 0..STEPSMODE {
        /*gamma_x, k_x calculation*/
        k_x = (i as f64) * stepsize;
        gamma_x = E_y::gamma_x(n_0, n_1, k, k_x);
        p = w as f64 + 1.0;
        comparator = E_y::find_k_x(n_0, n_1, p, a, PI, gamma_x, k_x);
        if comparator < 1e-3 && comparator < newlow &&
        k_x > 10.0 {
            vectors.k_x[w] = k_x;
            vectors.gamma_x[w] = gamma_x;
            newlow = comparator;
        }
    }
}

for w in 0..Q {
    newlow = 1000000.0;
    for i in 0..STEPSMODE {
        /*gamma_y, k_y calculation*/
        k_y = (i as f64) * stepsize;
        gamma_y = E_y::gamma_y(n_0, n_1, k, k_y);
        q = w as f64 + 1.0;
        comparator = E_y::find_k_y(n_0, n_1, q, d, PI, gamma_y, k_y);
        if comparator < 1e-3 && comparator < newlow &&
        k_y > 10.0 {
            vectors.k_y[w] = k_y;
            vectors.gamma_y[w] = gamma_y;
            newlow = comparator;
        }
    }
}

'first:for w in 0..P {
    'second:for e in 0..Q{
        k_x = vectors.k_x[w];
        k_y = vectors.k_y[e];
        gamma_x = vectors.gamma_x[w];
        gamma_y = vectors.gamma_y[e];
        beta = E_y::find_beta(n_0, n_1, k, k_x, k_y);
        n_eff = E_y::find_n_eff(beta, k);
        if (n_eff > n_0 && n_eff < n_1) && (k_x > 1.0 && k_y > 1.0) {
            vectors.beta[w][e] = beta;
            vectors.sols_Ey[index] = [(w+1) as f64, (e+1) as f64, a, d, n_0,
                n_1, k_x, k_y, gamma_x, gamma_y, beta, n_eff, 0.0];
            index += 1;
            //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
        }
    }
}

index = 0;
//start of block
{
let mut counter: i32 = 0;
for i in vectors.sols_Ey.iter() {
    if i[0] != 0.0 {
        counter += 1
    }

}
if counter > 0 {
    for t in 0..vectors.sols_Ey.len() {
        if vectors.sols_Ey[t][6] > 0.5 {
            let c1: f64 = (params.w_ang * mu_0 * (vectors.sols_Ey[t][10].powi(2))
            - (vectors.sols_Ey[t][7].powi(2)))/(vectors.sols_Ey[t][10].powi(3));
            let c2: f64 = c1.clone();
            let c3: f64 = (params.w_ang * mu_0 * (vectors.sols_Ey[t][10].powi(2))
            + (vectors.sols_Ey[t][9].powi(2)))/(vectors.sols_Ey[t][10].powi(3));


            let mut sum_x: f64 = 0.0;
            let mut sum_y: f64 = 0.0;
            let mut x: f64 = 0.0;
            let mut y: f64 = 0.0;
        //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
            let a: f64 = vectors.sols_Ey[t][2];
            let d: f64 = vectors.sols_Ey[t][3];
            let k_x: f64 = vectors.sols_Ey[t][6];
            let k_y: f64 = vectors.sols_Ey[t][7];
            let gamma_x: f64 = vectors.sols_Ey[t][8];
            let gamma_y: f64 = vectors.sols_Ey[t][9];
            let beta: f64 = vectors.sols_Ey[t][11];
            let phi: f64 = E_y::phi(vectors.sols_Ey[t][0], PI);
            let psi: f64 = E_y::psi(vectors.sols_Ey[t][1], PI);

            let range_x  = Range::new(-1.0 * a, a);
            let range_y  = Range::new(-1.0 * d, d);
            let mut rng = rand::thread_rng();
            let N: i32 = 1000;
            let cycles: i32 = 1;
            //Region 1, -a<=x<=a, -d<=y<=d
            for z in 0..cycles {
                for i in 0..N {
                    y = range_y.ind_sample(&mut rng);
                    for j in 0..N {
                        x = range_x.ind_sample(&mut rng);
                        sum_x += (E_y::Hx1(1.0, k_x, k_y, phi, psi, x, y)).powi(2)

                    }
                    sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                    sum_x = 0.0;
                }
            }
            let P1: f64 = c1 * (2.0 * d) * sum_y / ((N as f64) * (cycles as f64));
            sum_x = 0.0;
            sum_y = 0.0;
            let range_x  = Range::new(a, 3.0 * a);
            let range_y  = Range::new(-1.0 * d, d);
            let mut rng = rand::thread_rng();
            //Region 2, a<=x<=3a, -d<=y<=d
            for z in 0..cycles {
                for i in 0..N {
                    y = range_y.ind_sample(&mut rng);
                    for j in 0..N {
                        x = range_x.ind_sample(&mut rng);
                        sum_x += (E_y::Hx2(1.0, k_x, k_y, gamma_x, phi, psi, x, y, a)).powi(2)

                    }
                    sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                    sum_x = 0.0;
                }
            }

            let P2: f64 = c2 * (2.0 * d) * sum_y / ((N as f64) * (cycles as f64));
            let range_x  = Range::new(-1.0 * a, a);
            let range_y  = Range::new(d, 3.0 * d);
            let mut rng = rand::thread_rng();
            sum_x = 0.0;
            sum_y = 0.0;
            //Region 3, -a<=x<=a, d<=y<=3d
            for z in 0..cycles {
                for i in 0..N {
                    y = range_y.ind_sample(&mut rng);
                    for j in 0..N {
                        x = range_x.ind_sample(&mut rng);
                        sum_x += (E_y::Hx3(1.0, k_x, k_y, gamma_y, phi, psi, x, y, d)).powi(2)

                    }
                    sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                    sum_x = 0.0;
                }
            }
            let P3: f64 = c3 * (2.0 * d) * sum_y / ((N as f64) * (cycles as f64));
            let containment = P1 / (P1+P2+P3);
            vectors.sols_Ey[t][12] = containment;
        }
    }
/*
println!("vectors.sols_Ex[0] = {:?}", vectors.sols_Ex);
println!("0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x");
*/

}
write!(result_file, "E_y mode,{},nm\n", (g * 1e9).round());
write!(result_file,
    "P, Q, a/2, d/2, n_0, n_1, kx, ky, gamma_x, gamma_y, beta, n_eff, cont_x\n");

for elements in vectors.sols_Ey.iter() {
    if elements[0] > 0.0 {
        write!(result_file, "{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
            elements[0], elements[1], elements[2], elements[3], elements[4],
            elements[5], elements[6], elements[7], elements[8], elements[9],
            elements[10], elements[11], elements[12]
        );
    }
}
//0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
//End of block
}
completion += 1;
//End of for loop values_to_iterate
}
//End of if false block
drop(result_file);
println!("Mode calculation completed")
} else {

let mut min_wavelength: f64 = 200000000.0;
let mut outer_counter: i64 = 1;
for elements in values_to_iterate.iter() {
    if *elements < min_wavelength {
        min_wavelength = elements.clone();
    }
}
let mut result_file = match File::create("mode_results.csv") {
    Ok(content) => content,
    Err(why)    => panic!("Couldn't create file, {}", why.description()),
};
let mut completion: i32 = 0;
let mut n_x: f64 = 0.0;
let mut n_y: f64 = 0.0;
let mut end_cond: [i32; 2] = [0;2];

'outer:loop {
/////////////////////////////////PARAMETER SET UP///////////////////////////

let mut params: mode_params = mode_params {
    n_0: val3[14],                                                          // Refractive index of the cladding
    n_1: val3[14] + ((outer_counter as f64)/10000.0),                                  // Refractive index of the core
    wavelength: min_wavelength,                                              // Wavelength in m
    frequency: 0.0,                                                         // Frequency in Hz
    w_ang: 0.0,                                                             // Angular frequency in rad/s
    k: 0.0,                                                                 // Wavenumber in Vacuum
    a: val3[4] * 1.0e-6 * 0.5,                                                  // width
    d: val3[3] * 1.0e-6 * 0.5,                                                  // depth
};

let mut vectors: mode_vectors = mode_vectors {
    k_x: [0.0; P],
    k_y: [0.0; Q],
    gamma_x: [0.0; P],
    gamma_y: [0.0; Q],
    beta:[[0.0; P]; Q],
    sols_Ex: [[0.0; 13]; P*Q], //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
    sols_Ey: [[0.0; 13]; P*Q], //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
};

params.frequency = find_freq(params.wavelength);
params.w_ang = find_ang_freq(params.frequency);
params.k = params.w_ang / c;

let mut a: f64 = params.a;
let mut d: f64 = params.d;
let mut k: f64 = params.k;
let mut n_0: f64 = params.n_0;
let mut n_1: f64 = params.n_1;
let mut stepsize: f64 = k * n_1 / (STEPSMODE as f64);
let mut p: f64 = 0.0;
let mut q: f64 = 0.0;
let mut beta: f64 = 0.0;

let mut gamma_x: f64 = 0.0;
let mut k_x: f64 = 0.0;
let mut gamma_y: f64 = 0.0;
let mut k_y: f64= 0.0;
let mut n_eff: f64 = 0.0;
let mut comparator: f64 = 10.0;
let mut index: usize = 0;
let mut newlow: f64 = 1000000.0;
////////////////////////////////CALCULATION/////////////////////////////////
//////////////////////////////////E_X mode/////////////////////////////////
for w in 0..P {
    newlow = 1000000.0;
    for i in 0..STEPSMODE {
        /*gamma_x, k_x calculation*/
        k_x = (i as f64) * stepsize;
        gamma_x = E_x::gamma_x(n_0, n_1, k, k_x);
        p = w as f64 + 1.0;
        comparator = E_x::find_k_x(n_0, n_1, p, a, PI, gamma_x, k_x);
        if comparator < 1e-3 && comparator < newlow &&
        k_x > 10.0 {
            vectors.k_x[w] = k_x;
            vectors.gamma_x[w] = gamma_x;
            newlow = comparator;
        }
    }
}

for w in 0..Q {
    newlow = 1000000.0;
    for i in 0..STEPSMODE {
        /*gamma_y, k_y calculation*/
        k_y = (i as f64) * stepsize;
        gamma_y = E_x::gamma_y(n_0, n_1, k, k_y);
        q = w as f64 + 1.0;
        comparator = E_x::find_k_y(n_0, n_1, q, d, PI, gamma_y, k_y);
        if comparator < 1e-3 && comparator < newlow &&
        k_y > 10.0 {
            vectors.k_y[w] = k_y;
            vectors.gamma_y[w] = gamma_y;
            newlow = comparator;
        }
    }
}

'first:for w in 0..P {
    'second:for e in 0..Q{
        k_x = vectors.k_x[w];
        k_y = vectors.k_y[e];
        gamma_x = vectors.gamma_x[w];
        gamma_y = vectors.gamma_y[e];
        beta = E_x::find_beta(n_0, n_1, k, k_x, k_y);
        n_eff = E_x::find_n_eff(beta, k);
        if (n_eff > n_0 && n_eff < n_1) && (k_x > 1.0 && k_y > 1.0) {
            vectors.beta[w][e] = beta;
            vectors.sols_Ex[index] = [(w+1) as f64, (e+1) as f64, a, d, n_0,
                n_1, k_x, k_y, gamma_x, gamma_y, beta, n_eff, 0.0];
            index += 1;
            //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
        }
        if index >= 2 {
            n_x = params.n_1 - ((1.0/10000.0));
            end_cond[0] = 1;
            break 'first
        }
    }
}

////////////////////////////////RESET////////////////////////////////////////
vectors.k_x = [0.0; P];
vectors.k_y = [0.0; Q];
vectors.gamma_x = [0.0; P];
vectors.gamma_y = [0.0; Q];
vectors.beta = [[0.0; P]; Q];

a = params.a;
d = params.d;
k = params.k;
n_0 = params.n_0;
n_1 = params.n_1;
stepsize = k * n_1 / (STEPSMODE as f64);
p = 0.0;
q = 0.0;
beta = 0.0;
gamma_x = 0.0;
k_x = 0.0;
gamma_y = 0.0;
k_y = 0.0;
n_eff = 0.0;
comparator = 10.0;
index = 0;
newlow = 1000000.0;


////////////////////////////////CALCULATION/////////////////////////////////
//////////////////////////////////E_y mode/////////////////////////////////
for w in 0..P {
newlow = 1000000.0;
for i in 0..STEPSMODE {
    /*gamma_x, k_x calculation*/
    k_x = (i as f64) * stepsize;
    gamma_x = E_y::gamma_x(n_0, n_1, k, k_x);
    p = w as f64 + 1.0;
    comparator = E_y::find_k_x(n_0, n_1, p, a, PI, gamma_x, k_x);
    if comparator < 1e-3 && comparator < newlow &&
    k_x > 10.0 {
        vectors.k_x[w] = k_x;
        vectors.gamma_x[w] = gamma_x;
        newlow = comparator;
    }
}
}

for w in 0..Q {
newlow = 1000000.0;
for i in 0..STEPSMODE {
    /*gamma_y, k_y calculation*/
    k_y = (i as f64) * stepsize;
    gamma_y = E_y::gamma_y(n_0, n_1, k, k_y);
    q = w as f64 + 1.0;
    comparator = E_y::find_k_y(n_0, n_1, q, d, PI, gamma_y, k_y);
    if comparator < 1e-3 && comparator < newlow &&
    k_y > 10.0 {
        vectors.k_y[w] = k_y;
        vectors.gamma_y[w] = gamma_y;
        newlow = comparator;
    }
}
}

'first:for w in 0..P {
    'second:for e in 0..Q{
        k_x = vectors.k_x[w];
        k_y = vectors.k_y[e];
        gamma_x = vectors.gamma_x[w];
        gamma_y = vectors.gamma_y[e];
        beta = E_y::find_beta(n_0, n_1, k, k_x, k_y);
        n_eff = E_y::find_n_eff(beta, k);
        if (n_eff > n_0 && n_eff < n_1) && (k_x > 1.0 && k_y > 1.0) {
            vectors.beta[w][e] = beta;
            vectors.sols_Ey[index] = [(w+1) as f64, (e+1) as f64, a, d, n_0,
                n_1, k_x, k_y, gamma_x, gamma_y, beta, n_eff, 0.0];
                index += 1;
                //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
        }
        if index >= 2 {
            n_y = params.n_1 - ((1.0/10000.0));
            end_cond[1] = 1;
            break 'first
        }
    }
}

index = 0;
if end_cond[0] == 1 && end_cond[1] == 1 {
    break 'outer
}
outer_counter += 1;
//end of 'outer loop
println!("testing n_1 = {:.4}", (params.n_1));
}

for g in values_to_iterate.into_iter() {
println!("{} signal out of {:?} completed", completion + 1,size_of_values);
println!("Calculating modes for {} nm", (g * 1e9).round());
println!("__________________________________");
/////////////////////////////////PARAMETER SET UP///////////////////////////
let mut params: mode_params = mode_params {
    n_0: val3[14],                                                          // Refractive index of the cladding
    n_1: {if n_x < n_y {n_y} else {n_x}},                                   // Refractive index of the core
    wavelength: g,                                                          // Wavelength in m
    frequency: 0.0,                                                         // Frequency in Hz
    w_ang: 0.0,                                                             // Angular frequency in rad/s
    k: 0.0,                                                                 // Wavenumber in Vacuum
    a: val3[4] * 1.0e-6 * 0.5,                                                  // width
    d: val3[3] * 1.0e-6 * 0.5,                                                  // depth
};

let mut vectors: mode_vectors = mode_vectors {
    k_x: [0.0; P],
    k_y: [0.0; Q],
    gamma_x: [0.0; P],
    gamma_y: [0.0; Q],
    beta:[[0.0; P]; Q],
    sols_Ex: [[0.0; 13]; P*Q], //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
    sols_Ey: [[0.0; 13]; P*Q], //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
};

params.frequency = find_freq(params.wavelength);
params.w_ang = find_ang_freq(params.frequency);
params.k = params.w_ang / c;

let mut a: f64 = params.a;
let mut d: f64 = params.d;
let mut k: f64 = params.k;
let mut n_0: f64 = params.n_0;
let mut n_1: f64 = params.n_1;
let mut stepsize: f64 = k * n_1 / (STEPSMODE as f64);
let mut p: f64 = 0.0;
let mut q: f64 = 0.0;
let mut beta: f64 = 0.0;

let mut gamma_x: f64 = 0.0;
let mut k_x: f64 = 0.0;
let mut gamma_y: f64 = 0.0;
let mut k_y: f64= 0.0;
let mut n_eff: f64 = 0.0;
let mut comparator: f64 = 10.0;
let mut index: usize = 0;
let mut newlow: f64 = 1000000.0;

////////////////////////////////CALCULATION/////////////////////////////////
//////////////////////////////////E_X mode/////////////////////////////////
for w in 0..P {
    newlow = 1000000.0;
    for i in 0..STEPSMODE {
        /*gamma_x, k_x calculation*/
        k_x = (i as f64) * stepsize;
        gamma_x = E_x::gamma_x(n_0, n_1, k, k_x);
        p = w as f64 + 1.0;
        comparator = E_x::find_k_x(n_0, n_1, p, a, PI, gamma_x, k_x);
        if comparator < 1e-3 && comparator < newlow &&
        k_x > 10.0 {
            vectors.k_x[w] = k_x;
            vectors.gamma_x[w] = gamma_x;
            newlow = comparator;
        }
    }
}

for w in 0..Q {
    newlow = 1000000.0;
    for i in 0..STEPSMODE {
        /*gamma_y, k_y calculation*/
        k_y = (i as f64) * stepsize;
        gamma_y = E_x::gamma_y(n_0, n_1, k, k_y);
        q = w as f64 + 1.0;
        comparator = E_x::find_k_y(n_0, n_1, q, d, PI, gamma_y, k_y);
        if comparator < 1e-3 && comparator < newlow &&
        k_y > 10.0 {
            vectors.k_y[w] = k_y;
            vectors.gamma_y[w] = gamma_y;
            newlow = comparator;
        }
    }
}

'first:for w in 0..P {
    'second:for e in 0..Q{
        k_x = vectors.k_x[w];
        k_y = vectors.k_y[e];
        gamma_x = vectors.gamma_x[w];
        gamma_y = vectors.gamma_y[e];
        beta = E_x::find_beta(n_0, n_1, k, k_x, k_y);
        n_eff = E_x::find_n_eff(beta, k);
        if (n_eff > n_0 && n_eff < n_1) && (k_x > 1.0 && k_y > 1.0) {
            vectors.beta[w][e] = beta;
            vectors.sols_Ex[index] = [(w+1) as f64, (e+1) as f64, a, d, n_0,
                n_1, k_x, k_y, gamma_x, gamma_y, beta, n_eff, 0.0];
            index += 1;
            //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
        }
    }
}

index = 0;
{
let mut counter: i32 = 0;

for i in vectors.sols_Ex.iter() {
    if i[0] != 0.0 {
        counter += 1
    }

}
if counter > 0 {
    for t in 0..vectors.sols_Ex.len() {
        if vectors.sols_Ex[t][6] > 0.5 {
            let c1: f64 = (params.w_ang * mu_0 * (vectors.sols_Ex[t][10].powi(2))
            - (vectors.sols_Ex[t][6].powi(2)))/(vectors.sols_Ex[t][10].powi(3));
            let c2: f64 = (params.w_ang * mu_0 * (vectors.sols_Ex[t][10].powi(2))
            + (vectors.sols_Ex[t][8].powi(2)))/(vectors.sols_Ex[t][10].powi(3));
            let c3: f64 = c1.clone();
            let mut sum_x: f64 = 0.0;
            let mut sum_y: f64 = 0.0;
            let mut x: f64 = 0.0;
            let mut y: f64 = 0.0;
        //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
            let a: f64 = vectors.sols_Ex[t][2];
            let d: f64 = vectors.sols_Ex[t][3];
            let k_x: f64 = vectors.sols_Ex[t][6];
            let k_y: f64 = vectors.sols_Ex[t][7];
            let gamma_x: f64 = vectors.sols_Ex[t][8];
            let gamma_y: f64 = vectors.sols_Ex[t][9];
            let beta: f64 = vectors.sols_Ex[t][11];
            let phi: f64 = E_x::phi(vectors.sols_Ex[t][0], PI);
            let psi: f64 = E_x::psi(vectors.sols_Ex[t][1], PI);

            let range_x  = Range::new(-1.0 * a, a);
            let range_y  = Range::new(-1.0 * d, d);
            let mut rng = rand::thread_rng();
            let N: i32 = 1000;
            let cycles: i32 = 1;
            //Region 1, -a<=x<=a, -d<=y<=d
            for z in 0..cycles {
                for i in 0..N {
                    y = range_y.ind_sample(&mut rng);
                    for j in 0..N {
                        x = range_x.ind_sample(&mut rng);
                        sum_x += (E_x::Hy1(1.0, k_x, k_y, phi, psi, x, y)).powi(2)

                    }
                    sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                    sum_x = 0.0;
                }
            }
            let P1: f64 = c1 * (2.0 * d) * sum_y / ((N as f64) * (cycles as f64));
            sum_x = 0.0;
            sum_y = 0.0;
            let range_x  = Range::new(a, 3.0 * a);
            let range_y  = Range::new(-1.0 * d, d);
            let mut rng = rand::thread_rng();
            //Region 2, a<=x<=3a, -d<=y<=d
            for z in 0..cycles {
                for i in 0..N {
                    y = range_y.ind_sample(&mut rng);
                    for j in 0..N {
                        x = range_x.ind_sample(&mut rng);
                        sum_x += (E_x::Hy2(1.0, k_x, k_y, gamma_x, phi, psi, x, y, a)).powi(2)

                    }
                    sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                    sum_x = 0.0;
                }
            }

            let P2: f64 = c2 * (2.0 * d) * sum_y / ((N as f64) * (cycles as f64));
            let range_x  = Range::new(-1.0 * a, a);
            let range_y  = Range::new(d, 3.0 * d);
            let mut rng = rand::thread_rng();
            sum_x = 0.0;
            sum_y = 0.0;
            //Region 3, -a<=x<=a, d<=y<=3d
            for z in 0..cycles {
                for i in 0..N {
                    y = range_y.ind_sample(&mut rng);
                    for j in 0..N {
                        x = range_x.ind_sample(&mut rng);
                        sum_x += (E_x::Hy3(1.0, k_x, k_y, gamma_y, phi, psi, x, y, d)).powi(2)

                    }
                    sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                    sum_x = 0.0;
                }
            }
            let P3: f64 = c3 * (2.0 * d) * sum_y / ((N as f64) * (cycles as f64));
            let containment = P1 / (P1+P2+P3);
            vectors.sols_Ex[t][12] = containment;
        }
    }
/*
println!("vectors.sols_Ex[0] = {:?}", vectors.sols_Ex);
println!("0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x");
*/

}

write!(result_file, "E_x mode,{},nm\n", (g * 1e9).round());
write!(result_file,
    "P, Q, a/2, d/2, n_0, n_1, kx, ky, gamma_x, gamma_y, beta, n_eff, cont_x\n");
for elements in vectors.sols_Ex.iter() {
    if elements[0] > 0.0 {
        write!(result_file, "{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
            elements[0], elements[1], elements[2], elements[3], elements[4],
            elements[5], elements[6], elements[7], elements[8], elements[9],
            elements[10], elements[11], elements[12]
        );

    }
}

//end of block
//0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
}

////////////////////////////////RESET////////////////////////////////////////
vectors.k_x = [0.0; P];
vectors.k_y = [0.0; Q];
vectors.gamma_x = [0.0; P];
vectors.gamma_y = [0.0; Q];
vectors.beta = [[0.0; P]; Q];

a = params.a;
d = params.d;
k = params.k;
n_0 = params.n_0;
n_1 = params.n_1;
stepsize = k * n_1 / (STEPSMODE as f64);
p = 0.0;
q = 0.0;
beta = 0.0;
gamma_x = 0.0;
k_x = 0.0;
gamma_y = 0.0;
k_y = 0.0;
n_eff = 0.0;
comparator = 10.0;
index = 0;
newlow = 1000000.0;

////////////////////////////////CALCULATION/////////////////////////////////
//////////////////////////////////E_y mode/////////////////////////////////
for w in 0..P {
newlow = 1000000.0;
for i in 0..STEPSMODE {
    /*gamma_x, k_x calculation*/
    k_x = (i as f64) * stepsize;
    gamma_x = E_y::gamma_x(n_0, n_1, k, k_x);
    p = w as f64 + 1.0;
    comparator = E_y::find_k_x(n_0, n_1, p, a, PI, gamma_x, k_x);
    if comparator < 1e-3 && comparator < newlow &&
    k_x > 10.0 {
        vectors.k_x[w] = k_x;
        vectors.gamma_x[w] = gamma_x;
        newlow = comparator;
    }
}
}

for w in 0..Q {
    newlow = 1000000.0;
    for i in 0..STEPSMODE {
        /*gamma_y, k_y calculation*/
        k_y = (i as f64) * stepsize;
        gamma_y = E_y::gamma_y(n_0, n_1, k, k_y);
        q = w as f64 + 1.0;
        comparator = E_y::find_k_y(n_0, n_1, q, d, PI, gamma_y, k_y);
        if comparator < 1e-3 && comparator < newlow &&
        k_y > 10.0 {
            vectors.k_y[w] = k_y;
            vectors.gamma_y[w] = gamma_y;
            newlow = comparator;
        }
    }
}

'first:for w in 0..P {
    'second:for e in 0..Q{
        k_x = vectors.k_x[w];
        k_y = vectors.k_y[e];
        gamma_x = vectors.gamma_x[w];
        gamma_y = vectors.gamma_y[e];
        beta = E_y::find_beta(n_0, n_1, k, k_x, k_y);
        n_eff = E_y::find_n_eff(beta, k);
        if (n_eff > n_0 && n_eff < n_1) && (k_x > 1.0 && k_y > 1.0) {
            vectors.beta[w][e] = beta;
            vectors.sols_Ey[index] = [(w+1) as f64, (e+1) as f64, a, d, n_0,
                n_1, k_x, k_y, gamma_x, gamma_y, beta, n_eff, 0.0];
                index += 1;
                //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
        }
    }
}

index = 0;
//start of block
{
let mut counter: i32 = 0;
for i in vectors.sols_Ey.iter() {
if i[0] != 0.0 {
    counter += 1
}

}
if counter > 0 {
for t in 0..vectors.sols_Ey.len() {
    if vectors.sols_Ey[t][6] > 0.5 {
        let c1: f64 = (params.w_ang * mu_0 * (vectors.sols_Ey[t][10].powi(2))
        - (vectors.sols_Ey[t][7].powi(2)))/(vectors.sols_Ey[t][10].powi(3));
        let c2: f64 = c1.clone();
        let c3: f64 = (params.w_ang * mu_0 * (vectors.sols_Ey[t][10].powi(2))
        + (vectors.sols_Ey[t][9].powi(2)))/(vectors.sols_Ey[t][10].powi(3));
        let mut sum_x: f64 = 0.0;
        let mut sum_y: f64 = 0.0;
        let mut x: f64 = 0.0;
        let mut y: f64 = 0.0;
    //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
        let a: f64 = vectors.sols_Ey[t][2];
        let d: f64 = vectors.sols_Ey[t][3];
        let k_x: f64 = vectors.sols_Ey[t][6];
        let k_y: f64 = vectors.sols_Ey[t][7];
        let gamma_x: f64 = vectors.sols_Ey[t][8];
        let gamma_y: f64 = vectors.sols_Ey[t][9];
        let beta: f64 = vectors.sols_Ey[t][11];
        let phi: f64 = E_y::phi(vectors.sols_Ey[t][0], PI);
        let psi: f64 = E_y::psi(vectors.sols_Ey[t][1], PI);

        let range_x  = Range::new(-1.0 * a, a);
        let range_y  = Range::new(-1.0 * d, d);
        let mut rng = rand::thread_rng();
        let N: i32 = 1000;
        let cycles: i32 = 1;
        //Region 1, -a<=x<=a, -d<=y<=d
        for z in 0..cycles {
            for i in 0..N {
                y = range_y.ind_sample(&mut rng);
                for j in 0..N {
                    x = range_x.ind_sample(&mut rng);
                    sum_x += (E_y::Hx1(1.0, k_x, k_y, phi, psi, x, y)).powi(2)

                }
                sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                sum_x = 0.0;
            }
        }
        let P1: f64 = c1 * (2.0 * d) * sum_y / ((N as f64) * (cycles as f64));
        sum_x = 0.0;
        sum_y = 0.0;
        let range_x  = Range::new(a, 3.0 * a);
        let range_y  = Range::new(-1.0 * d, d);
        let mut rng = rand::thread_rng();
        //Region 2, a<=x<=3a, -d<=y<=d
        for z in 0..cycles {
            for i in 0..N {
                y = range_y.ind_sample(&mut rng);
                for j in 0..N {
                    x = range_x.ind_sample(&mut rng);
                    sum_x += (E_y::Hx2(1.0, k_x, k_y, gamma_x, phi, psi, x, y, a)).powi(2)

                }
                sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                sum_x = 0.0;
            }
        }

        let P2: f64 = c2 * (2.0 * d) * sum_y / ((N as f64) * (cycles as f64));
        let range_x  = Range::new(-1.0 * a, a);
        let range_y  = Range::new(d, 3.0 * d);
        let mut rng = rand::thread_rng();
        sum_x = 0.0;
        sum_y = 0.0;
        //Region 3, -a<=x<=a, d<=y<=3d
        for z in 0..cycles {
            for i in 0..N {
                y = range_y.ind_sample(&mut rng);
                for j in 0..N {
                    x = range_x.ind_sample(&mut rng);
                    sum_x += (E_y::Hx3(1.0, k_x, k_y, gamma_y, phi, psi, x, y, d)).powi(2)

                }
                sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                sum_x = 0.0;
            }
        }
        let P3: f64 = c3 * (2.0 * d) * sum_y/((N as f64) * (cycles as f64));
        let containment = P1 / (P1+P2+P3);
        vectors.sols_Ey[t][12] = containment;
    }
//End of for t in 0..vectors.sols_Ey.len() block
}
//End of if counter > 0 block
}
write!(result_file, "E_y mode,{},nm\n", (g * 1e9).round());
write!(result_file,
"P, Q, a/2, d/2, n_0, n_1, kx, ky, gamma_x, gamma_y, beta, n_eff, cont_x\n");

for elements in vectors.sols_Ey.iter() {
if elements[0] > 0.0 {
    write!(result_file, "{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
        elements[0], elements[1], elements[2], elements[3], elements[4],
        elements[5], elements[6], elements[7], elements[8], elements[9],
        elements[10], elements[11], elements[12]
    );
}
}
//End of block
}
//0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
completion += 1;
//End of for loop values_to_iterate
}
drop(result_file);
println!("Mode calculation completed")
//End of else block
}
drop(size_of_values);
drop(file);
//End of mode block
}

if val4[2] == true
//Start of gain block
{
/////////////////////////////////Reading the file///////////////////////////////
println!("Starting gain calculation");
println!("________________________________________________");
let path = Path::new("mode_results.csv");
let display = path.display();
let mut results_file = match File::open(path) {
    Ok(file) => file,
    Err(why) => panic!("Unable to open file {:?}, {}", display, why.description()),
};

let mut results_contents: String = String::new();
match results_file.read_to_string(&mut results_contents) {
    Ok(_) => println!("Results file read successfully!"),
    Err(why) => panic!("Couldn't read file {:?}, {}", display, why.description()),
}

let val5: Vec<&str> = results_contents.split("\n").collect();

let mut val6: Vec<f64> = vec![];
let mut val7: Vec<&str> = vec![];
for elements in val5.into_iter() {
    for parts in elements.split(",") {
        match parts.trim().parse::<f64>() {
            Ok(value) => val6.push(value),
            Err(_) => val6.push(-1.0),
        }
        match parts.trim().parse::<f64>() {
            Ok(value) => continue,
            Err(_) => val7.push(parts.trim()),
        }

    }
}
val6.push(-1.0);
let mut containment_vector: Vec<f64> = vec![];
let mut fundamental_containment: Vec<f64> = vec![];
let mut fundamental_containment_wavelength: Vec<f64> = vec![];
//Start of containment block
{
let mut start_position: usize = 0;
let mut end_position: usize = 0;
for i in 0..val6.len()-1 {
    if i > 0 && (val6[i].round() as i32) == -1 {
        if val6[i-1].round() as i32 != -1 {
            if start_position == 0 {
                start_position = i;
            } else {
                end_position = i-1;
                containment_vector.push(val6[start_position-1]);
                if (start_position + 13) < val6.len() {
                    start_position += 13;
                    while start_position < end_position {
                        start_position += 13;
                        containment_vector.push(val6[start_position])
                    }
                }
                start_position = 0;
                end_position = 0
            }


            //println!("val6[{}] = {}", i-1, val6[i-1] )
        }
    }
}
let mut comparator: [f64; 2] = [0.0,0.0];
for i in 0..containment_vector.len()-1 {
    if containment_vector[i] > 1.1 && (comparator[0].round() as i32) == 0 {
        comparator[0] = containment_vector[i+1];
        fundamental_containment_wavelength.push(containment_vector[i])
    } else if containment_vector[i] > 1.1 && comparator[0] > 0.0 {
        comparator[1] = containment_vector[i+1];
        if comparator [0] < comparator [1] {
            fundamental_containment.push(comparator[0]);
            comparator = [0.0;2]
        } else {
            fundamental_containment.push(comparator[1]);
            comparator = [0.0;2]
        }
    }
}


//end of containment block
}

drop(val6);
drop(val7);


//Finds the containment factor
let find_containment = |wavelength: f64| -> f64 {
    let wavel: f64 = wavelength * 1.0e9;
    let size: f64 = fundamental_containment.len() as f64;
    let size0: usize = ((size / 2.0).round()) as usize;
    let size1: usize = (size0 + 1) as usize;
    let x0: f64 = fundamental_containment_wavelength[size0];
    let x1: f64 = fundamental_containment_wavelength[size1];
    let y0: f64 = fundamental_containment[size0];
    let y1: f64 = fundamental_containment[size1];
    let m: f64 = (y1 - y0) / (x1 - x0);
    m * (wavel - x0) + y0
};

////////////////////////////////Initializing parameters for gain/////////////////////
let mut params: starting_params = starting_params {
    s_freq: vec![0.0; SIGNALS],
    s_start: val3[0] * 1.0e-9,
    s_wavelength: vec![0.0; SIGNALS],                                     //All signals of this wavelength
    s_sig: vec![vec![0.0, 0.0];SIGNALS],
    s_gam: vec![0.0; SIGNALS],
    s_alp: vec![0.0; SIGNALS],
    fp_freq: vec![0.0; FPUMPS],
    fp_wavelength: vec![val3[1] * 1.0e-9;FPUMPS],
    fp_sig: vec![vec![0.0, 0.0]; FPUMPS],
    fp_gam: vec![0.0; FPUMPS],
    fp_alp: vec![0.0; FPUMPS],
    bp_freq: vec![0.0; BPUMPS],
    bp_wavelength: vec![val3[2] * 1.0e-9;BPUMPS],
    bp_sig: vec![vec![0.0, 0.0]; BPUMPS],
    bp_gam: vec![0.0; BPUMPS],
    bp_alp: vec![0.0; BPUMPS],
    width: val3[3] * 1.0e-6,                                                          //width in m
    depth: val3[4] * 1.0e-6,                                                          //depth in m
    length: val3[5] * 1.0e-2,                                                        //length in m
    area: 0.0,
    N: val3[6],
    tau: vec![val3[7] * 1.0e-3; LIFETIMES],
    fp0: vec![val3[8];FPUMPS],                                                   // Forward Pump power at length 0, Watts
    bp0: vec![val3[9];BPUMPS],                                                      // Backward Pump power at length L, Watts
    s0: vec![val3[10];SIGNALS],                                                   // Signal power at length 0, Watts
};

for i in 0..SIGNALS {
    let wavelength = params.s_start;
    params.s_wavelength[i] = params.s_start +
    ((i as f64)) * 1.0e-9;
    params.s_freq[i] = find_freq(params.s_wavelength[i]);
    params.s_sig[i][0] = find_absorption(params.s_wavelength[i]);
    params.s_sig[i][1] = find_emission(params.s_wavelength[i]);
    params.s_gam[i] = fundamental_containment[i+2];
}


for i in 0..FPUMPS {
    params.fp_freq[i] = find_freq(params.fp_wavelength[i]);
    params.fp_sig[i][0] = find_absorption(params.fp_wavelength[i]);
    params.fp_sig[i][1] = find_emission(params.fp_wavelength[i]);
    params.fp_gam[i] = fundamental_containment[0];
}

for i in 0..BPUMPS {
    params.bp_freq[i] = find_freq(params.bp_wavelength[i]);
    params.bp_sig[i][0] = find_absorption(params.bp_wavelength[i]);
    params.bp_sig[i][1] = find_emission(params.bp_wavelength[i]);
    params.bp_gam[i] = fundamental_containment[1];
}

params.area = find_area(params.width, params.depth);                        //Calculate the area of the cross section

let mut ase: ASE_signals = ASE_signals {
    fASEstart_wavelength: val3[11] * 1.0e-9,
    channel_spacing: val3[12] * 1.0,                                        //125 GHz = 1 nm,  1/s
    fASE_alp: vec![0.0; FASE],
    fASE_sig: vec![vec![0.0, 0.0]; FASE],                                           //Assumed same for all signals
    fASE_wavelength: vec![0.0; FASE],
    fASE_freq: vec![0.0; FASE],
    fASE_gam: vec![0.0; FASE],                                                  // Assumed same overlap as signal
    fASE0: vec![0.0; FASE],                                                     //ASE power at z = 0
    bASEstart_wavelength: val3[13] * 1.0e-9,
    bASE_alp: vec![0.0; BASE],                                                  //No attenuation
    bASE_sig: vec![vec![0.0, 0.0]; BASE],
    bASE_wavelength: vec![0.0; BASE],
    bASE_freq: vec![0.0; BASE],
    bASE_gam: vec![0.0; BASE],
    bASE0: vec![0.0; BASE],                                                     //ASE power a z = l
};

for i in 0..FASE {
    ase.fASE_wavelength[i] =
    ase.fASEstart_wavelength + ((i as f64) * ase.channel_spacing/125.0e9) * 1.0e-9;
    ase.fASE_freq[i] = find_freq(ase.fASE_wavelength[i]);
    ase.fASE_sig[i][0] = find_absorption(ase.fASE_wavelength[i]);
    ase.fASE_sig[i][1] = find_emission(ase.fASE_wavelength[i]);
    for w in 0..fundamental_containment_wavelength.len() {
        if ase.fASE_wavelength[i] * 1.0e9
        - fundamental_containment_wavelength[w] < 0.01 {
            ase.fASE_gam[i] = fundamental_containment[w]
        } else {
            ase.fASE_gam[i] = find_containment(ase.fASE_wavelength[i]);
        }
    }
    ase.fASE0[i] = 0.0;
    //ase.fASE[i] = 2.0 * h * ase.fASE_freq[i] * ase.channel_spacing;        //Seed noise
}


for i in 0..BASE {
    ase.bASE_wavelength[i] =
    ase.bASEstart_wavelength + ((i as f64) * ase.channel_spacing/125.0e9) * 1.0e-9;
    ase.bASE_freq[i] = find_freq(ase.bASE_wavelength[i]);
    ase.bASE_sig[i][0] = find_absorption(ase.bASE_wavelength[i]);
    ase.bASE_sig[i][1] = find_emission(ase.bASE_wavelength[i]);
    for w in 0..fundamental_containment_wavelength.len() {
        if ase.bASE_wavelength[i] * 1.0e9
        - fundamental_containment_wavelength[w] < 0.01 {
            ase.bASE_gam[i] = fundamental_containment[w]
        } else {
            ase.bASE_gam[i] = find_containment(ase.bASE_wavelength[i]);
        }
    }
    ase.bASE0[i] = 0.0;
    //ase.bASE[i] = 2.0 * h * ase.bASE_freq[i] * ase.channel_spacing;        //Seed noise
}

let mut vectors: vectors = vectors {
    signals:vec![vec![0.0; SIGNALS]; STEPS],
    fpumps: vec![vec![0.0; FPUMPS]; STEPS],
    bpumps: vec![vec![0.0; BPUMPS]; STEPS],
    fASE: vec![vec![0.0; FASE]; STEPS],
    bASE: vec![vec![0.0; BASE]; STEPS],
    sums: vec![vec![0.0; 6]; STEPS],
    position: vec![0.0; STEPS],
    rsignals: vec![vec![0.0; SIGNALS]; STEPS],
    rfpumps: vec![vec![0.0; FPUMPS]; STEPS],
    rbpumps: vec![vec![0.0; BPUMPS]; STEPS],
    rfASE: vec![vec![0.0; FASE]; STEPS],
    rbASE: vec![vec![0.0; BASE]; STEPS],
    rsums: vec![vec![0.0; 6]; STEPS],

};

for i in 0..SIGNALS {
    vectors.signals[0][i] = params.s0[i];
    vectors.sums[0][1] += params.s0[i];
}

for i in 0..FPUMPS {
    vectors.fpumps[0][i] = params.fp0[i];
    vectors.sums[0][2] += params.fp0[i];
}

for i in 0..BPUMPS {
    vectors.bpumps[STEPS - 1][i] = params.bp0[i];
    vectors.sums[STEPS - 1][3] += params.bp0[i];
}

for i in 0..FASE {
    vectors.fASE[0][i] = ase.fASE0[i];
    vectors.sums[0][4] += ase.fASE0[i];
}

for i in 0..BASE {
    vectors.bASE[STEPS - 1][i] = ase.bASE0[i];
    vectors.sums[STEPS - 1][5] += ase.bASE0[i];
}

let mut ct: constants = constants {
    A: vec![0.0; SIGNALS],
    D: vec![0.0; SIGNALS],
    H: vec![0.0; SIGNALS],
    fC: vec![0.0; FPUMPS],
    fF: vec![0.0; FPUMPS],
    fG: vec![0.0; FPUMPS],
    bC: vec![0.0; BPUMPS],
    bF: vec![0.0; BPUMPS],
    bG: vec![0.0; BPUMPS],
    fB: vec![0.0; FASE],
    fE: vec![0.0; FASE],
    fI: vec![0.0; FASE],
    fJ: vec![0.0; FASE],
    bB: vec![0.0; BASE],
    bE: vec![0.0; BASE],
    bI: vec![0.0; BASE],
    bJ: vec![0.0; BASE],
};

for i in 0..SIGNALS {
    ct.A[i] = params.tau[0] * params.s_sig[i][0] * params.s_gam[i] /
            (params.area * h * params.s_freq[i]);
    ct.D[i] = params.tau[0] * (params.s_sig[i][0] + params.s_sig[i][1]) *
            params.s_gam[i] / (params.area * h * params.s_freq[i]);
    ct.H[i] = params.s_sig[i][0] + params.s_sig[i][1];
}

for i in 0..FPUMPS {
    ct.fC[i] = params.tau[0] * params.fp_sig[i][0] * params.fp_gam[i] /
            (params.area * h * params.fp_freq[i]);
    ct.fF[i] = params.tau[0] * (params.fp_sig[i][0] + params.fp_sig[i][1]) *
            params.fp_gam[i] / (params.area * h * params.fp_freq[i]);
    ct.fG[i] = params.fp_sig[i][0] + params.fp_sig[i][1];
}

for i in 0..BPUMPS {
    ct.bC[i] = params.tau[0] * params.bp_sig[i][0] * params.bp_gam[i] /
            (params.area * h * params.bp_freq[i]);
    ct.bF[i] = params.tau[0] * (params.bp_sig[i][0] + params.bp_sig[i][1]) *
            params.bp_gam[i] / (params.area * h * params.bp_freq[i]);
    ct.bG[i] = params.bp_sig[i][0] + params.bp_sig[i][1];
}

for i in 0..FASE {
    ct.fB[i] = params.tau[0] * ase.fASE_sig[i][0] * ase.fASE_gam[i] /
            (params.area * h * ase.fASE_freq[i]);
    ct.fE[i] = params.tau[0] * (ase.fASE_sig[i][0] + ase.fASE_sig[i][1]) *
            ase.fASE_gam[i] / (params.area * h * ase.fASE_freq[i]);
    ct.fI[i] = ase.fASE_sig[i][1] * ase.fASE_gam[i] * ase.channel_spacing *
            h * ase.fASE_freq[i];
    ct.fJ[i] = ase.fASE_sig[i][0] + ase.fASE_sig[i][1];
    /*
    assert!(ct.fB[i]>0.0);
    assert!(ct.fE[i]>0.0);
    assert!(ct.fI[i]>0.0);
    assert!(ct.fJ[i]>0.0);
    */
}

for i in 0..BASE {
    ct.bB[i] = params.tau[0] * ase.bASE_sig[i][0] * ase.bASE_gam[i] /
            (params.area * h * ase.bASE_freq[i]);
    ct.bE[i] = params.tau[0] * (ase.bASE_sig[i][0] + ase.bASE_sig[i][1]) *
            ase.bASE_gam[i] / (params.area * h * ase.bASE_freq[i]);
    ct.bI[i] = ase.bASE_sig[i][1] * ase.bASE_gam[i] * ase.channel_spacing *
            h * ase.bASE_freq[i];
    ct.bJ[i] = ase.bASE_sig[i][0] + ase.bASE_sig[i][1];
}

let mut rk: rk = rk {
    c1: 1.0/6.0,
    c2: 1.0/3.0,
    c3: 1.0/3.0,
    c4: 1.0/6.0,
};

let stepsize: f64 = params.length / (STEPS as f64 - 1.0);
let mut ztot: f64 = 0.0;
let mut APS: f64 = 0.0;
let mut BPA: f64 = 0.0;
let mut CPP: f64 = 0.0;
let mut DPS: f64 = 0.0;
let mut EPA: f64 = 0.0;
let mut FPP: f64 = 0.0;
let mut N20: f64 = 0.0;
let mut dir: f64 = 1.0;
let mut k1: f64 = 0.0;
let mut k2: f64 = 0.0;
let mut k3: f64 = 0.0;
let mut k4: f64 = 0.0;
let mut rel_err: f64 = 0.0;

for z in 0..SIGNALS {
    APS += ct.A[z] * vectors.signals[0][z];
    DPS += ct.D[z] * vectors.signals[0][z];
}

for z in 0..FASE {
    BPA += ct.fB[z] * vectors.fASE[0][z];
    EPA += ct.fE[z] * vectors.fASE[0][z];
}

for z in 0..BASE {
    BPA += ct.bB[z] * vectors.bASE[0][z];
    EPA += ct.bE[z] * vectors.bASE[0][z];
}

for z in 0..FPUMPS {
    CPP += ct.fC[z] * vectors.fpumps[0][z];
    FPP += ct.fF[z] * vectors.fpumps[0][z];
}

for z in 0..BPUMPS {
    CPP += ct.bC[z] * vectors.bpumps[0][z];
    FPP += ct.bF[z] * vectors.bpumps[0][z];
}

N20 = N2(params.N, APS, BPA, CPP, DPS, EPA, FPP);
vectors.sums[0][0] = N20;
let mut N21: f64 = N20.clone();
vectors.sums[0][0] = N20;
vectors.sums[1][0] = N21;
println!("Starting N20 = {:?}", N20/params.N);

for i in 0..(STEPS-1) {
    vectors.position[i+1] = stepsize + i as f64 * stepsize;
}


///////////////////////////////CALCULATION LOOP/////////////////////////////////
for n in 0..CYCLES{
    println!("cycle # {}", n+1);
    for i in 0..(STEPS-1) {
        for j in 0..MATCHING {
            N20 = vectors.sums[i][0];
            N21 = vectors.sums[i+1][0];

            for l in 0..SIGNALS {

                k1 = dir * dPs(
                    params.N, ct.H[l], params.s_sig[l][0], params.s_gam[l],
                    params.s_alp[l], N21,
                    vectors.signals[i][l]
                );
                k2 = dir * dPs(
                    params.N, ct.H[l], params.s_sig[l][0], params.s_gam[l],
                    params.s_alp[l], N21,
                    (vectors.signals[i][l] + k1 * 0.5 * stepsize)
                );
                k3 = dir * dPs(
                    params.N, ct.H[l], params.s_sig[l][0], params.s_gam[l],
                    params.s_alp[l], N21,
                    (vectors.signals[i][l] + k2 * 0.5 * stepsize)
                );
                k4 = dir * dPs(
                    params.N, ct.H[l], params.s_sig[l][0], params.s_gam[l],
                    params.s_alp[l], N21,
                    (vectors.signals[i][l] + k3 * stepsize)
                );
                vectors.signals[i+1][l] = vectors.signals[i][l] +
                (rk.c1 * k1 + rk.c2 * k2 + rk.c3 * k3 + rk.c4 * k4) * stepsize;

                if vectors.signals[i+1][l] < 0.0 {vectors.signals[i+1][l] = 0.0}


            }

            for l in 0..FPUMPS {

                k1 = dir * dPp(
                    params.N, ct.fG[l], params.fp_sig[l][0], params.fp_gam[l],
                    params.fp_alp[l], N21,
                    vectors.fpumps[i][l]
                );

                k2 = dir * dPp(
                    params.N, ct.fG[l], params.fp_sig[l][0], params.fp_gam[l],
                    params.fp_alp[l], N21,
                    (vectors.fpumps[i][l] + k1 * 0.5 * stepsize)
                );

                k3 = dir * dPp(
                    params.N, ct.fG[l], params.fp_sig[l][0], params.fp_gam[l],
                    params.fp_alp[l], N21,
                    (vectors.fpumps[i][l] + k2 * 0.5 * stepsize)
                );

                k4 = dir * dPp(
                    params.N, ct.fG[l], params.fp_sig[l][0], params.fp_gam[l],
                    params.fp_alp[l], N21,
                    (vectors.fpumps[i][l] + k3 * stepsize)
                );

                vectors.fpumps[i+1][l] = vectors.fpumps[i][l] +
                (rk.c1 * k1 + rk.c2 * k2 + rk.c3 * k3 + rk.c4 * k4) * stepsize;

                if vectors.fpumps[i+1][l] < 0.0 {vectors.fpumps[i+1][l] = 0.0}

            }

            for l in 0..FASE {

                k1 = dir * dPASE(
                    params.N, ct.fJ[l], ct.fI[l], ase.fASE_sig[l][0], ase.fASE_gam[l],
                    ase.fASE_alp[l], N21,
                    vectors.fASE[i][l]
                );

                k2 = dir * dPASE(
                    params.N, ct.fJ[l], ct.fI[l], ase.fASE_sig[l][0], ase.fASE_gam[l],
                    ase.fASE_alp[l], N21,
                    (vectors.fASE[i][l] + k1 * 0.5 * stepsize)
                );

                k3 = dir * dPASE(
                    params.N, ct.fJ[l], ct.fI[l], ase.fASE_sig[l][0], ase.fASE_gam[l],
                    ase.fASE_alp[l], N21,
                    (vectors.fASE[i][l] + k2 * 0.5 * stepsize)
                );

                k4 = dir * dPASE(
                    params.N, ct.fJ[l], ct.fI[l], ase.fASE_sig[l][0], ase.fASE_gam[l],
                    ase.fASE_alp[l], N21,
                    (vectors.fASE[i][l] + k3 * stepsize)
                );

                vectors.fASE[i+1][l] = vectors.fASE[i][l] +
                (rk.c1 * k1 + rk.c2 * k2 + rk.c3 * k3 + rk.c4 * k4) * stepsize;

                if vectors.fASE[i+1][l] < 0.0 {vectors.fASE[i+1][l] = 0.0}

            }

            APS = 0.0;
            DPS = 0.0;
            BPA = 0.0;
            EPA = 0.0;
            CPP = 0.0;
            FPP = 0.0;

            for z in 0..SIGNALS {
                APS += ct.A[z] * vectors.signals[i+1][z];
                DPS += ct.D[z] * vectors.signals[i+1][z];
            }

            for z in 0..FASE {
                BPA += ct.fB[z] * vectors.fASE[i+1][z];
                EPA += ct.fE[z] * vectors.fASE[i+1][z];
            }

            for z in 0..BASE {
                BPA += ct.bB[z] * vectors.bASE[i+1][z];
                EPA += ct.bE[z] * vectors.bASE[i+1][z];
            }

            for z in 0..FPUMPS {
                CPP += ct.fC[z] * vectors.fpumps[i+1][z];
                FPP += ct.fF[z] * vectors.fpumps[i+1][z];
            }

            for z in 0..BPUMPS {
                CPP += ct.bC[z] * vectors.bpumps[i+1][z];
                FPP += ct.bF[z] * vectors.bpumps[i+1][z];
            }


            N21 = N2(params.N, APS, BPA, CPP, DPS, EPA, FPP);

            vectors.sums[i][0] = N20;
            vectors.sums[i+1][0] = N2(params.N, APS, BPA, CPP, DPS, EPA, FPP);

        }


        }

    vectors.sums.reverse();
    vectors.signals.reverse();
    vectors.fpumps.reverse();
    vectors.fASE.reverse();
    vectors.bpumps.reverse();
    vectors.bASE.reverse();


    vectors.bpumps[0] = params.bp0.clone();
    vectors.bASE[0] = vec![0.0; BASE];
    //println!("vectors.base = {:?}", vectors.bASE);
    APS = 0.0;
    DPS = 0.0;
    BPA = 0.0;
    EPA = 0.0;
    CPP = 0.0;
    FPP = 0.0;

    for z in 0..SIGNALS {
        APS += ct.A[z] * vectors.signals[0][z];
        DPS += ct.D[z] * vectors.signals[0][z];
    }

    for z in 0..FASE {
        BPA += ct.fB[z] * vectors.fASE[0][z];
        EPA += ct.fE[z] * vectors.fASE[0][z];
    }

    for z in 0..BASE {
        BPA += ct.bB[z] * vectors.bASE[0][z];
        EPA += ct.bE[z] * vectors.bASE[0][z];
    }

    for z in 0..FPUMPS {
        CPP += ct.fC[z] * vectors.fpumps[0][z];
        FPP += ct.fF[z] * vectors.fpumps[0][z];
    }

    for z in 0..BPUMPS {
        CPP += ct.bC[z] * vectors.bpumps[0][z];
        FPP += ct.bF[z] * vectors.bpumps[0][z];

    }

    N20 = N2(params.N, APS, BPA, CPP, DPS, EPA, FPP);
    vectors.sums[0][0] = N20;

    APS = 0.0;
    DPS = 0.0;
    BPA = 0.0;
    EPA = 0.0;
    CPP = 0.0;
    FPP = 0.0;

    for z in 0..SIGNALS {
        APS += ct.A[z] * vectors.signals[1][z];
        DPS += ct.D[z] * vectors.signals[1][z];
    }

    for z in 0..FASE {
        BPA += ct.fB[z] * vectors.fASE[1][z];
        EPA += ct.fE[z] * vectors.fASE[1][z];
    }

    for z in 0..BASE {
        BPA += ct.bB[z] * vectors.bASE[1][z];
        EPA += ct.bE[z] * vectors.bASE[1][z];
    }

    for z in 0..FPUMPS {
        CPP += ct.fC[z] * vectors.fpumps[1][z];
        FPP += ct.fF[z] * vectors.fpumps[1][z];
    }

    for z in 0..BPUMPS {
        CPP += ct.bC[z] * vectors.bpumps[1][z];
        FPP += ct.bF[z] * vectors.bpumps[1][z];
    }

    N21 = N2(params.N, APS, BPA, CPP, DPS, EPA, FPP);
    vectors.sums[1][0] = N21;

    for i in 0..(STEPS-1) {
        for j in 0..MATCHING {
            N20 = vectors.sums[i][0];
            N21 = vectors.sums[i+1][0];

            for l in 0..BPUMPS {

                k1 = dir * dPp(
                    params.N, ct.bG[l], params.bp_sig[l][0], params.bp_gam[l],
                    params.bp_alp[l], N21,
                    vectors.bpumps[i][l]
                );

                k2 = dir * dPp(
                    params.N, ct.bG[l], params.bp_sig[l][0], params.bp_gam[l],
                    params.bp_alp[l], N21,
                    (vectors.bpumps[i][l] + k1 * 0.5 * stepsize)
                );

                k3 = dir * dPp(
                    params.N, ct.bG[l], params.bp_sig[l][0], params.bp_gam[l],
                    params.bp_alp[l], N21,
                    (vectors.bpumps[i][l] + k2 * 0.5 * stepsize)
                );

                k4 = dir * dPp(
                    params.N, ct.bG[l], params.bp_sig[l][0], params.bp_gam[l],
                    params.bp_alp[l], N21,
                    (vectors.bpumps[i][l] + k3 * stepsize)
                );

                vectors.bpumps[i+1][l] = vectors.bpumps[i][l] +
                (rk.c1 * k1 + rk.c2 * k2 + rk.c3 * k3 + rk.c4 * k4) * stepsize;

                if vectors.bpumps[i+1][l] < 0.0 {vectors.bpumps[i+1][l] = 0.0}

            }

            for l in 0..BASE {

                k1 = dir * dPASE(
                    params.N, ct.bJ[l], ct.bI[l], ase.bASE_sig[l][0], ase.bASE_gam[l],
                    ase.bASE_alp[l], N21,
                    vectors.bASE[i][l]
                );

                k2 = dir * dPASE(
                    params.N, ct.bJ[l], ct.bI[l], ase.bASE_sig[l][0], ase.bASE_gam[l],
                    ase.bASE_alp[l], N21,
                    (vectors.bASE[i][l] + k1 * 0.5 * stepsize)
                );

                k3 = dir * dPASE(
                    params.N, ct.bJ[l], ct.bI[l], ase.bASE_sig[l][0], ase.bASE_gam[l],
                    ase.bASE_alp[l], N21,
                    (vectors.bASE[i][l] + k2 * 0.5 * stepsize)
                );

                k4 = dir * dPASE(
                    params.N, ct.bJ[l], ct.bI[l], ase.bASE_sig[l][0], ase.bASE_gam[l],
                    ase.bASE_alp[l], N21,
                    (vectors.bASE[i][l] + k3 * stepsize)
                );

                vectors.bASE[i+1][l] = vectors.bASE[i][l] +
                (rk.c1 * k1 + rk.c2 * k2 + rk.c3 * k3 + rk.c4 * k4) * stepsize;
                if vectors.bASE[i+1][l] < 0.0 {vectors.bASE[i+1][l] = 0.0}
            }

            APS = 0.0;
            DPS = 0.0;
            BPA = 0.0;
            EPA = 0.0;
            CPP = 0.0;
            FPP = 0.0;

            for z in 0..SIGNALS {
                APS += ct.A[z] * vectors.signals[i+1][z];
                DPS += ct.D[z] * vectors.signals[i+1][z];
            }

            for z in 0..FASE {
                BPA += ct.fB[z] * vectors.fASE[i+1][z];
                EPA += ct.fE[z] * vectors.fASE[i+1][z];
            }

            for z in 0..BASE {
                BPA += ct.bB[z] * vectors.bASE[i+1][z];
                EPA += ct.bE[z] * vectors.bASE[i+1][z];
            }

            for z in 0..FPUMPS {
                CPP += ct.fC[z] * vectors.fpumps[i+1][z];
                FPP += ct.fF[z] * vectors.fpumps[i+1][z];
            }

            for z in 0..BPUMPS {
                CPP += ct.bC[z] * vectors.bpumps[i+1][z];
                FPP += ct.bF[z] * vectors.bpumps[i+1][z];
            }

            N21 = N2(params.N, APS, BPA, CPP, DPS, EPA, FPP);
            vectors.sums[i][0] = N20;
            vectors.sums[i+1][0] = N2(params.N, APS, BPA, CPP, DPS, EPA, FPP);

        }

    }


    vectors.sums.reverse();
    vectors.signals.reverse();
    vectors.fpumps.reverse();
    vectors.bpumps.reverse();
    vectors.fASE.reverse();
    vectors.bASE.reverse();

    vectors.signals[0] = params.s0.clone();
    vectors.fpumps[0] = params.fp0.clone();
    vectors.fASE[0] = vec![0.0; FASE];

    APS = 0.0;
    DPS = 0.0;
    BPA = 0.0;
    EPA = 0.0;
    CPP = 0.0;
    FPP = 0.0;

    for z in 0..SIGNALS {
        APS += ct.A[z] * vectors.signals[0][z];
        DPS += ct.D[z] * vectors.signals[0][z];
    }

    for z in 0..FASE {
        BPA += ct.fB[z] * vectors.fASE[0][z];
        EPA += ct.fE[z] * vectors.fASE[0][z];
    }

    for z in 0..BASE {
        BPA += ct.bB[z] * vectors.bASE[0][z];
        EPA += ct.bE[z] * vectors.bASE[0][z];
    }

    for z in 0..FPUMPS {
        CPP += ct.fC[z] * vectors.fpumps[0][z];
        FPP += ct.fF[z] * vectors.fpumps[0][z];
    }

    for z in 0..BPUMPS {
        CPP += ct.bC[z] * vectors.bpumps[0][z];
        FPP += ct.bF[z] * vectors.bpumps[0][z];
    }

    N20 = N2(params.N, APS, BPA, CPP, DPS, EPA, FPP);
    vectors.sums[0][0] = N20;

    APS = 0.0;
    DPS = 0.0;
    BPA = 0.0;
    EPA = 0.0;
    CPP = 0.0;
    FPP = 0.0;

    for z in 0..SIGNALS {
        APS += ct.A[z] * vectors.signals[1][z];
        DPS += ct.D[z] * vectors.signals[1][z];
    }

    for z in 0..FASE {
        BPA += ct.fB[z] * vectors.fASE[1][z];
        EPA += ct.fE[z] * vectors.fASE[1][z];
    }

    for z in 0..BASE {
        BPA += ct.bB[z] * vectors.bASE[1][z];
        EPA += ct.bE[z] * vectors.bASE[1][z];
    }

    for z in 0..FPUMPS {
        CPP += ct.fC[z] * vectors.fpumps[1][z];
        FPP += ct.fF[z] * vectors.fpumps[1][z];
    }

    for z in 0..BPUMPS {
        CPP += ct.bC[z] * vectors.bpumps[1][z];
        FPP += ct.bF[z] * vectors.bpumps[1][z];
    }

    N21 = N2(params.N, APS, BPA, CPP, DPS, EPA, FPP);
    vectors.sums[1][0] = N21;

}

/////////////////////////////////PLOTS//////////////////////////////////////////


let mut N2_sum: Vec<f64> = vec![0.0; STEPS];
let mut signal_input: f64 = 0.0;
for i in 0..SIGNALS {
    signal_input += params.s0[i]
}

let mut signal_sum_lowest: Vec<f64> = vec![0.0; STEPS];
let mut signal_sum: Vec<f64> = vec![0.0; STEPS];
let mut fpump_sum: Vec<f64> = vec![0.0; STEPS];
let mut bpump_sum: Vec<f64> = vec![0.0; STEPS];
let mut fASE_sum: Vec<f64> = vec![0.0; STEPS];
let mut bASE_sum: Vec<f64> = vec![0.0; STEPS];
let mut gain: Vec<f64> = vec![0.0; STEPS];
let mut rolling_gain: Vec<f64> = vec![0.0; STEPS];
let mut fASE_noise: Vec<f64> = vec![0.0; FASE];
let mut fASE_noise_new: Vec<f64> = vec![0.0; FASE];
let mut fASE_power: Vec<f64> = vec![0.0; FASE];
let mut fASE_power_new: Vec<f64> = vec![0.0; FASE];
let mut fASE_spectra: Vec<f64> = vec![0.0; FASE];                            //Spectra in nm
let mut signal_spectra: Vec<f64> = vec![0.0; SIGNALS];
let mut signal_for_noise: Vec<f64> = vec![0.0;FASE];
let mut spectral_gain: Vec<f64> = vec![0.0; SIGNALS];
let mut fASE_spectral_gain: Vec<f64> = vec![];
let mut low_signal: usize = 0;
let mut comparator: f64 = 100000.0;


for i in 0..SIGNALS {
    signal_spectra[i] = params.s_wavelength[i] * 1.0e9;
    spectral_gain[i] = find_gain(vectors.signals[STEPS-1][i], params.s0[i]);
    if spectral_gain[i] < comparator {
        comparator = spectral_gain[i];
        low_signal = i;
    }
}

for i in 0..STEPS {
    for w in 1..6 {
        vectors.sums[i][w] = 0.0
    }
}


for i in 0..STEPS {
    for w in 0..SIGNALS {
        if i == 0 {vectors.sums[i][1] += params.s0[w]
        } else {vectors.sums[i][1] += vectors.signals[i][w]}
    }
    for w in 0..FPUMPS {
        if i == 0 {vectors.sums[i][2] = vectors.fpumps[i][w]
        } else {vectors.sums[i][2] += vectors.fpumps[i][w]}
    }
    for w in 0..BPUMPS {
        if i == 0 {vectors.sums[i][3] = vectors.bpumps[i][w]
        } else {vectors.sums[i][3] += vectors.bpumps[i][w]}
    }
    for w in 0..FASE {
        vectors.sums[i][4] += vectors.fASE[i][w]
    }
    for w in 0..BASE {
        vectors.sums[i][5] += vectors.bASE[i][w]
    }

}

println!("vectors.sums[0][1] = {:?}", vectors.sums[0][1]);
println!("low_signal gain={:?}", comparator);
println!("low_signal index={:?}, wavelength = {:?} nm", low_signal,
        params.s_wavelength[low_signal] * 1.0e9 );

for i in 0..STEPS {
    N2_sum[i] = vectors.sums[i][0] * 100.0 / params.N;
    signal_sum[i] = vectors.sums[i][1] * 1000.0;
    signal_sum_lowest[i] = vectors.signals[i][low_signal as usize] * 1000.0;
    fpump_sum[i] = vectors.sums[i][2] * 1000.0;
    bpump_sum[i] = vectors.sums[i][3] * 1000.0;
    fASE_sum[i] = vectors.sums[i][4] * 1000.0;
    bASE_sum[i] = vectors.sums[i][5] * 1000.0;
    vectors.position[i] = vectors.position[i] * 100.0;
    gain[i] = find_gain(vectors.signals[i][low_signal], vectors.signals[0][low_signal]);
    if i == 0 {rolling_gain[i] = 0.0} else {
    rolling_gain[i] = find_gain(vectors.signals[i][low_signal], vectors.signals[i-1][low_signal])}
}

/*Noise Figure Calculation */
let mut index: usize = 3;
for i in 0..FASE {
    fASE_spectra[i] = ase.fASE_wavelength[i] * 1.0e9;
    if i < 4 {
        fASE_spectral_gain.push(spectral_gain[i]);
        signal_for_noise[i] = vectors.signals[STEPS-1][i]
    } else if i % 3 != 0 {
        fASE_spectral_gain.push(spectral_gain[index]);
        signal_for_noise[i] = vectors.signals[STEPS-1][index]
    } else {
        fASE_spectral_gain.push(spectral_gain[index]);
        signal_for_noise[i] = vectors.signals[STEPS-1][index];
        index += 1
    }
}

for i in 0..FASE {
    fASE_power[i] = vectors.fASE[STEPS-1][i] * 1000.0;
    let v: f64 = find_freq(fASE_spectra[i] * 1.0e-9);
    let deltav: f64 = ase.channel_spacing;
    //let deltav: f64 = ase.channel_spacing * 0.2; //1 Angstrom = 0.1 nm filter therefore (62.5 * 2) = 1 nm / 10 page 231 c7
    let hvdeltav: f64 = h * v * deltav;
    let N2: f64 = vectors.sums[STEPS-1][0] ;
    let N1: f64 = params.N - N2;
    let sigma_e: f64 = find_emission(fASE_spectra[i] * 1.0e-9);
    let sigma_a: f64 = find_absorption(fASE_spectra[i] * 1.0e-9);
    let G: f64 = 10_f64.powf(fASE_spectral_gain[i] / 10.0);
    let nsp: f64 =  N2 / (N2 - (sigma_a/sigma_e) * N1);
    let PASE: f64 = 2.0 * nsp * hvdeltav * (G-1.0);
    fASE_power_new[i] = PASE * 1000.0;
}

//New way to calculate noise
let mut Pout: f64 = 0.0;
for i in 0..FASE {
    let N2: f64 = vectors.sums[STEPS-1][0] ;
    let N1: f64 = params.N - N2;
    let sigma_e: f64 = find_emission(fASE_spectra[i] * 1.0e-9);
    let sigma_a: f64 = find_absorption(fASE_spectra[i] * 1.0e-9);
    let nsp: f64 =  N2 / (N2 - (sigma_a/sigma_e) * N1);
    let G: f64 = 10_f64.powf(fASE_spectral_gain[i] / 10.0);
    let NF: f64 = 10.0 * ((2.0 * nsp * (G-1.0)/G) + 1.0/G).log(10.0); /*Seems to be ok,
    currently using. G>>1, 10 dB gain or higher*/
    fASE_noise_new[i] = NF;
    let v: f64 = find_freq(fASE_spectra[i] * 1.0e-9);
    let deltav: f64 = ase.channel_spacing;
    let hvdeltav: f64 = h * v * deltav;
    let Pin: f64 = vectors.signals[0][0];
    Pout = signal_for_noise[i];
    if Pout < 1.0e-9 {Pout = 1.0e-9}
    //let PASE: f64 = fASE_sum[STEPS-1];
    let PASE: f64 = 2.0 * fASE_power[i] * 1.0e-3; //New way to calc power (wrong)
    let NF: f64 = 10.0 * (((PASE /(hvdeltav)) + 1.0) * (Pin/Pout)).log(10.0);
    fASE_noise[i] = NF;
}

let mut power_old: f64 = 0.0;
let mut power_new: f64 = 0.0;
for elements in fASE_power.iter() {
    power_old += *elements
}
for elements in fASE_power_new.iter() {
    power_new += *elements
}

if val4[3] == true {
/*
/*ASE PLOT*/
let mut fg1 = Figure::new();
fg1.axes2d()
.lines(vectors.position.iter(), fASE_sum.iter(), &[Caption("Forwards ASE sum"), Color("blue")])
.lines(vectors.position.iter(), bASE_sum.iter(), &[Caption("Backwards ASE sum"), Color("green")])
.set_x_label("Z position in cm", &[])
.set_y_label("Power in mW",&[]);
//fg1.set_terminal("pngcairo", "ase_Power.png");
fg1.show();
*/
/*Forward ASE spectral noise*/
let mut fg2 = Figure::new();
fg2.axes2d()
//.lines(fASE_spectra.iter().skip(20).take(70),
//fASE_noise.iter().skip(20).take(70), &[Caption("Spectral noise, valid for gain>10 dB"),
//Color("blue")])
.lines(fASE_spectra.iter().skip(20).take(70),
fASE_noise_new.iter().skip(20).take(70), &[Caption("Spectral noise"), //Valid for G>>1, if G = 10 then gain = 7 dB
Color("green")]) //Modified to take values for the range 1530-1565, valid for 0.5 nm spacing
.set_x_label("Wavelength in nm", &[])
.set_y_label("Noise in dB",&[]);
//fg2.set_terminal("pngcairo", "ase_spectral_noise.png");
fg2.show();

/*
/*Forward ASE spectral power*/
let mut fg9 = Figure::new();
fg9.axes2d()
.lines(fASE_spectra.iter(), fASE_power.iter(), &[Caption("Spectral Forward ASE"), Color("blue")])
//.lines(fASE_spectra.iter(), fASE_power_new.iter(), &[Caption("Spectral Forward ASE new"), Color("green")])
.set_x_label("Wavelength in nm", &[])
.set_y_label("Power in mW",&[]);
//fg9.set_terminal("pngcairo", "ase_spectral_power.png");
fg9.show();
*/

/*PUMP PLOT*/
let mut fg3 = Figure::new();
fg3.axes2d()
.lines(vectors.position.iter(), fpump_sum.iter(), &[Caption("Forward Pump Power"), Color("black")])
.lines(vectors.position.iter(), bpump_sum.iter(), &[Caption("Backward Pump Power"), Color("green")])
.set_x_label("Z position in cm", &[])
.set_y_label("Power in mW",&[]);
//fg3.set_terminal("pngcairo", "pump_power.png");
fg3.show();
/*
/*LOWEST SIGNAL PLOT*/
let mut fg4 = Figure::new();
fg4.axes2d()
.lines(vectors.position.iter(), signal_sum_lowest.iter(), &[Caption("Lowest gain Signal Power"), Color("red")])
.set_x_label("Z position in cm", &[])
.set_y_label("Power in mW",&[]);
fg4.set_terminal("pngcairo", "low_signal_power.png");
fg4.show();

/*TOTAL SIGNAL PLOT*/
let mut fg8 = Figure::new();
fg8.axes2d()
.lines(vectors.position.iter(), signal_sum.iter(), &[Caption("Total Signal Power"), Color("red")])
.set_x_label("Z position in cm", &[])
.set_y_label("Power in mW",&[]);
fg8.set_terminal("pngcairo", "total_signal_power.png");
fg8.show();
*/
/*N2 PLOT*/
let mut fg5 = Figure::new();
fg5.axes2d()
.lines(vectors.position.iter(), N2_sum.iter(), &[Caption("Excited ions"), Color("yellow")])
.set_x_label("Z position in cm", &[])
.set_y_label("Fractional Population Inversion",&[]);
//fg5.set_terminal("pngcairo", "population_inversion.png");
fg5.show();

/*Lowest Gain plot*/
let mut fg6 = Figure::new();
fg6.axes2d()
.lines(vectors.position.iter(), gain.iter(), &[Caption("Lowest Gain evolution"), Color("red")])
.set_x_label("Z position in cm", &[])
.set_y_label("Gain in dB",&[]);
//fg6.set_terminal("pngcairo", "lowest_gain_evolution.png");
fg6.show();

/*Spectral Gain plot*/
let mut fg7 = Figure::new();
fg7.axes2d()
.lines(signal_spectra.iter(), spectral_gain.iter(), &[Caption("Spectral gain"), Color("red")])
.set_x_label("Wavelength in nm", &[])
.set_y_label("Gain in dB",&[]);
//fg7.set_terminal("pngcairo", "spectral_gain.png");
fg7.show();

} else {
    println!("______________________________________________");
    println!("Note: Gain Plots disabled, writing data to file gain_data.csv");
    let mut gain_data = match File::create("gain_data.csv") {
        Ok(contents) => contents,
        Err(why)     => panic!("Couldn't create gain_data.csv, {}", why.description()),
    };

    write!(gain_data, "Spectral gain,, Spectral noise,, Inverted population,,\
     Forward Pump power,,Backward Pump Power,,Lowest signal gain,,Forward ASE,,\
     Backwards ASE\n");
    write!(gain_data, "Wavelength (nm), Gain (dB), Wavelength (nm), Noise (dB),\
    Position (cm), Inverted Population (%), Position (cm), Power (mW), \
    Position(cm), Power (mW),Position (cm), Gain (dB),Position (cm), \
    Forwards ASE sum (mW), Position(cm), Backwards ASE sum (mW)\n");

    for i in 0..vectors.position.len() {
        if i < signal_spectra.len() {
            write!(gain_data,"{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
            signal_spectra[i], spectral_gain[i], fASE_spectra[i], fASE_noise_new[i],
            vectors.position[i], N2_sum[i], vectors.position[i], fpump_sum[i],
            vectors.position[i], bpump_sum[i], vectors.position[i], gain[i],
            vectors.position[i], fASE_sum[i], vectors.position[i], bASE_sum[i]
        ).map(|value| value);
        } else if i <  fASE_spectra.len() {
            write!(gain_data,",,{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n",
            fASE_spectra[i], fASE_noise_new[i],
            vectors.position[i], N2_sum[i], vectors.position[i], fpump_sum[i],
            vectors.position[i], bpump_sum[i], vectors.position[i], gain[i],
            vectors.position[i], fASE_sum[i], vectors.position[i], bASE_sum[i]
        )
            .map(|value| value);
        } else {
            write!(gain_data,",,,,{},{},{},{},{},{},{},{},{},{},{},{}\n",
            vectors.position[i], N2_sum[i], vectors.position[i], fpump_sum[i],
            vectors.position[i], bpump_sum[i], vectors.position[i], gain[i],
            vectors.position[i], fASE_sum[i], vectors.position[i], bASE_sum[i]
        )
            .map(|value| value);
        }

    }

    println!("signal_spectra[2] = {}", signal_spectra[2]);

}
//End of gain block
}

//End of main()
}
