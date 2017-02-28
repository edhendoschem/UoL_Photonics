extern crate gnuplot;
use gnuplot::*;
use std::{f64, i64};
use std::f64::consts;

struct StartingParams {
    spot_size: f64, // m^2
    width: f64, // m
    depth: f64, // m
    length: f64, // m
    erbium_ions_per_volume: f64, // ions/m^3
    wavelength_pump: f64, // m
    wavelength_signal: f64, // m
    pump_absorption_crosssection: f64, // m^2
    signal_absorption_crosssection: f64, // m^2
    signal_emission_crosssection: f64, // m^2
    lifetime: f64, // seconds
    planck_const: f64, // J.s
    frequency_pump: f64, // 1/s
    frequency_signal: f64, // 1/s
}

const PLANCK: f64 = 6.626070e-34; //Planck's constant J.s
const C: f64 =  299792458.0; //Speed of light m/s

fn find_freq(wavelength: &f64) -> f64 {
    //Wavelength in metres
    C / wavelength
}

//Threshold intensity for gain in the signal field Watts/m^2
fn Ith(frequency: &f64, cross_section: &f64, lifetime: &f64) -> f64 {
    (PLANCK * frequency) / (cross_section * lifetime)
}

fn find_area(width: &f64, depth: &f64) -> f64 {
    width * depth
}

fn find_volume(area: &f64, length: &f64) -> f64 {
    area * length
}

//Converts dBm to Watts
fn convert_dBm(dBm: &f64) -> f64 {
    1.0 * 10_f64.powf((dBm - 30.0) / 10.0)
}

//Finds the normalized signal intensity with respect to z
fn norm_signal_z(
    norm_signal_int: &f64, // Normalized signal intensity I_signal_norm
    gain_coefficient: &f64, // Gain coefficient alpha_pump
    z_position: &f64 // position in mm
    ) -> f64 {
        // Note this is equivalent to:
        //norm_signal_int * e^(gain_coefficient * z_position)
        norm_signal_int * (gain_coefficient * (z_position / 1000.0)).exp()
}

fn n_quantity(
    pump_freq: &f64,
    pump_crosssection: &f64,
    signal_freq: &f64,
    signal_crosssection: &f64
    ) -> f64 {
        ((PLANCK * pump_freq) / (PLANCK * signal_freq))
        * (signal_crosssection / pump_crosssection)
}

//Equivalent radius: Calculate the area equivalent radius of a rectangle
fn equiv_rad(a: &f64, b: &f64) -> f64 {
    (a * b / consts::PI).sqrt()
}

//Saturation intensity
fn I_sat(norm_pump_int_z: &f64, n_quantity: &f64) -> f64 {
    (1.0 + norm_pump_int_z) / (2.0 * n_quantity)
}

//Gain calculation
fn gain(Signal_int_end: &f64, Signal_int_start: &f64) -> f64 {
    10.0 * (Signal_int_end / Signal_int_start).log(10.0)
}

//Overlap factor Gamma for gaussian profile
fn overlap_factor(spot_size: &f64, equivalent_radius: &f64) -> f64 {
    1.0 - (-1.0 * equivalent_radius.powf(2.0) / spot_size.powf(2.0)).exp()
}

//Saturation power calculation
fn P_sat(
    area: &f64,
    absorption_crosssection: &f64,
    emission_crosssection: &f64,
    lifetime: &f64,
    overlap_factor: &f64
) -> f64 {
    area /
    ((absorption_crosssection + emission_crosssection) *
    lifetime * overlap_factor)
}

//alpha_factor for power calculation
fn alpha(
    erbium_conc: &f64,
    overlap_factor: &f64,
    signal_emission_crosssection: &f64
) -> f64 {
    erbium_conc * overlap_factor * signal_emission_crosssection
}

//Convert dBm into P in terms of photons/s (1/s units) assuming 1 second
fn dBm_to_photons(dBm: &f64) -> f64 {
    let x: f64 = convert_dBm(dBm);
    x / PLANCK
}


fn test() -> f64 {
    (-1.0 * 3_f64.powf(2.0)/2_f64.powf(3.0)).exp()
}

fn main() {

    let params: StartingParams = StartingParams {
        width: 5.0e-6, // m
        depth: 5.0e-6, // m
        length: 1.0e-2, //m
        erbium_ions_per_volume: 4.45e26, // atoms/m^3
        spot_size: 1.963e-11, // m^2
        wavelength_pump: 980e-9, // m
        wavelength_signal: 1536e-9, // m
        pump_absorption_crosssection: 3.12e-25, // m^2
        signal_absorption_crosssection: 5.8e-25, // m^2
        signal_emission_crosssection: 7.27e-25, // m^2
        lifetime: 12.94e-3, // seconds
        planck_const: PLANCK, // J.s
        frequency_pump: find_freq(&980e-9), // 1/s
        frequency_signal: find_freq(&1536e-9), // 1/s
    };

    //Signal vector in dBm
    let signal_dBm: Vec<f64>  = vec![-60.0, -40.0, -20.0, -10.0, 0.0];

    //Populate the vector with values from -20 dBm to 20 dBm in 0.1 intervals
    let mut pump_dBm: Vec<f64> = vec![];
    for elements in -200..201 {
        let x: f64 = elements as f64;
        pump_dBm.push(x / 10.0);
    };

    //Threshold intensity in Watts/m^2
    let I_th: f64 = Ith(
        &params.frequency_pump,
        &params.pump_absorption_crosssection,
        &params.lifetime);

    //Area of the cross section in m^2
    let area: f64 = find_area(&params.width, &params.depth);
    let volume: f64 = find_volume(&area, &params.length);

    //Threshold power in Watts
    let threshold_power: f64 = &I_th * &area;

    //Convert signal dBm to Watts/m^2
    let mut I_signal: Vec<f64> = vec![];
    for i in signal_dBm.iter() {
        I_signal.push(
            convert_dBm(i as &f64) / &area
        )
    };

    //Convert pump dBm to watts/m^2
    let mut I_pump: Vec<f64> = vec![];
    for i in pump_dBm.iter() {
        I_pump.push(
            convert_dBm(i as &f64) / &area
        )
    };

    //Convert pump dBm to miliWatts
    let mut I_pump_mW: Vec<f64> = vec![];
    for i in pump_dBm.iter() {
        I_pump_mW.push(
            convert_dBm(i as &f64) * 1000.0
        )
    };

    //Normalized signal intensity
    let mut I_signal_norm: Vec<f64> = vec![];
    for i in I_signal.iter() {
        I_signal_norm.push(
            i / &I_th
        )
    };

    //Normalized pump intensity
    let mut I_pump_norm: Vec<f64> = vec![];
    for i in I_pump.iter() {
        I_pump_norm.push(
            i / &I_th
        )
    };

    //n_quantity determination
    let n_quantity: f64 = n_quantity(
        &params.frequency_pump,
        &params.signal_absorption_crosssection,
        &params.frequency_signal,
        &params.signal_absorption_crosssection
    );
    //Gain coefficient
    let mut alpha_pump: Vec<f64> = vec![];
    for i in I_pump_norm.iter() {
        alpha_pump.push(
        ((i - 1.0) / (i + 1.0))
        * &params.signal_emission_crosssection
        * &params.erbium_ions_per_volume
        )
    };

    //Calculating I_sat using I_pump_norm at z = 0
    let mut I_sat_z0: Vec<f64> = vec![];
    for elements in I_pump_norm.iter() {
        I_sat_z0.push(
            I_sat(elements, &n_quantity)
        )
    };

    //Z is the length in milimetres
    let mut z: Vec<f64> = vec![];
    for i in 0..11 {
        z.push(i as f64)
    };

    //Equivalent Radius in metres
    let equivalent_radius: f64 = equiv_rad(
            &params.width,
            &params.depth
    );

    //Overlap Factor, adimensional units
    let gamma: f64 = overlap_factor(
        &params.spot_size,
        &equivalent_radius
    );

    /* -60 dBm signal */
    //storage of final values
    let mut I_result_05: Vec<f64> = vec![];
    let mut I_result_1: Vec<f64> = vec![];
    let mut I_result_2: Vec<f64> = vec![];
    //storage of gain values
    let mut final_gain_05: Vec<f64> = vec![];
    let mut final_gain_1: Vec<f64> = vec![];
    let mut final_gain_2: Vec<f64> = vec![];
    //Signal that is going to be used (in this case -40 dBm)
    let signal_60: f64 = I_signal_norm[0]; // -60 dBm
    let signal_40: f64 = I_signal_norm[1]; // -40 dBm
    let signal_20: f64 = I_signal_norm[2]; // -20 dBm
    let signal_10: f64 = I_signal_norm[3]; // -10 dBm
    let signal_00: f64 = I_signal_norm[4]; //  00 dBm

    //println!("{:?}", alpha_pump);
    //Final values of I_signal_norm
    for elements in alpha_pump.iter() {
        I_result_05.push(norm_signal_z(&signal_40, elements, &5.0));
        I_result_1.push(norm_signal_z(&signal_40, elements, &10.0));
        I_result_2.push(norm_signal_z(&signal_40, elements, &20.0));
    };

    //Gain of 0.5 cm waveguide
    for elements in I_result_05.iter() {
        final_gain_05.push(gain(elements, &signal_40))
    };
    //Gain of 1 cm waveguide
    for elements in I_result_1.iter() {
        final_gain_1.push(gain(elements, &signal_40))
    };
    //Gain of 2 cm waveguide
    for elements in I_result_2.iter() {
        final_gain_2.push(gain(elements, &signal_40))
    };

    //Plot of gain vs power
    //pump power in mW
    let ref x = I_pump_mW;
    //Pump power in dBm
    //let ref x = pump_dBm;
    let ref y_05 = final_gain_05;
    let ref y_1 = final_gain_1;
    let ref y_2 = final_gain_2;
    let mut fg = Figure::new();
    fg.axes2d()
    .lines(x, y_2, &[Caption("2.0 cm waveguide"), Color("blue")])
    .lines(x, y_1, &[Caption("1.0 cm waveguide"), Color("green")])
    .lines(x, y_05, &[Caption("0.5 cm waveguide"), Color("red")])
    .set_x_label("Pump power in mW", &[])
    .set_y_label("Gain in dB",&[]);
    fg.show();

    /*Population inversion*/
    let mut phi_p: Vec<f64> = vec![];
    for elements in I_pump.iter() {
        phi_p.push( elements / (PLANCK * &params.frequency_pump ))
    }
    let phi_th: f64 = 1.0 / (&params.lifetime * &params.pump_absorption_crosssection);

    let mut phi_p_norm: Vec<f64> = vec![];
    for elements in phi_p.iter() {
        phi_p_norm.push(elements / &phi_th)
    }

    let mut frac_inv: Vec<f64> = vec![];
    for elements in phi_p_norm.iter() {
        frac_inv.push((elements - 1.0) / (elements + 1.0))
    }

    let mut zero_vector: Vec<f64> = vec![0.0; x.len()];
    let mut fg2 = Figure::new();
    fg2.axes2d()
    .lines(x, frac_inv, &[Caption("Fractional Inversion"), Color("blue")])
    .lines(x, zero_vector, &[Caption("Inversion Threshold"), Color("black")])
    .set_x_label("Pump power in mW", &[])
    .set_y_label("Fractional inversion",&[]);
    fg2.show();
    println!("Threshold power = {}", threshold_power * 1000.0);
    println!("Frequency pump {}, Frequency signal {}", params.frequency_pump, params.frequency_signal);
}
