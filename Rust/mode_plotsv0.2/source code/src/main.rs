extern crate rand;
extern crate gnuplot;

use rand::distributions::{IndependentSample, Range};
use gnuplot::*;
use std::{f64, i64};
use std::f64::consts;
use std::error::Error;
use std::str::FromStr;
use std::fs::File;
use std::fs::create_dir_all;
use std::fs::remove_dir_all;
use std::io::prelude::*;
use std::path::Path;
use std::{thread, time};



/* /////////////////////////////CONSTANTS///////////////////////////// */
const PI: f64 = consts::PI;
const h: f64 = 6.626070e-34; // Planck's constant, J.s
const c: f64 =  299792458.0; // Speed of light, m/s
const mu_0: f64 = (4.0 * PI) * 1e-7;  // Permeability of free space in H/m
const epsilon_0: f64 = 8.854187817 * 1e-12;

/* /////////////////////////////STRUCTS///////////////////////////// */
//#[derive(Clone)]
struct text_file {
    values: Vec<f64>,
    booleans: Vec<bool>,
}

struct waveguide {
    n_up: f64, //Upper cladding refractive index
    n_dn: f64, //Lower cladding refractive index
    n_l: f64, //Left cladding refractive index
    n_r: f64, //Right cladding refractive index
    n_co: f64, //Core refractive index
    a: f64, //Half of the width in meters
    b: f64, //Half of the depth in meters
    waveln: f64, //Wavelength in meters
    freq: f64, //Frequency in Hz
    w_ang: f64, //Angular frequency
    k_0: f64, //Wavenumber in vacuum in radians/meters
    k_up: f64, //Wavenumber upper cladding in radians/meters
    k_dn: f64, //Wavenumber in lower cladding in radians/meters
    k_l: f64, //Wavenumber in left cladding in radians/meters
    k_r: f64, //Wavenumber in right cladding in radians/meters
    k_co: f64, //Wavenumber in core in radians/meters
}

struct results {
    m: f64, //mode in the x direction
    n: f64, //mode in the y direction
    k_x: f64, //traverse wavenumber in the x direction in radians / meter
    k_y: f64, //traverse wavenumber in the y direction in radians / meter
    beta: f64, //Wavenumber in the medium
    n_eff: f64, //Effective refractive index
    phi_x: f64, //Phase shift to compensate for asymmetric modes
    phi_y: f64, //Phase shift to compensate for asymmetric modes
    E: Vec<f64>, //Intensity of the field (Electric??)
}



/* /////////////////////////////FUNCTIONS///////////////////////////// */
fn read_values(a: &mut text_file, s: &str) -> () {
    let path = Path::new(s);
    let display = path.display();
    let mut file = match File::open(path) {
        Ok(contents) => contents,
        Err(why)     => panic!("Couldn't open file in {}, {}",
                        display, why.description()),
    };

    let mut contents: String = String::new();
    match file.read_to_string(&mut contents) {
        Ok(_)       => println!("File read successfully."),
        Err(why)    => panic!("Unable to read file, {}", why.description()),
    }

    let text_to_parse: Vec<&str> = contents.split("\n").collect();
    for elements in text_to_parse.into_iter() {
        let intermediate: Vec<&str> = elements.split("=").collect();
        for values in intermediate.into_iter() {
            match values.trim().parse::<f64>() {
                Ok(val) => a.values.push(val),
                Err(_)  => match FromStr::from_str(values.trim()){
                    Ok(aa) => a.booleans.push(aa),
                    Err(_) => continue,
                }
            }
        }
    }
    drop(file);

} //End of read values
/******************************************************************************/
/******************************************************************************/
fn initialize_waveguide(a: &text_file, b: &mut waveguide) -> () {
    b.n_up = a.values[0];
    b.n_dn = a.values[1];
    b.n_l = a.values[2];
    b.n_r = a.values[3];
    b.n_co = a.values[4];
    b.a = a.values[5] * 1.0e-6/ 2.0;
    b.b = a.values[6] * 1.0e-6/ 2.0;
    b.waveln = a.values[7] * 1.0e-9;
    b.freq = c / b.waveln;
    b.w_ang = 2.0 * PI * b.freq;
    b.k_0 = b.w_ang * (epsilon_0 * mu_0).sqrt();
    b.k_up = b.k_0 * b.n_up;
    b.k_dn = b.k_0 * b.n_dn;
    b.k_l = b.k_0 * b.n_l;
    b.k_r = b.k_0 * b.n_r;
    b.k_co = b.k_0 * b.n_co;
} //End of initialize_waveguide
/******************************************************************************/
/******************************************************************************/
fn find_k_x(a: &waveguide, b: &mut results) -> () {
    let mut lowest_val: f64 = a.k_co;
    let mut new_val: f64 = 0.0;
    let mut try_val: f64 = 0.0;
    let c_1a: f64 = a.n_co.powi(2)/a.n_l.powi(2);
    let c_1b: f64 = a.n_co.powi(2)/a.n_r.powi(2);
    let mut c2: f64 = 0.0;
    let mut c3: f64 = 0.0;
    let mut k_x: f64 = 0.0;
    let mut upper_limit: f64 = 0.0;
    if a.k_co.powi(2) - a.k_l.powi(2) < a.k_co.powi(2) - a.k_r.powi(2) {
        upper_limit = (a.k_co.powi(2) - a.k_l.powi(2)).sqrt();
    } else {
        upper_limit = (a.k_co.powi(2) - a.k_r.powi(2)).sqrt();
    }

    for i in 1..(upper_limit as i64) {
        try_val = i as f64;
        c2 = (a.k_co.powi(2) - a.k_l.powi(2) - try_val.powi(2)).sqrt() / try_val;
        c3 = (a.k_co.powi(2) - a.k_r.powi(2) - try_val.powi(2)).sqrt() / try_val;
        new_val = ((c_1a * c2).atan() + (c_1b * c3).atan() - 2.0 * a.a * try_val + b.m * PI).abs();
        if new_val < lowest_val {
            lowest_val = new_val;
            k_x = try_val;
        } else {
            break;
        }

        if (k_x as i64) == (upper_limit as i64)-1 {
            b.k_x = -1.0;
            return
        }
    }

    b.k_x = k_x;
} //End of find_k_x

fn find_k_y(a: &waveguide, b: &mut results) -> () {
    let mut lowest_val: f64 = a.k_co;
    let mut new_val: f64 = 0.0;
    let mut try_val: f64 = 0.0;
    let mut c2: f64 = 0.0;
    let mut c3: f64 = 0.0;
    let mut k_y: f64 = 0.0;
    let mut upper_limit: f64 = 0.0;
    if a.k_co.powi(2) - a.k_up.powi(2) < a.k_dn.powi(2) - a.k_r.powi(2) {
        upper_limit = (a.k_co.powi(2) - a.k_up.powi(2)).sqrt();
    } else {
        upper_limit = (a.k_co.powi(2) - a.k_dn.powi(2)).sqrt();
    }

    for i in 1..(upper_limit as i64) {
        try_val = i as f64;
        c2 = (a.k_co.powi(2) - a.k_up.powi(2) - try_val.powi(2)).sqrt() / try_val;
        c3 = (a.k_co.powi(2) - a.k_dn.powi(2) - try_val.powi(2)).sqrt() / try_val;
        new_val = ((c2).atan() + (c3).atan() - 2.0 * a.b * try_val + b.n * PI).abs();
        if new_val < lowest_val {
            lowest_val = new_val;
            k_y = try_val;
        } else {
            break;
        }

        if (k_y as i64) == (upper_limit as i64)-1 {
            b.k_y = -1.0;
            return
        }
    }

    b.k_y = k_y;
} //End of find_k_y

fn find_beta(a: &waveguide, b: &mut results) -> () {
    let mut beta_sqr: f64 = 0.0;
    beta_sqr = a.k_co.powi(2) - b.k_x.powi(2) - b.k_y.powi(2);
    if beta_sqr < 0.0 {
        b.beta = -1.0;
    } else {
        b.beta = beta_sqr.sqrt();
        b.n_eff = b.beta / a.k_0
    }
} //End of find_beta
/******************************************************************************/
/******************************************************************************/
fn set_phase(a: &mut results) -> () {
    a.phi_x = 0.5 * a.m * PI;
    a.phi_y = 0.5 * a.n * PI;
} //End of set_phase
/******************************************************************************/
/******************************************************************************/
fn find_field(a: &waveguide, b: &mut results, x: f64, y: f64) -> () {
    let mut field: f64 = 0.0;
    let c_up: f64 = (b.beta.powi(2) - a.k_up.powi(2) + b.k_x.powi(2)).sqrt();
    let c_dn: f64 = (b.beta.powi(2) - a.k_dn.powi(2) + b.k_x.powi(2)).sqrt();
    let c_l: f64 = (b.beta.powi(2) - a.k_l.powi(2) + b.k_y.powi(2)).sqrt();
    let c_r: f64 = (b.beta.powi(2) - a.k_r.powi(2) + b.k_y.powi(2)).sqrt();
    if (x >= -1.0 * a.a && x <= a.a) && (y >= -1.0 * a.b && y <= a.b) {
        //Core
        field = (b.k_x * x + b.phi_x).cos() * (b.k_y * y + b.phi_y).cos();
        b.E.push(field);
    } else if (x > a.a) && (y >= -1.0 * a.b && y <= a.b) {
        //right
        field = (b.k_y * y + b.phi_y).cos() * (-1.0 * c_r * (x-a.a)).exp() * (b.k_x * a.a + b.phi_x).cos();
        b.E.push(field);
    } else if (x >= -1.0 * a.a && x <= a.a) && (y > a.b) {
        //up
        field = (b.k_x * x + b.phi_x).cos() * (-1.0 * c_up * (y-a.b)).exp() * (b.k_y * a.b + b.phi_y).cos();
        b.E.push(field);
    } else if (x < -1.0 * a.a) && (y >= -1.0 * a.b && y <= a.b) {
        //left
        field = (b.k_y * y + b.phi_y).cos() * (1.0 * c_l * (x+a.a)).exp() * (b.k_x * -1.0 * a.a + b.phi_x).cos();
        b.E.push(field);
    } else if (x >= -1.0 * a.a && x <= a.a) && (y < -1.0 * a.b) {
        //down
        field = (b.k_x * x + b.phi_x).cos() * (1.0 * c_dn * (y+a.b)).exp() * (b.k_y * -1.0 * a.b + b.phi_y).cos();
        b.E.push(field);
    } else {
        field = f64::NAN;
        b.E.push(field);
    }
}//End of find_field

fn populate_field(a: &waveguide, b: &mut results, points: usize) -> () {
    let mut range_x: Vec<f64> = vec![0.0; points];
    let mut range_y: Vec<f64> = vec![0.0; points];
    let mut current_step_a: f64 = 0.0;
    let mut current_step_b: f64 = 0.0;

    for i in 0..points {
        current_step_a = -2.0 * a.a + a.a * 4.0 * (i as f64) / (points as f64);
        current_step_b = -2.0 * a.b + a.b * 4.0 * (i as f64) / (points as f64);
        range_x[i] = current_step_a;
        range_y[i] = current_step_b;
    }

    for y in range_y.into_iter() {
        for x in range_x.clone().into_iter() {
            find_field(a,b,x,y);
        }
    }
}//End populate_field
/******************************************************************************/
/******************************************************************************/
fn plot(a: &waveguide, b: &results, points: usize) -> bool {
    let lowest: f64 = find_lowest(a);
    if b.k_x > 0.0 && b.k_y > 0.0 && b.beta > 0.0 && (b.n_eff > lowest && b.n_eff <= a.n_co)  {
        let mut figure = Figure::new();
        figure.axes3d()
        .set_title(&format!("EM field for m = {}, n = {}",b.m, b.n), &[])
        .surface(b.E.iter(), points, points,
        Some((-2.0 * a.a * 1e6, -2.0 * a.b * 1e6, 2.0 * a.a * 1e6, 2.0 * a.b * 1e6)), &[])
        //.set_z_label("E_x", &[])
        .set_x_label("Width", &[])
        .set_y_label("Depth",&[])
        .set_view_map();
        figure.set_terminal("pngcairo", &format!("./mode_graphs/EM_field_m_{}_n_{}.png",
                    b.m as i32, b.n as i32));
        figure.show();
        return true
    } else {
        return false
    }
}//End of plot
/******************************************************************************/
/******************************************************************************/
fn find_lowest(a: &waveguide) -> f64 {
    let mut all_n: [f64;4] = [0.0;4];
    all_n[0] = a.n_up;
    all_n[1] = a.n_dn;
    all_n[2] = a.n_l;
    all_n[3] = a.n_r;
    let mut lowest: f64 = 10000.0;
    for i in all_n.into_iter() {
        if i < &lowest {
            lowest = *i;
        }
    }
    lowest
}//End of find lowest


/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
fn main() {
    //Removing previous directories, creating new ones
    println!("Removing old directories and files...");
    remove_dir_all("mode_graphs").map_err(|err| return);
    remove_dir_all("mode_data").map_err(|err| return);

    println!("Creating new directories...");
    create_dir_all("mode_graphs").map_err(|err|
        panic!("Unable to create mode_graphs directory: {}", err.description()));
    create_dir_all("mode_data").map_err(|err|
        panic!("Unable to create mode_data directory: {}", err.description()));

    println!("Creating data structures...");
    let resolution: usize = 261; //Plot resolution, 261 recommended

    let mut text_vals: text_file = text_file {
        values: vec![],
        booleans: vec![],
    };

    let mut data: waveguide = waveguide {
        n_up: 0.0,
        n_dn: 0.0,
        n_l: 0.0,
        n_r: 0.0,
        n_co: 0.0,
        a: 0.0,
        b: 0.0,
        waveln: 0.0,
        freq: 0.0,
        w_ang: 0.0,
        k_0: 0.0,
        k_up: 0.0,
        k_dn: 0.0,
        k_l: 0.0,
        k_r: 0.0,
        k_co: 0.0,
    };

    let mut E_x: results = results {
        m: 0.0,
        n: 0.0,
        k_x: 0.0,
        k_y: 0.0,
        beta: 0.0,
        n_eff: 0.0,
        phi_x: 0.0,
        phi_y: 0.0,
        E: vec![],
    };

    println!("Creating files...");
    //Creating a file to store the output data:
    let mut file = match File::create("./mode_data/mode_data.csv") {
        Ok(contents) => contents,
        Err(e) => panic!("Unable to create file: {}", e.description()),
    };

    match write!(file,"wavelength in nm,width in micron (2a),depth in micron (2b),n_up,n_dn,n_l,n_r,n_co,n_eff,\
    TE (m),TM (n),k_x,k_y,beta\n") {
        Ok(content) => content,
        Err(why) => panic!("Unable to write mode_data.csv: {}", why.description()),
    };

    //Main body
    println!("Reading config.txt...");
    read_values(&mut text_vals, "config.txt"); //read the value from the target file and stores it in the struct
    initialize_waveguide(&text_vals, &mut data); //Input values from the text file to the struct

    println!("Calculating modes...");
    //Loop through first 100 modes
    for m in 0..10 {
        for n in 0..10{
            E_x.m = m as f64;
            E_x.n = n as f64;
            E_x.E = vec![]; //This clears the vector so it can be populated again
            find_k_x(&data, &mut E_x);
            find_k_y(&data, &mut E_x);
            find_beta(&data,&mut E_x); //Also finds n_eff
            set_phase(&mut E_x);
            populate_field(&data, &mut E_x, resolution); //Populates the vector E inside result struct
            if plot(&data, &E_x, resolution) == false {
                break;
            };//This will plot and also terminate the inner loop if the condition fails
            match write!(file,"{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n", data.waveln * 1.0e9,
            data.a * 2.0e6, data.b * 2.0e6, data.n_up, data.n_dn, data.n_l, data.n_r, data.n_co,
            E_x.n_eff, E_x.m, E_x.n, E_x.k_x as i64, E_x.k_y as i64, E_x.beta as i64) {
                Ok(contents) => contents,
                Err(why) => panic!("Unable to write mode_data.csv: {}", why.description()),
            };
        }
    }

    println!("Program completed successfully");



}//End of main()
