//aquiaqui page 30 of Okamoto book, pending plot electric fields instead of magnetic
//Pendiente: 1) Eliminar Ez y d1Hx1, 2)Aplicar Ey a todo, 3) Examinar porque Ey se ve raro
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


/////////////////////////////CONSTANTS/////////////////////////////
const PI: f64 = consts::PI;
const h: f64 = 6.626070e-34; // Planck's constant, J.s
const c: f64 =  299792458.0; // Speed of light, m/s
const mu_0: f64 = (4.0 * PI) * 1e-7;                                            // Permeability of free space in H/m
const epsilon_0: f64 = 8.854187817 * 1e-12;
/////////////////////////////STRUCTS/////////////////////////////
struct mode_params {
    n_0: f64,                                                                   //Cladding index
    n_1: f64,                                                                   //Core index
    a: f64,                                                                     //Width
    d: f64,                                                                     //Depth
    wavelength: f64,                                                            //Wavelength
    P: usize,                                                                   //Modenumber in X direction
    Q: usize,                                                                   //Modenumber in Y direction
    STEPS: i64,                                                                 //Number of steps to be used
    frequency: f64,                                                             //Frequency
    w_ang: f64,                                                                 //Angular frequency
    k: f64,                                                                     //Wavenumber in vacuum
}

struct mode_solutions {
    sols_x: Vec<f64>,
    sols_y: Vec<f64>,
    E_x: Vec<f64>,
    E_y: Vec<f64>,
}
/////////////////////////////FUNCTIONS/////////////////////////////
fn find_ang_freq(frequency: f64) -> f64 {
    2.0 * PI * frequency
}

fn find_freq(wavelength: f64) -> f64 {
        c / wavelength
}
/////////////////////////////MODS/////////////////////////////
mod E_x {
    extern crate rand;
    use rand::thread_rng;
    use rand::distributions::{IndependentSample, Range};
    use std::{f64, i64};

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

    //////////////////////////Region 1, -a<= x <= a, -d<=y<=d////////////////////////////
    pub fn Hy1(A: f64, k_x: f64, k_y: f64, phi: f64, psi: f64, x: f64, y: f64)
    -> f64 {
        A * (k_x * x - phi).cos() * (k_y * y - psi).cos()
    }

    pub fn d2Hy1(A: f64, k_x: f64, k_y: f64, phi: f64, psi: f64, x: f64, y: f64)
    -> f64 {
        -1.0 * A * (k_x).powi(2) * (k_x * x - phi).cos() * (k_y * y - psi).cos()
    }

    pub fn Ex1(Hy1: f64, d2Hy1: f64, w_ang: f64, mu_0: f64, beta: f64, epsilon_0: f64, n: f64)
    -> f64 {
        (w_ang * mu_0 / beta) * Hy1 + (1.0 / (w_ang * epsilon_0 * (n).powi(2) * beta)) * d2Hy1
    }

    ///////////////////////////////Region 2, a> x, -d<=y<=d/////////////////////////////
    pub fn Hy2(A: f64, k_x: f64, k_y: f64, gamma_x: f64, phi: f64, psi: f64,
        x: f64, y: f64, a: f64) -> f64 {
        A * (k_x * a - phi).cos() * (-1.0 * gamma_x * (x - a)).exp()
        * (k_y * y - psi).cos()
    }

    pub fn d2Hy2(A: f64, k_x: f64, k_y: f64, gamma_x: f64, phi: f64, psi: f64,
        x: f64, y: f64, a: f64) -> f64 {
        A * (gamma_x).powi(2) * (k_x * a - phi).cos() * (-1.0 * gamma_x * (x - a)).exp()
        * (k_y * y - psi).cos()
    }

    pub fn Ex2(Hy2: f64, d2Hy2: f64, w_ang: f64, mu_0: f64, beta: f64, epsilon_0: f64, n: f64)
    -> f64 {
        (w_ang * mu_0 / beta) * Hy2 + (1.0 / (w_ang * epsilon_0 * (n).powi(2) * beta)) * d2Hy2
    }

    ///////////////////////////////Region 3, -a<=x<=a, y>d/////////////////////////////
    pub fn Hy3(A: f64, k_x: f64, k_y: f64, gamma_y: f64, phi: f64, psi: f64,
        x: f64, y: f64, d: f64) -> f64 {
        A * (k_x * x - phi).cos() * (-1.0 * gamma_y * (y - d)).exp()
        * (k_y * d - psi).cos()
    }

    pub fn d2Hy3(A: f64, k_x: f64, k_y: f64, gamma_y: f64, phi: f64, psi: f64,
        x: f64, y: f64, d: f64) -> f64 {
        -1.0 * (k_x).powi(2) * A * (k_x * x - phi).cos() * (-1.0 * gamma_y * (y - d)).exp()
        * (k_y * d - psi).cos()
    }

    pub fn Ex3(Hy3: f64, d2Hy3: f64, w_ang: f64, mu_0: f64, beta: f64, epsilon_0: f64, n: f64)
    -> f64 {
        (w_ang * mu_0 / beta) * Hy3 + (1.0 / (w_ang * epsilon_0 * (n).powi(2) * beta)) * d2Hy3
    }

    ///////////////////////////////Region 4,  x<-a, -d<=y<=d/////////////////////////////
    pub fn Hy4(A: f64, k_x: f64, k_y: f64, gamma_x: f64, phi: f64, psi: f64,
        x: f64, y: f64, a: f64) -> f64 {
        A * (k_x * a + phi).cos() * (1.0 * gamma_x * (x + a)).exp()
        * (k_y * y - psi).cos()
    }

    pub fn d2Hy4(A: f64, k_x: f64, k_y: f64, gamma_x: f64, phi: f64, psi: f64,
        x: f64, y: f64, a: f64) -> f64 {
        A * (gamma_x).powi(2) * (k_x * a + phi).cos() * (1.0 * gamma_x * (x + a)).exp()
        * (k_y * y - psi).cos()
    }

    pub fn Ex4(Hy4: f64, d2Hy4: f64, w_ang: f64, mu_0: f64, beta: f64, epsilon_0: f64, n: f64)
    -> f64 {
        (w_ang * mu_0 / beta) * Hy4 + (1.0 / (w_ang * epsilon_0 * (n).powi(2) * beta)) * d2Hy4
    }

    ///////////////////////////////Region 5, -a<=x<=a, y<-d/////////////////////////////
    pub fn Hy5(A: f64, k_x: f64, k_y: f64, gamma_y: f64, phi: f64, psi: f64,
        x: f64, y: f64, d: f64) -> f64 {
        A * (k_x * x - phi).cos() * (1.0 * gamma_y * (y + d)).exp()
        * (k_y * d + psi).cos()
    }

    pub fn d2Hy5(A: f64, k_x: f64, k_y: f64, gamma_y: f64, phi: f64, psi: f64,
        x: f64, y: f64, d: f64) -> f64 {
        -1.0 * A * (k_x).powi(2) * (k_x * x - phi).cos() * (1.0 * gamma_y * (y + d)).exp()
        * (k_y * d + psi).cos()
    }

    pub fn Ex5(Hy5: f64, d2Hy5: f64, w_ang: f64, mu_0: f64, beta: f64, epsilon_0: f64, n: f64)
    -> f64 {
        (w_ang * mu_0 / beta) * Hy5 + (1.0 / (w_ang * epsilon_0 * (n).powi(2) * beta)) * d2Hy5
    }

    pub fn find_containment(A: f64, a: f64, d: f64, w_ang: f64, beta: f64,
        k_x: f64, gamma_x: f64, k_y: f64, gamma_y: f64, mu_0: f64, phi: f64,
        psi: f64, iterations: i64) -> f64 {
        let c1: f64 = (w_ang * mu_0 * beta.powi(2) - k_x.powi(2)) / beta.powi(3);
        let c2: f64 = (w_ang * mu_0 * beta.powi(2) + gamma_x.powi(2)) / beta.powi(3);
        let c3: f64 = c1.clone();
        //Region 1, -a<=x<=a, -d<=y<=d
        let range_x = Range::new(-a, a);
        let range_y = Range::new(-d, d);
        let mut rng = rand::thread_rng();
        let mut x: f64 = 0.0;
        let mut y: f64 = 0.0;
        let mut sum_x: f64 = 0.0;
        let mut sum_y: f64 = 0.0;
        let mut total: f64 = 0.0;
        let averaging: i64 = 10;
        for i in 0..averaging {
            for j in 0..iterations {
                y = range_y.ind_sample(&mut rng);
                for k in 0..iterations {
                    x = range_x.ind_sample(&mut rng);
                    sum_x += Hy1(A, k_x, k_y, phi, psi, x, y).powi(2)
                }
                sum_y += (a - (-a)) * sum_x / (iterations as f64);
                sum_x = 0.0;
            }
            total += (d - (-d)) * sum_y / (iterations as f64);
            sum_y = 0.0;
        }
        let P1: f64 = c1 * total / (averaging as f64);

        //Region 2, a<=x<=3a, -d<=y<=d
        let range_x = Range::new(a, 3.0 * a);
        let range_y = Range::new(-d, d);
        x = 0.0;
        y = 0.0;
        sum_x = 0.0;
        sum_y = 0.0;
        total = 0.0;
        for i in 0..averaging {
            for j in 0..iterations {
                y = range_y.ind_sample(&mut rng);
                for k in 0..iterations {
                    x = range_x.ind_sample(&mut rng);
                    sum_x += Hy2(A, k_x, k_y, gamma_x, phi, psi, x, y, a).powi(2)
                }
                sum_y += (3.0 * a - a) * sum_x / (iterations as f64);
                sum_x = 0.0;
            }
            total += (d - (-d)) * sum_y / (iterations as f64);
            sum_y = 0.0;
        }
        let P2: f64 = c2 * total / (averaging as f64);

        //Region 3, -a<=x<=a, d<=y<=3d
        let range_x = Range::new(-a, a);
        let range_y = Range::new(d, 3.0 * d);
        x = 0.0;
        y = 0.0;
        sum_x = 0.0;
        sum_y = 0.0;
        total = 0.0;
        for i in 0..averaging {
            for j in 0..iterations {
                y = range_y.ind_sample(&mut rng);
                for k in 0..iterations {
                    x = range_x.ind_sample(&mut rng);
                    sum_x += Hy3(A, k_x, k_y, gamma_x, phi, psi, x, y, d).powi(2)
                }
                sum_y += (3.0 * a - a) * sum_x / (iterations as f64);
                sum_x = 0.0;
            }
            total += (d - (-d)) * sum_y / (iterations as f64);
            sum_y = 0.0;
        }
        let P3: f64 = c3 * total / (averaging as f64);
        P1 / (P1 + P2 + P3)
    //End of the function
    }
}

mod E_y {
    extern crate rand;
    use rand::thread_rng;
    use rand::distributions::{IndependentSample, Range};
    use std::{f64, i64};

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

    //////////////////////////Region 1, -a<= x <= a, -d<=y<=d///////////////////////////
    pub fn Hx1(A: f64, k_x: f64, k_y: f64, phi: f64, psi: f64, x: f64, y: f64)
    -> f64 {
        A * (k_x * x - phi).cos() * (k_y * y - psi).cos()
    }

    pub fn d2Hx1(Hx1: f64, k_y: f64) -> f64 {
        -1.0 * k_y.powi(2) * Hx1
    }

    pub fn Ey1(Hx1: f64, d2Hx1: f64, w_ang: f64, mu_0: f64, beta: f64, epsilon_0: f64, n: f64)
    -> f64 {
        -1.0 * ((w_ang * mu_0 / beta) * Hx1 + (1.0 / (w_ang * epsilon_0 * (n).powi(2) * beta)) * d2Hx1)
    }

    ///////////////////////////////Region 2, a> x, -d<=y<=d/////////////////////////////
    pub fn Hx2(A: f64, k_x: f64, k_y: f64, gamma_x: f64, phi: f64, psi: f64,
        x: f64, y: f64, a: f64) -> f64 {
        A * (k_x * a - phi).cos() * (-1.0 * gamma_x * (x - a)).exp()
        * (k_y * y - psi).cos()
    }

    pub fn d2Hx2(Hx2: f64, k_y: f64) -> f64 {
        -1.0 * k_y.powi(2) * Hx2
    }

    pub fn Ey2(Hx2: f64, d2Hx2: f64, w_ang: f64, mu_0: f64, beta: f64, epsilon_0: f64, n: f64)
    -> f64 {
        -1.0 * ((w_ang * mu_0 / beta) * Hx2 + (1.0 / (w_ang * epsilon_0 * (n).powi(2) * beta)) * d2Hx2)
    }

    ///////////////////////////////Region 3, -a<=x<=a, y>d//////////////////////////////
    pub fn Hx3(A: f64, k_x: f64, k_y: f64, gamma_y: f64, phi: f64, psi: f64,
        x: f64, y: f64, d: f64) -> f64 {
        A * (k_x * x - phi).cos() * (-1.0 * gamma_y * (y - d)).exp()
        * (k_y * d - psi).cos()
    }

   pub fn d2Hx3(Hx3: f64, gamma_y: f64) -> f64 {
        1.0 * gamma_y.powi(2) * Hx3
    }

    pub fn Ey3(Hx3: f64, d2Hx3: f64, w_ang: f64, mu_0: f64, beta: f64, epsilon_0: f64, n: f64)
    -> f64 {
        -1.0 * ((w_ang * mu_0 / beta) * Hx3 + (1.0 / (w_ang * epsilon_0 * (n).powi(2) * beta)) * d2Hx3)
    }

    ///////////////////////////////Region 4,  x<-a, -d<=y<=d////////////////////////////
    pub fn Hx4(A: f64, k_x: f64, k_y: f64, gamma_x: f64, phi: f64, psi: f64,
        x: f64, y: f64, a: f64) -> f64 {
        A * (k_x * a + phi).cos() * (1.0 * gamma_x * (x + a)).exp()
        * (k_y * y - psi).cos()
    }

    pub fn d2Hx4(Hx4: f64, k_y: f64) -> f64 {
        -1.0 * k_y.powi(2) * Hx4
    }

    pub fn Ey4(Hx4: f64, d2Hx4: f64, w_ang: f64, mu_0: f64, beta: f64, epsilon_0: f64, n: f64)
    -> f64 {
        -1.0 * ((w_ang * mu_0 / beta) * Hx4 + (1.0 / (w_ang * epsilon_0 * (n).powi(2) * beta)) * d2Hx4)
    }

    ///////////////////////////////Region 5, -a<=x<=a, y<-d/////////////////////////////
    pub fn Hx5(A: f64, k_x: f64, k_y: f64, gamma_y: f64, phi: f64, psi: f64,
        x: f64, y: f64, d: f64) -> f64 {
        A * (k_x * x - phi).cos() * (1.0 * gamma_y * (y + d)).exp()
        * (k_y * d + psi).cos()
    }

    pub fn d2Hx5(Hx5: f64, gamma_y: f64) -> f64 {
         1.0 * gamma_y.powi(2) * Hx5
     }

     pub fn Ey5(Hx5: f64, d2Hx5: f64, w_ang: f64, mu_0: f64, beta: f64, epsilon_0: f64, n: f64)
     -> f64 {
         -1.0 * ((w_ang * mu_0 / beta) * Hx5 + (1.0 / (w_ang * epsilon_0 * (n).powi(2) * beta)) * d2Hx5)
     }

    pub fn find_containment(A: f64, a: f64, d: f64, w_ang: f64, beta: f64,
        k_x: f64, gamma_x: f64, k_y: f64, gamma_y: f64, mu_0: f64, phi: f64,
        psi: f64, iterations: i64) -> f64 {
        let c1: f64 = (w_ang * mu_0 * beta.powi(2) - k_y.powi(2)) / beta.powi(3);
        let c2: f64 = c1.clone();
        let c3: f64 = (w_ang * mu_0 * beta.powi(2) + gamma_y.powi(2)) / beta.powi(3);
        //Region 1, -a<=x<=a, -d<=y<=d
        let range_x = Range::new(-a, a);
        let range_y = Range::new(-d, d);
        let mut rng = rand::thread_rng();
        let mut x: f64 = 0.0;
        let mut y: f64 = 0.0;
        let mut sum_x: f64 = 0.0;
        let mut sum_y: f64 = 0.0;
        let mut total: f64 = 0.0;
        let averaging: i64 = 10;
        for i in 0..averaging {
            for j in 0..iterations {
                y = range_y.ind_sample(&mut rng);
                for k in 0..iterations {
                    x = range_x.ind_sample(&mut rng);
                    sum_x += Hx1(A, k_x, k_y, phi, psi, x, y).powi(2)
                }
                sum_y += (a - (-a)) * sum_x / (iterations as f64);
                sum_x = 0.0;
            }
            total += (d - (-d)) * sum_y / (iterations as f64);
            sum_y = 0.0;
        }
        let P1: f64 = c1 * total / (averaging as f64);

        //Region 2, a<=x<=3a, -d<=y<=d
        let range_x = Range::new(a, 3.0 * a);
        let range_y = Range::new(-d, d);
        x = 0.0;
        y = 0.0;
        sum_x = 0.0;
        sum_y = 0.0;
        total = 0.0;
        for i in 0..averaging {
            for j in 0..iterations {
                y = range_y.ind_sample(&mut rng);
                for k in 0..iterations {
                    x = range_x.ind_sample(&mut rng);
                    sum_x += Hx2(A, k_x, k_y, gamma_x, phi, psi, x, y, a).powi(2)
                }
                sum_y += (3.0 * a - a) * sum_x / (iterations as f64);
                sum_x = 0.0;
            }
            total += (d - (-d)) * sum_y / (iterations as f64);
            sum_y = 0.0;
        }
        let P2: f64 = c2 * total / (averaging as f64);

        //Region 3, -a<=x<=a, d<=y<=3d
        let range_x = Range::new(-a, a);
        let range_y = Range::new(d, 3.0 * d);
        x = 0.0;
        y = 0.0;
        sum_x = 0.0;
        sum_y = 0.0;
        total = 0.0;
        for i in 0..averaging {
            for j in 0..iterations {
                y = range_y.ind_sample(&mut rng);
                for k in 0..iterations {
                    x = range_x.ind_sample(&mut rng);
                    sum_x += Hx3(A, k_x, k_y, gamma_x, phi, psi, x, y, d).powi(2)
                }
                sum_y += (3.0 * a - a) * sum_x / (iterations as f64);
                sum_x = 0.0;
            }
            total += (d - (-d)) * sum_y / (iterations as f64);
            sum_y = 0.0;
        }
        let P3: f64 = c3 * total / (averaging as f64);
        P1 / (P1 + P2 + P3)
    //End of the function
    }
}

fn main() {
    //////////////////////////Reading config file///////////////
    let path = Path::new("config.txt");
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
    let mut data: Vec<f64> = vec![];
    let mut booleans: Vec<bool> = vec![];
    for elements in text_to_parse.into_iter() {
        let mut intermediate: Vec<&str> = elements.split("=").collect();
        for values in intermediate.into_iter() {
            match values.trim().parse::<f64>() {
                Ok(val) => data.push(val),
                Err(_)  => match FromStr::from_str(values.trim()){
                    Ok(aa) => booleans.push(aa),
                    Err(_) => continue,
                }
            }
        }
    }

    println!("data = {:?}", data);
    println!("booleans = {:?}", booleans);

    let mut params: mode_params = mode_params {
        n_0: 0.0,
        n_1: 0.0,
        a: 0.0,
        d: 0.0,
        wavelength: 0.0,
        P: 0,
        Q: 0,
        STEPS: 0,
        frequency: 0.0,
        w_ang: 0.0,
        k: 0.0,
    };

    params.n_0 = data[0];
    params.n_1 = data[1];
    params.a = data[2] * 1e-6 * 0.5;
    params.d = data[3] * 1e-6 * 0.5;
    params.wavelength = data[4] * 1.0e-9;
    params.P = 7;//data[5] as usize;
    params.Q = 7;//data[6] as usize;
    params.STEPS = 1000000;
    params.frequency = find_freq((data[4] * 1.0e-9));
    params.w_ang = find_ang_freq(params.frequency);
    params.k = params.w_ang / c;
    drop(data);

    let mut solutions: mode_solutions = mode_solutions {
        sols_x: vec![0.0; params.P * 4],
        sols_y: vec![0.0; params.Q * 4],
//0 = P or Q, 1 = k_x or k_y, 2 = gamma_x or gamma_y
        E_x: vec![0.0; params.P * params.Q * 13],
        E_y: vec![0.0; params.P * params.Q * 13],
//0=P, 1=Q, 2=a, 3=d, 4=n_0, 5=n_1, 6=n_eff, 7=k_x, 8=gamma_x, 9=k_y, 10=gamma_y, 11=beta, 12=containment
    };
    let mut comparator: f64 = 0.0;
    let mut newlow: f64 = 1.0e20;

    let mut file = match File::create("Mode_results.csv") {
        Ok(contents) => contents,
        Err(why)     => panic!("Unable to create file: {}", why.description()),
    };

    //Creating directories
    remove_dir_all("mode_graphs");
    remove_dir_all("mode_plot_data");
    create_dir_all("mode_graphs");
    create_dir_all("mode_plot_data");



//Start of calculation block
{
    if booleans[0] == true {
        let n_0: f64 = params.n_0;
        let P: usize = 2;
        let Q: usize = 2;
        let a: f64 = params.a;
        let d: f64 = params.d;
        let wavelength: f64 = params.wavelength;
        let frequency: f64 = params.frequency;
        let k: f64 = params.k;
        let STEPS: i64 = params.STEPS;
        let mut index: usize = 0;
        let mut min_n_1: [f64; 2] = [0.0,0.0];
        let mut solutions: mode_solutions = mode_solutions {
            sols_x: vec![0.0; 2 * 4],
            sols_y: vec![0.0; 2 * 4],
    //0 = P or Q, 1 = k_x or k_y, 2 = gamma_x or gamma_y
            E_x: vec![0.0; 2 * 2 * 13],
            E_y: vec![0.0; 2 * 2 * 13],
    //0=P, 1=Q, 2=a, 3=d, 4=n_0, 5=n_1, 6=n_eff, 7=k_x, 8=gamma_x, 9=k_y, 10=gamma_y, 11=beta, 12=containment
        };


        'outer: for i in 0..20000 {
            let mut n_1: f64 = params.n_0 + (i as f64) / 10000.0;
            println!("Evaluating n_1 = {:.4}", n_1);
            let STEPSIZE: f64 = k * n_1 / (STEPS as f64);
            println!("STEPSIZE = {}", STEPSIZE);
            //EX calculation//
            for p in 0..P {
                newlow = 1.0e20;
                for j in 0..STEPS {
                    let k_x: f64 = (j as f64) * STEPSIZE;
                    let phi: f64 = E_x::phi((p+1) as f64, PI);
                    let gamma_x: f64 = E_x::gamma_x(n_0, n_1, k, k_x);
                    comparator = E_x::find_k_x(n_0, n_1, (p+1) as f64, a, PI, gamma_x, k_x);
                    if comparator < 1.0e-3 && comparator < newlow && k_x > 10.0  {
                        solutions.sols_x[0 + 4 * p] = (p+1) as f64;
                        solutions.sols_x[1 + 4 * p] = k_x;
                        solutions.sols_x[2 + 4 * p] = gamma_x;
                        solutions.sols_x[3 + 4 * p] = comparator;
                        newlow = comparator;
                    }
                }
            }

            for q in 0..Q {
                newlow = 1e20;
                for j in 0..STEPS {
                    let k_y: f64 = (j as f64) * STEPSIZE;
                    let psi: f64 = E_x::psi((q+1) as f64, PI);
                    let gamma_y: f64 = E_x::gamma_y(n_0, n_1, k, k_y);
                    comparator = E_x::find_k_y(n_0, n_1, (q+1) as f64, d, PI, gamma_y, k_y);
                    if comparator < 1.0e-3 && comparator < newlow && k_y > 10.0  {
                        solutions.sols_y[0 + 4 * q] = (q+1) as f64;
                        solutions.sols_y[1 + 4 * q] = k_y;
                        solutions.sols_y[2 + 4 * q] = gamma_y;
                        solutions.sols_y[3 + 4 * q] = comparator;
                        newlow = comparator;
                    }
                }
            }

            index = 0;
            for p in 0..P {
                for q in 0..Q {
                    let k_x: f64 = solutions.sols_x[1 + 4 * p];
                    let gamma_x: f64 = solutions.sols_x[2 + 4 * p];
                    let k_y: f64 = solutions.sols_y[1 + 4 * q];
                    let gamma_y: f64 = solutions.sols_y[2 + 4 * q];
                    let beta: f64 = E_x::find_beta(n_0, n_1, k, k_x, k_y);
                    let n_eff: f64 = E_x::find_n_eff(beta, k);
                    let comparator_2: f64 = ((n_eff / n_0) - 1.0) * 100.0;
                    if (n_eff > n_0 && n_eff < n_1 /*&& comparator_2 > 0.1*/) && (k_x > 10.0 && k_y > 10.0) {
                        if min_n_1[0] < 0.1 {
                            min_n_1[0] = n_1;
                        }
                        solutions.E_x[0 + index * 13] = (p+1) as f64;
                        solutions.E_x[1 + index * 13] = (q+1) as f64;
                        solutions.E_x[2 + index * 13] = a;
                        solutions.E_x[3 + index * 13] = d;
                        solutions.E_x[4 + index * 13] = n_0;
                        solutions.E_x[5 + index * 13] = n_1;
                        solutions.E_x[6 + index * 13] = n_eff;
                        solutions.E_x[7 + index * 13] = k_x;
                        solutions.E_x[8 + index * 13] = gamma_x;
                        solutions.E_x[9 + index * 13] = k_y;
                        solutions.E_x[10 + index * 13] = gamma_y;
                        solutions.E_x[11 + index * 13] = beta;
                        solutions.E_x[12 + index * 13] = 0.0;
                        index += 1;
                        }
                    }
                }

            //Resetting vectors
            solutions.sols_x = vec![0.0; 2 * 4];
            solutions.sols_y = vec![0.0; 2 * 4];
            //EY calculation//
            for p in 0..P {
                newlow = 1.0e20;
                for j in 0..STEPS {
                    let k_x: f64 = (j as f64) * STEPSIZE;
                    let phi: f64 = E_y::phi((p+1) as f64, PI);
                    let gamma_x: f64 = E_y::gamma_x(n_0, n_1, k, k_x);
                    comparator = E_y::find_k_x(n_0, n_1, (p+1) as f64, a, PI, gamma_x, k_x);
                    if comparator < 1.0e-3 && comparator < newlow && k_x > 10.0  {
                        solutions.sols_x[0 + 4 * p] = (p+1) as f64;
                        solutions.sols_x[1 + 4 * p] = k_x;
                        solutions.sols_x[2 + 4 * p] = gamma_x;
                        solutions.sols_x[3 + 4 * p] = comparator;
                        newlow = comparator;
                    }
                }
            }

            for q in 0..Q {
                newlow = 1e20;
                for j in 0..STEPS {
                    let k_y: f64 = (j as f64) * STEPSIZE;
                    let psi: f64 = E_y::psi((q+1) as f64, PI);
                    let gamma_y: f64 = E_y::gamma_y(n_0, n_1, k, k_y);
                    comparator = E_y::find_k_y(n_0, n_1, (q+1) as f64, d, PI, gamma_y, k_y);
                    if comparator < 1.0e-3 && comparator < newlow && k_y > 10.0  {
                        solutions.sols_y[0 + 4 * q] = (q+1) as f64;
                        solutions.sols_y[1 + 4 * q] = k_y;
                        solutions.sols_y[2 + 4 * q] = gamma_y;
                        solutions.sols_y[3 + 4 * q] = comparator;
                        newlow = comparator;
                    }
                }
            }

            index = 0;
            for p in 0..P {
                for q in 0..Q {
                    let k_x: f64 = solutions.sols_x[1 + 4 * p];
                    let gamma_x: f64 = solutions.sols_x[2 + 4 * p];
                    let k_y: f64 = solutions.sols_y[1 + 4 * q];
                    let gamma_y: f64 = solutions.sols_y[2 + 4 * q];
                    let beta: f64 = E_y::find_beta(n_0, n_1, k, k_x, k_y);
                    let n_eff: f64 = E_y::find_n_eff(beta, k);
                    let comparator_2: f64 = ((n_eff / n_0) - 1.0) * 100.0;
                    if (n_eff > n_0 && n_eff < n_1 /*&& comparator_2 > 0.1*/) && (k_x > 10.0 && k_y > 10.0) {
                        if min_n_1[1] < 0.1 {
                            min_n_1[1] = n_1;
                            let mut n1comp: f64 = 1000.0;
                            if min_n_1[1] > min_n_1[0]{
                                min_n_1[0] = min_n_1[1];
                            } else {
                                min_n_1[1] = min_n_1[0];
                            }
                        }
                        solutions.E_y[0 + index * 13] = (p+1) as f64;
                        solutions.E_y[1 + index * 13] = (q+1) as f64;
                        solutions.E_y[2 + index * 13] = a;
                        solutions.E_y[3 + index * 13] = d;
                        solutions.E_y[4 + index * 13] = n_0;
                        solutions.E_y[5 + index * 13] = n_1;
                        solutions.E_y[6 + index * 13] = n_eff;
                        solutions.E_y[7 + index * 13] = k_x;
                        solutions.E_y[8 + index * 13] = gamma_x;
                        solutions.E_y[9 + index * 13] = k_y;
                        solutions.E_y[10 + index * 13] = gamma_y;
                        solutions.E_y[11 + index * 13] = beta;
                        solutions.E_y[12 + index * 13] = 0.0;
                        index += 1;
                        }
                    }
                }


            if solutions.E_y[13] > 0.1 && solutions.E_x[13] > 0.1 {
                params.n_1 = n_1 - 1.0 / 10000.0;
                solutions.E_y[5] = params.n_1;
                solutions.E_x[5] = params.n_1;
                println!("Min n_1 for a guided mode = {:.4}", min_n_1[0]);
                println!("Max n_1 before multimode = {:.4}", params.n_1);
                break 'outer;
            }
        //End of for i in 0..20000 loop
        }

        //Containment for E_y

        let mut p: f64 = solutions.E_y[0];
        let mut q: f64 = solutions.E_y[1];
        let mut a: f64 = solutions.E_y[2];
        let mut d: f64 = solutions.E_y[3];
        let mut k_x: f64 = solutions.E_y[7];
        let mut gamma_x: f64 = solutions.E_y[8];
        let mut k_y: f64 = solutions.E_y[9];
        let mut gamma_y: f64 = solutions.E_y[10];
        let mut beta: f64 = solutions.E_y[11];
        let w_ang: f64 = params.w_ang;
        let mut phi: f64 = E_y::phi(p, PI);
        let mut psi: f64 = E_y::psi(q, PI);
        let iterations: i64 = (STEPS as f64).sqrt() as i64;
        solutions.E_y[12] = E_y::find_containment(1.0, a, d, w_ang, beta, k_x,
            gamma_x, k_y, gamma_y, mu_0, phi, psi, iterations);
        println!("Solutions.E_y[0..13] = {:?}", &solutions.E_y[0..13]);

        //Containment for E_x
        p = solutions.E_x[0];
        q = solutions.E_x[1];
        a = solutions.E_x[2];
        d = solutions.E_x[3];
        k_x = solutions.E_x[7];
        gamma_x = solutions.E_x[8];
        k_y = solutions.E_x[9];
        gamma_y = solutions.E_x[10];
        beta = solutions.E_x[11];
        phi = E_x::phi(p, PI);
        psi = E_x::psi(q, PI);
        solutions.E_x[12] = E_x::find_containment(1.0, a, d, w_ang, beta, k_x,
            gamma_x, k_y, gamma_y, mu_0, phi, psi, iterations);
        println!("Solutions.E_x[0..13] = {:?}", &solutions.E_x[0..13]);
        write!(file, "E_x modes\n");
        write!(file, "p,q,width/2 (a),depth/2 (d),cladding ref. index (n_0),\
        core ref. index (n_1),n_eff,Containment,k_x,gamma_x,k_y,gamma_y,\
        Propagation constant (beta)\n");
        write!(file, "{0},{1},{2},{3},{4},{5},{6},{12},{7},{8},{9},{10},{11}\n",
        solutions.E_x[0], solutions.E_x[1], solutions.E_x[2], solutions.E_x[3],
        solutions.E_x[4], solutions.E_x[5], solutions.E_x[6], solutions.E_x[7],
        solutions.E_x[8], solutions.E_x[9], solutions.E_x[10], solutions.E_x[11],
        solutions.E_x[12]);
        write!(file, "E_y modes\n");
        write!(file, "p,q,width/2 (a),depth/2 (d),cladding ref. index (n_0),\
        core ref. index (n_1),n_eff,Containment,k_x,gamma_x,k_y,gamma_y,\
        Propagation constant (beta)\n");
        write!(file, "{0},{1},{2},{3},{4},{5},{6},{12},{7},{8},{9},{10},{11}\n",
        solutions.E_y[0], solutions.E_y[1], solutions.E_y[2], solutions.E_y[3],
        solutions.E_y[4], solutions.E_y[5], solutions.E_y[6], solutions.E_y[7],
        solutions.E_y[8], solutions.E_y[9], solutions.E_y[10], solutions.E_y[11],
        solutions.E_y[12]);
        drop(file);

        //Plots
        if solutions.E_y[11] > 1.0 {
        //Plot for E_x
        let npoints: usize = 261;
        let mut sum_x: f64 = 0.0;
        let mut sum_y: f64 = 0.0;
        let mut x: f64 = 0.0;
        let mut y: f64 = 0.0;
        let mut Hx1: f64 = 0.0;
        let mut Hx2: f64 = 0.0;
        let mut Hx3: f64 = 0.0;
        let mut Hx4: f64 = 0.0;
        let mut Hx5: f64 = 0.0;
        let mut d2Hx1: f64 = 0.0;
        let mut d2Hx2: f64 = 0.0;
        let mut d2Hx3: f64 = 0.0;
        let mut d2Hx4: f64 = 0.0;
        let mut d2Hx5: f64 = 0.0;
        let mut Ey1: f64 = 0.0;
        let mut Ey2: f64 = 0.0;
        let mut Ey3: f64 = 0.0;
        let mut Ey4: f64 = 0.0;
        let mut Ey5: f64 = 0.0;
        let p: f64 = solutions.E_y[0];
        let q: f64 = solutions.E_y[1];
        let a: f64 = solutions.E_y[2];
        let d: f64 = solutions.E_y[3];
        let n_0: f64 = solutions.E_y[4];
        let n_1: f64 = solutions.E_y[5];
        let k_x: f64 = solutions.E_y[7];
        let gamma_x: f64 = solutions.E_y[8];
        let k_y: f64 = solutions.E_y[9];
        let gamma_y: f64 = solutions.E_y[10];
        let beta: f64 = solutions.E_y[11];
        let phi: f64 = E_y::phi(p, PI);
        let psi: f64 = E_y::psi(q, PI);
        let mut range_x: Vec<f64> = vec![0.0;npoints];
        let mut range_y: Vec<f64> = vec![0.0;npoints];
        let mut current_step_a: f64 = 0.0;
        let mut current_step_d: f64 = 0.0;
        let mut E_y_data: Vec<f64> = vec![0.0;(npoints * npoints)];
        let mut H_x_data: Vec<f64> = vec![0.0;(npoints * npoints)];
        let mut data_counter: usize = 0;

        for i in 0..npoints {
            current_step_a = -2.0 * a + a * 4.0 * (i as f64) / (npoints as f64);
            current_step_d = -2.0 * d + d * 4.0 * (i as f64) / (npoints as f64);
            range_x[i] = current_step_a;
            range_y[i] = current_step_d;
        }

        let name: &str = &format!("mode_plot_data/E_y_data_P{}Q{}.csv", p, q);
        let mut plot_file = match File::create(name) {
            Ok(contents) => contents,
            Err(why)     => panic!("Unable to create file: {}", why.description()),
        };
        write!(plot_file, "E_y data\n");
        write!(plot_file, "p, q, x (width), y (depth), z (Electric field intensity)\n");
        for x in range_x.clone().into_iter() { //Changes made in this whole block
            for y in range_y.clone().into_iter() {
                if (x >= -a && x <= a) && (y >= -d && y <= d){
                    H_x_data[data_counter] = E_y::Hx1(1.0, k_x, k_y, phi, psi, x, y);
                    Hx1 = H_x_data[data_counter];
                    d2Hx1 = E_y::d2Hx1(Hx1, k_y);
                    Ey1 = E_y::Ey1(Hx1, d2Hx1, w_ang, mu_0, beta, epsilon_0, n_1);//aquiaqui
                    E_y_data[data_counter] = Ey1;
                    write!(plot_file, "{}, {}, {}, {}, {}\n",
                    p, q, x, y, E_y_data[data_counter]);
                    data_counter += 1;
                } else if (x > a) && (y >= -d && y <= d) {
                    H_x_data[data_counter] = E_y::Hx2(1.0, k_x, k_y, gamma_x,
                    phi, psi, x, y, a);
                    Hx2 = H_x_data[data_counter];
                    d2Hx2 = E_y::d2Hx2(Hx2, k_y);
                    Ey2 = E_y::Ey2(Hx2, d2Hx2, w_ang, mu_0, beta, epsilon_0, n_0);
                    E_y_data[data_counter] = Ey2;
                    write!(plot_file, "{}, {}, {}, {}, {}\n",
                    p, q, x, y, E_y_data[data_counter]);
                    data_counter += 1;
                } else if (x >= -a && x <= a) && (y > d) {
                    H_x_data[data_counter] = E_y::Hx3(1.0, k_x, k_y, gamma_x,
                    phi, psi, x, y, d);
                    Hx3 = H_x_data[data_counter];
                    d2Hx3 = E_y::d2Hx3(Hx3, gamma_y);
                    Ey3 = E_y::Ey3(Hx3, d2Hx3, w_ang, mu_0, beta, epsilon_0, n_0);
                    E_y_data[data_counter] = Ey3;
                    write!(plot_file, "{}, {}, {}, {}, {}\n",
                    p, q, x, y, E_y_data[data_counter]);
                    data_counter += 1;
                } else if (x < -a) && (y >= -d && y <= d) {
                    H_x_data[data_counter] = E_y::Hx4(1.0, k_x, k_y, gamma_x,
                    phi, psi, x, y, a);
                    Hx4 = H_x_data[data_counter];
                    d2Hx4 = E_y::d2Hx4(Hx4, k_y);
                    Ey4 = E_y::Ey4(Hx4, d2Hx4, w_ang, mu_0, beta, epsilon_0, n_0);
                    E_y_data[data_counter] = Ey4;
                    write!(plot_file, "{}, {}, {}, {}, {}\n",
                    p, q, x, y, E_y_data[data_counter]);
                    data_counter += 1;
                } else if (x >= -a && x <= a) && (y < -d) {
                    H_x_data[data_counter] = E_y::Hx5(1.0, k_x, k_y, gamma_x,
                    phi, psi, x, y, d);
                    Hx5 = H_x_data[data_counter];
                    d2Hx5 = E_y::d2Hx5(Hx5, gamma_y);
                    Ey5 = E_y::Ey5(Hx5, d2Hx5, w_ang, mu_0, beta, epsilon_0, n_0);
                    E_y_data[data_counter] = Ey5;
                    write!(plot_file, "{}, {}, {}, {}, {}\n",
                    p, q, x, y, E_y_data[data_counter]);
                    data_counter += 1;
                } else {
                    E_y_data[data_counter] = 0.0;
                    write!(plot_file, "{}, {}, {}, {}, {}\n",
                    p, q, x, y, E_y_data[data_counter]);
                    data_counter += 1;
                }
            }
        }
        drop(plot_file);

        if booleans[1] == true {
        //Plot for E_x
        let mut figure = Figure::new();
        let mut figure2 = Figure::new();
        if booleans[3] == true {
            figure.axes3d()
            .set_title(&format!("Electric field E_y for P={}, Q ={}",p, q), &[])
            .surface(E_y_data.iter(), npoints, npoints,
            Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
            .set_z_label("E_y", &[])
            .set_x_label("Width", &[])
            .set_y_label("Depth",&[])
            .set_view(45.0, 45.0); //3D display
            if booleans[2] == false {
                figure.set_terminal("pngcairo", &format!("mode_graphs/E_y_field_P_{}_Q_{}.png", p, q));
            }
            figure.show();

            figure2.axes3d()
            .set_title(&format!("Magnetic field H_x for P={}, Q ={}",p, q), &[])
            .surface(H_x_data.iter(), npoints, npoints,
            Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
            .set_z_label("H_x", &[])
            .set_x_label("Width", &[])
            .set_y_label("Depth",&[])
            .set_view(45.0, 45.0); //3D display
            if booleans[2] == false {
                figure2.set_terminal("pngcairo", &format!("mode_graphs/H_x_field_P_{}_Q_{}.png", p, q));
            }
            figure2.show();

        } else {
            figure.axes3d()
            .set_title(&format!("Electric field E_y for P={}, Q ={}",p, q), &[])
            .surface(E_y_data.iter(), npoints, npoints,
            Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
            .set_z_label("E_y", &[])
            .set_x_label("Width", &[])
            .set_y_label("Depth",&[])
            .set_view_map();
            if booleans[2] == false {
                figure.set_terminal("pngcairo", &format!("mode_graphs/E_y_field_P_{}_Q_{}.png", p, q));
            }
            figure.show();

            figure2.axes3d()
            .set_title(&format!("Magnetic field H_x for P={}, Q ={}",p, q), &[])
            .surface(H_x_data.iter(), npoints, npoints,
            Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
            .set_z_label("H_x", &[])
            .set_x_label("Width", &[])
            .set_y_label("Depth",&[])
            .set_view_map();
            if booleans[2] == false {
                figure2.set_terminal("pngcairo", &format!("mode_graphs/H_x_field_P_{}_Q_{}.png", p, q));
            }
            figure2.show();
        }
        //End of if booleans[1] == true
        }
        //End if solutions.E_x[11] > 1.0
        }

        if solutions.E_x[11] > 1.0 {
        //Plot for E_x
        let npoints: usize = 261;
        let mut sum_x: f64 = 0.0;
        let mut sum_y: f64 = 0.0;
        let mut x: f64 = 0.0;
        let mut y: f64 = 0.0;
        let mut Hy1: f64 = 0.0;
        let mut Hy2: f64 = 0.0;
        let mut Hy3: f64 = 0.0;
        let mut Hy4: f64 = 0.0;
        let mut Hy5: f64 = 0.0;
        let mut d2Hy1: f64 = 0.0;
        let mut d2Hy2: f64 = 0.0;
        let mut d2Hy3: f64 = 0.0;
        let mut d2Hy4: f64 = 0.0;
        let mut d2Hy5: f64 = 0.0;
        let mut Ex1: f64 = 0.0;
        let mut Ex2: f64 = 0.0;
        let mut Ex3: f64 = 0.0;
        let mut Ex4: f64 = 0.0;
        let mut Ex5: f64 = 0.0;
        let p: f64 = solutions.E_x[0];
        let q: f64 = solutions.E_x[1];
        let a: f64 = solutions.E_x[2];
        let d: f64 = solutions.E_x[3];
        let n_0: f64 = solutions.E_x[4];
        let n_1: f64 = solutions.E_x[5];
        let k_x: f64 = solutions.E_x[7];
        let gamma_x: f64 = solutions.E_x[8];
        let k_y: f64 = solutions.E_x[9];
        let gamma_y: f64 = solutions.E_x[10];
        let beta: f64 = solutions.E_x[11];
        let phi: f64 = E_x::phi(p, PI);
        let psi: f64 = E_x::psi(q, PI);
        let mut range_x: Vec<f64> = vec![0.0;npoints];
        let mut range_y: Vec<f64> = vec![0.0;npoints];
        let mut current_step_a: f64 = 0.0;
        let mut current_step_d: f64 = 0.0;
        let mut E_x_data: Vec<f64> = vec![0.0;(npoints * npoints)];
        let mut H_y_data: Vec<f64> = vec![0.0;(npoints * npoints)];
        let mut data_counter: usize = 0;

        for i in 0..npoints {
            current_step_a = -2.0 * a + a * 4.0 * (i as f64) / (npoints as f64);
            current_step_d = -2.0 * d + d * 4.0 * (i as f64) / (npoints as f64);
            range_x[i] = current_step_a;
            range_y[i] = current_step_d;
        }

        let name: &str = &format!("mode_plot_data/E_x_data_P{}Q{}.csv", p, q);
        let mut plot_file = match File::create(name) {
            Ok(contents) => contents,
            Err(why)     => panic!("Unable to create file: {}", why.description()),
        };
        write!(plot_file, "E_y data\n");
        write!(plot_file, "p, q, x (width), y (depth), z (Electric field intensity)\n");
        for x in range_x.clone().into_iter() { //Changes made in this whole block
            for y in range_y.clone().into_iter() {
                if (x >= -a && x <= a) && (y >= -d && y <= d){
                    H_y_data[data_counter] = E_x::Hy1(1.0, k_x, k_y, phi, psi, x, y);
                    Hy1 = E_x_data[data_counter];
                    d2Hy1 = E_x::d2Hy1(1.0 , k_x, k_y, phi, psi, x, y);
                    Ex1 = E_x::Ex1(Hy1, d2Hy1, w_ang, mu_0, beta, epsilon_0, n_1);
                    E_x_data[data_counter] = Ex1;
                    write!(plot_file, "{}, {}, {}, {}, {}\n",
                    p, q, x, y, E_x_data[data_counter]);
                    data_counter += 1;
                } else if (x > a) && (y >= -d && y <= d) {
                    H_y_data[data_counter] = E_x::Hy2(1.0, k_x, k_y, gamma_x,
                    phi, psi, x, y, a);
                    Hy2 = E_x_data[data_counter];
                    d2Hy2 = E_x::d2Hy2(1.0 , k_x, k_y, gamma_x, phi, psi, x, y, a);
                    Ex2 = E_x::Ex2(Hy2, d2Hy2, w_ang, mu_0, beta, epsilon_0, n_0);
                    E_x_data[data_counter] = Ex2;
                    write!(plot_file, "{}, {}, {}, {}, {}\n",
                    p, q, x, y, E_x_data[data_counter]);
                    data_counter += 1;
                } else if (x >= -a && x <= a) && (y > d) {
                    H_y_data[data_counter] = E_x::Hy3(1.0, k_x, k_y, gamma_x,
                    phi, psi, x, y, d);
                    Hy3 = E_x_data[data_counter];
                    d2Hy3 = E_x::d2Hy3(1.0, k_x, k_y, gamma_y, phi, psi, x, y, d);
                    Ex3 = E_x::Ex3(Hy3, d2Hy3, w_ang, mu_0, beta, epsilon_0, n_0);
                    E_x_data[data_counter] = Ex3;
                    write!(plot_file, "{}, {}, {}, {}, {}\n",
                    p, q, x, y, E_x_data[data_counter]);
                    data_counter += 1;
                } else if (x < -a) && (y >= -d && y <= d) {
                    H_y_data[data_counter] = E_x::Hy4(1.0, k_x, k_y, gamma_x,
                    phi, psi, x, y, a);
                    Hy4 = E_x_data[data_counter];
                    d2Hy4 = E_x::d2Hy4(1.0, k_x, k_y, gamma_x, phi, psi, x, y, a);
                    Ex4 = E_x::Ex4(Hy4, d2Hy4, w_ang, mu_0, beta, epsilon_0, n_0);
                    E_x_data[data_counter] = Ex4;
                    write!(plot_file, "{}, {}, {}, {}, {}\n",
                    p, q, x, y, E_x_data[data_counter]);
                    data_counter += 1;
                } else if (x >= -a && x <= a) && (y < -d) {
                    H_y_data[data_counter] = E_x::Hy5(1.0, k_x, k_y, gamma_x,
                    phi, psi, x, y, d);
                    Hy5 = E_x_data[data_counter];
                    d2Hy5 = E_x::d2Hy5(1.0, k_x, k_y, gamma_y, phi, psi, x, y, d);
                    Ex5 = E_x::Ex5(Hy5, d2Hy5, w_ang, mu_0, beta, epsilon_0, n_0);
                    E_x_data[data_counter] = Ex5;
                    write!(plot_file, "{}, {}, {}, {}, {}\n",
                    p, q, x, y, E_x_data[data_counter]);
                    data_counter += 1;
                } else {
                    E_x_data[data_counter] = 0.0;
                    write!(plot_file, "{}, {}, {}, {}, {}\n",
                    p, q, x, y, E_x_data[data_counter]);
                    data_counter += 1;
                }
            }
        }
        drop(plot_file);

        if booleans[1] == true {
        //Plot for E_x
        let mut figure = Figure::new();
        let mut figure2 = Figure::new();
        if booleans[3] == true {
            figure.axes3d()
            .set_title(&format!("Electric field E_x for P={}, Q ={}",p, q), &[])
            .surface(E_x_data.iter(), npoints, npoints,
            Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
            .set_z_label("E_x", &[])
            .set_x_label("Width", &[])
            .set_y_label("Depth",&[])
            .set_view(45.0, 45.0); //3D display
            if booleans[2] == false {
                figure.set_terminal("pngcairo", &format!("mode_graphs/E_x_field_P_{}_Q_{}.png", p, q));
            }
            figure.show();

            figure2.axes3d()
            .set_title(&format!("Magnetic field H_y for P={}, Q ={}",p, q), &[])
            .surface(H_y_data.iter(), npoints, npoints,
            Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
            .set_z_label("H_y", &[])
            .set_x_label("Width", &[])
            .set_y_label("Depth",&[])
            .set_view(45.0, 45.0); //3D display
            if booleans[2] == false {
                figure2.set_terminal("pngcairo", &format!("mode_graphs/H_y_field_P_{}_Q_{}.png", p, q));
            }
            figure2.show();

        } else {
            figure.axes3d()
            .set_title(&format!("Electric field E_x for P={}, Q ={}",p, q), &[])
            .surface(E_x_data.iter(), npoints, npoints,
            Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
            .set_z_label("E_x", &[])
            .set_x_label("Width", &[])
            .set_y_label("Depth",&[])
            .set_view_map();
            if booleans[2] == false {
                figure.set_terminal("pngcairo", &format!("mode_graphs/E_x_field_P_{}_Q_{}.png", p, q));
            }
            figure.show();

            figure2.axes3d()
            .set_title(&format!("Magnetic field H_y for P={}, Q ={}",p, q), &[])
            .surface(H_y_data.iter(), npoints, npoints,
            Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
            .set_z_label("H_y", &[])
            .set_x_label("Width", &[])
            .set_y_label("Depth",&[])
            .set_view_map();
            if booleans[2] == false {
                figure2.set_terminal("pngcairo", &format!("mode_graphs/H_y_field_P_{}_Q_{}.png", p, q));
            }
            figure2.show();
        }
        //End of if booleans[1] == true
        }
        //End if solutions.E_x[11] > 1.0
        }

    //End of first part of if block, else block next
    } else {

    let STEPS: i64 = params.STEPS;
    let a: f64 = params.a;
    let d: f64 = params.d;
    let w_ang: f64 = params.w_ang;
    let k: f64 = params.k;
    let P: usize = params.P;
    let Q: usize = params.Q;
    let iterations: i64 = (STEPS as f64).sqrt() as i64;
    let n_0: f64 = params.n_0;
    let n_1: f64 = params.n_1;
    let STEPSIZE: f64 = k * n_1 / (STEPS as f64);
    let mut index: usize = 0;


    //EX calculation//
    for p in 0..P {
        newlow = 1.0e20;
        for j in 0..STEPS {
            let k_x: f64 = (j as f64) * STEPSIZE;
            let phi: f64 = E_x::phi((p+1) as f64, PI);
            let gamma_x: f64 = E_x::gamma_x(n_0, n_1, k, k_x);
            comparator = E_x::find_k_x(n_0, n_1, (p+1) as f64, a, PI, gamma_x, k_x);
            if comparator < 1.0e-3 && comparator < newlow && k_x > 10.0  {
                solutions.sols_x[0 + 4 * p] = (p+1) as f64;
                solutions.sols_x[1 + 4 * p] = k_x;
                solutions.sols_x[2 + 4 * p] = gamma_x;
                solutions.sols_x[3 + 4 * p] = comparator;
                newlow = comparator;
            }
        }
    }
    println!("E_x solutions.sols_x = {:?}", solutions.sols_x);

    for q in 0..Q {
        newlow = 1.0e20;
        for j in 0..STEPS {
            let k_y: f64 = (j as f64) * STEPSIZE;
            let psi: f64 = E_x::psi((q+1) as f64, PI);
            let gamma_y: f64 = E_x::gamma_y(n_0, n_1, k, k_y);
            comparator = E_x::find_k_y(n_0, n_1, (q+1) as f64, d, PI, gamma_y, k_y);
            if comparator < 1.0e-3 && comparator < newlow && k_y > 10.0  {
                solutions.sols_y[0 + 4 * q] = (q+1) as f64;
                solutions.sols_y[1 + 4 * q] = k_y;
                solutions.sols_y[2 + 4 * q] = gamma_y;
                solutions.sols_y[3 + 4 * q] = comparator;
                newlow = comparator;
            }
        }
    }
    println!("E_x solutions.sols_y = {:?}", solutions.sols_y);

    index = 0;
    for p in 0..P {
        for q in 0..Q {
            let k_x: f64 = solutions.sols_x[1 + 4 * p];
            let gamma_x: f64 = solutions.sols_x[2 + 4 * p];
            let k_y: f64 = solutions.sols_y[1 + 4 * q];
            let gamma_y: f64 = solutions.sols_y[2 + 4 * q];
            let beta: f64 = E_x::find_beta(n_0, n_1, k, k_x, k_y);
            let n_eff: f64 = E_x::find_n_eff(beta, k);
            let phi: f64 = E_x::phi((p+1) as f64, PI);
            let psi: f64 = E_x::phi((q+1) as f64, PI);
            let comparator_2: f64 = ((n_eff / n_0) - 1.0) * 100.0;
            println!("Finding E_x solutions P={}, Q={}", p+1,q+1);
            if (n_eff > n_0 && n_eff < n_1 /*&& comparator_2 > 0.1*/) && (k_x > 10.0 && k_y > 10.0) {
                solutions.E_x[0 + index * 13] = (p+1) as f64;
                solutions.E_x[1 + index * 13] = (q+1) as f64;
                solutions.E_x[2 + index * 13] = a;
                solutions.E_x[3 + index * 13] = d;
                solutions.E_x[4 + index * 13] = n_0;
                solutions.E_x[5 + index * 13] = n_1;
                solutions.E_x[6 + index * 13] = n_eff;
                solutions.E_x[7 + index * 13] = k_x;
                solutions.E_x[8 + index * 13] = gamma_x;
                solutions.E_x[9 + index * 13] = k_y;
                solutions.E_x[10 + index * 13] = gamma_y;
                solutions.E_x[11 + index * 13] = beta;
                solutions.E_x[12 + index * 13] = E_x::find_containment(1.0, a,
                    d, w_ang, beta, k_x, gamma_x, k_y, gamma_y, mu_0, phi, psi,
                    iterations);
                index += 1;
                }
            }
        }

    //Resetting vectors
    solutions.sols_x = vec![0.0; params.P * 4];
    solutions.sols_y = vec![0.0; params.Q * 4];

    //EY calculation//
    for p in 0..P {
        newlow = 1.0e20;
        for j in 0..STEPS {
            let k_x: f64 = (j as f64) * STEPSIZE;
            let phi: f64 = E_y::phi((p+1) as f64, PI);
            let gamma_x: f64 = E_y::gamma_x(n_0, n_1, k, k_x);
            comparator = E_y::find_k_x(n_0, n_1, (p+1) as f64, a, PI, gamma_x, k_x);
            if comparator < 1.0e-3 && comparator < newlow && k_x > 10.0  {
                solutions.sols_x[0 + 4 * p] = (p+1) as f64;
                solutions.sols_x[1 + 4 * p] = k_x;
                solutions.sols_x[2 + 4 * p] = gamma_x;
                solutions.sols_x[3 + 4 * p] = comparator;
                newlow = comparator;
            }
        }
    }

    println!("E_y solutions.sols_x = {:?}", solutions.sols_x);

    for q in 0..Q {
        newlow = 1.0e20;
        for j in 0..STEPS {
            let k_y: f64 = (j as f64) * STEPSIZE;
            let psi: f64 = E_y::psi((q+1) as f64, PI);
            let gamma_y: f64 = E_y::gamma_y(n_0, n_1, k, k_y);
            comparator = E_y::find_k_y(n_0, n_1, (q+1) as f64, d, PI, gamma_y, k_y);
            if comparator < 1.0e-3 && comparator < newlow && k_y > 10.0  {
                solutions.sols_y[0 + 4 * q] = (q+1) as f64;
                solutions.sols_y[1 + 4 * q] = k_y;
                solutions.sols_y[2 + 4 * q] = gamma_y;
                solutions.sols_y[3 + 4 * q] = comparator;
                newlow = comparator;
            }
        }
    }

    println!("E_y solutions.sols_y = {:?}", solutions.sols_y);

    index = 0;
    for p in 0..P {
        for q in 0..Q {
            let k_x: f64 = solutions.sols_x[1 + 4 * p];
            let gamma_x: f64 = solutions.sols_x[2 + 4 * p];
            let k_y: f64 = solutions.sols_y[1 + 4 * q];
            let gamma_y: f64 = solutions.sols_y[2 + 4 * q];
            let beta: f64 = E_y::find_beta(n_0, n_1, k, k_x, k_y);
            let n_eff: f64 = E_y::find_n_eff(beta, k);
            let phi: f64 = E_y::phi((p+1) as f64, PI);
            let psi: f64 = E_y::phi((q+1) as f64, PI);
            let comparator_2: f64 = ((n_eff / n_0) - 1.0) * 100.0;
            println!("Finding E_y solutions P={}, Q={}", p+1,q+1);
            if (n_eff > n_0 && n_eff < n_1 /*&& comparator_2 > 0.1*/) && (k_x > 10.0 && k_y > 10.0) {
                println!("comparator_2 = {}, n_eff = {}, n_0 = {}", comparator_2, n_eff, n_0);
                solutions.E_y[0 + index * 13] = (p+1) as f64;
                solutions.E_y[1 + index * 13] = (q+1) as f64;
                solutions.E_y[2 + index * 13] = a;
                solutions.E_y[3 + index * 13] = d;
                solutions.E_y[4 + index * 13] = n_0;
                solutions.E_y[5 + index * 13] = n_1;
                solutions.E_y[6 + index * 13] = n_eff;
                solutions.E_y[7 + index * 13] = k_x;
                solutions.E_y[8 + index * 13] = gamma_x;
                solutions.E_y[9 + index * 13] = k_y;
                solutions.E_y[10 + index * 13] = gamma_y;
                solutions.E_y[11 + index * 13] = beta;
                solutions.E_y[12 + index * 13] = E_y::find_containment(1.0, a,
                    d, w_ang, beta, k_x, gamma_x, k_y, gamma_y, mu_0, phi, psi,
                    iterations);
                index += 1;
                }
            }
        }

        write!(file, "E_x modes\n");
        write!(file, "p,q,width/2 (a),depth/2 (d),cladding ref. index (n_0),\
        core ref. index (n_1),n_eff,Containment,k_x,gamma_x,k_y,gamma_y,\
        Propagation constant (beta)\n");
        let vec_len: usize = solutions.E_x.len() / 13;
        for i in 0..vec_len {
            if solutions.E_x[0 + i * 13] > 0.0 {
                write!(file, "{0},{1},{2},{3},{4},{5},{6},{12},{7},{8},{9},{10},{11}\n",
                solutions.E_x[0 + i * 13], solutions.E_x[1+ i * 13],
                solutions.E_x[2 + i * 13], solutions.E_x[3 + i * 13],
                solutions.E_x[4 + i * 13], solutions.E_x[5 + i * 13],
                solutions.E_x[6 + i * 13], solutions.E_x[7 + i * 13],
                solutions.E_x[8 + i * 13], solutions.E_x[9 + i * 13],
                solutions.E_x[10 + i * 13], solutions.E_x[11 + i * 13],
                solutions.E_x[12 + i * 13]);
            }

        }

        write!(file, "E_y modes\n");
        write!(file, "p,q,width/2 (a),depth/2 (d),cladding ref. index (n_0),\
        core ref. index (n_1),n_eff,Containment,k_x,gamma_x,k_y,gamma_y,\
        Propagation constant (beta)\n");
        let vec_len: usize = solutions.E_x.len() / 13;
        for i in 0..vec_len {
            if solutions.E_y[0 + i * 13] > 0.0 {
                write!(file, "{0},{1},{2},{3},{4},{5},{6},{12},{7},{8},{9},{10},{11}\n",
                solutions.E_y[0 + i * 13], solutions.E_y[1 + i * 13],
                solutions.E_y[2 + i * 13], solutions.E_y[3 + i * 13],
                solutions.E_y[4 + i * 13], solutions.E_y[5 + i * 13],
                solutions.E_y[6 + i * 13], solutions.E_y[7 + i * 13],
                solutions.E_y[8 + i * 13], solutions.E_y[9 + i * 13],
                solutions.E_y[10 + i * 13], solutions.E_y[11 + i * 13],
                solutions.E_y[12 + i * 13]);
            }

        }
        drop(file);

        if booleans[1] == true {
        //Plots
        index = 0;
        for _ in 0..P{
            for _ in 0..Q{
                if solutions.E_y[11 + index * 13] > 1.0 {

                //Plot for E_y
                let npoints: usize = 261;
                let mut sum_x: f64 = 0.0;
                let mut sum_y: f64 = 0.0;
                let mut x: f64 = 0.0;
                let mut y: f64 = 0.0;
                let p: f64 = solutions.E_y[0 + index * 13];
                let q: f64 = solutions.E_y[1 + index * 13];
                let a: f64 = solutions.E_y[2 + index * 13];
                let d: f64 = solutions.E_y[3 + index * 13];
                let k_x: f64 = solutions.E_y[7 + index * 13];
                let gamma_x: f64 = solutions.E_y[8 + index * 13];
                let k_y: f64 = solutions.E_y[9 + index * 13];
                let gamma_y: f64 = solutions.E_y[10 + index * 13];
                let beta: f64 = solutions.E_y[11 + index * 13];
                let phi: f64 = E_y::phi(p, PI);
                let psi: f64 = E_y::psi(q, PI);
                let mut range_x: Vec<f64> = vec![0.0;npoints];
                let mut range_y: Vec<f64> = vec![0.0;npoints];
                let mut current_step_a: f64 = 0.0;
                let mut current_step_d: f64 = 0.0;
                let mut E_y_data: Vec<f64> = vec![0.0;(npoints * npoints)];
                let mut data_counter: usize = 0;

                for i in 0..npoints {
                    current_step_a = -2.0 * a + a * 4.0 * (i as f64) / (npoints as f64);
                    current_step_d = -2.0 * d + d * 4.0 * (i as f64) / (npoints as f64);
                    range_x[i] = current_step_a;
                    range_y[i] = current_step_d;
                }

                let name: &str = &format!("mode_plot_data/E_y_data_P{}Q{}.csv", p, q);
                let mut plot_file = match File::create(name) {
                    Ok(contents) => contents,
                    Err(why)     => panic!("Unable to create file: {}", why.description()),
                };
                write!(plot_file, "E_y data\n");
                write!(plot_file, "p, q, x (width), y (depth), z (Electric field intensity)\n");
                for x in range_x.clone().into_iter() {
                    for y in range_y.clone().into_iter() {
                        if (x >= -a && x <= a) && (y >= -d && y <= d){
                            E_y_data[data_counter] = E_y::Hx1(1.0, k_x, k_y, phi, psi, x, y);
                            write!(plot_file, "{}, {}, {}, {}, {}\n",
                            p, q, x, y, E_y_data[data_counter]);
                            data_counter += 1;
                        } else if (x > a) && (y >= -d && y <= d) {
                            E_y_data[data_counter] = E_y::Hx2(1.0, k_x, k_y, gamma_x,
                            phi, psi, x, y, a);
                            write!(plot_file, "{}, {}, {}, {}, {}\n",
                            p, q, x, y, E_y_data[data_counter]);
                            data_counter += 1;
                        } else if (x >= -a && x <= a) && (y > d) {
                            E_y_data[data_counter] = E_y::Hx3(1.0, k_x, k_y, gamma_x,
                            phi, psi, x, y, d);
                            write!(plot_file, "{}, {}, {}, {}, {}\n",
                            p, q, x, y, E_y_data[data_counter]);
                            data_counter += 1;
                        } else if (x < -a) && (y >= -d && y <= d) {
                            E_y_data[data_counter] = E_y::Hx4(1.0, k_x, k_y, gamma_x,
                            phi, psi, x, y, a);
                            write!(plot_file, "{}, {}, {}, {}, {}\n",
                            p, q, x, y, E_y_data[data_counter]);
                            data_counter += 1;
                        } else if (x >= -a && x <= a) && (y < -d) {
                            E_y_data[data_counter] = E_y::Hx5(1.0, k_x, k_y, gamma_x,
                            phi, psi, x, y, d);
                            write!(plot_file, "{}, {}, {}, {}, {}\n",
                            p, q, x, y, E_y_data[data_counter]);
                            data_counter += 1;
                        } else {
                            E_y_data[data_counter] = 0.0;
                            write!(plot_file, "{}, {}, {}, {}, {}\n",
                            p, q, x, y, E_y_data[data_counter]);
                            data_counter += 1;
                        }
                    }
                }
                drop(plot_file);

                if booleans[1] == true {
                //Plot for E_y
                let mut figure = Figure::new();
                let mut figure2 = Figure::new();
                if booleans[3] == true {
                    figure.axes3d()
                    .set_title(&format!("Electric field E_y for P={}, Q ={}",p, q), &[])
                    .surface(E_y_data.iter(), npoints, npoints,
                    Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
                    .set_z_label("E_y", &[])
                    .set_x_label("Width", &[])
                    .set_y_label("Depth",&[])
                    .set_view(45.0, 45.0); //3D display
                    if booleans[2] == false {
                        figure.set_terminal("pngcairo", &format!("mode_graphs/E_y_field_P_{}_Q_{}.png", p, q));
                    }
                    figure.show();
                } else {
                    figure.axes3d()
                    .set_title(&format!("Electric field E_y for P={}, Q ={}",p, q), &[])
                    .surface(E_y_data.iter(), npoints, npoints,
                    Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
                    .set_z_label("E_y", &[])
                    .set_x_label("Width", &[])
                    .set_y_label("Depth",&[])
                    .set_view_map();
                    if booleans[2] == false {
                        figure.set_terminal("pngcairo", &format!("mode_graphs/E_y_field_P_{}_Q_{}.png", p, q));
                    }
                    figure.show();
                }
                //end of if booleans[1] == true
                }
                //end of E_y plots
                }

            if solutions.E_x[11 + index * 13] > 1.0 {
                //Plot for E_x
                let npoints: usize = 261;
                let mut sum_x: f64 = 0.0;
                let mut sum_y: f64 = 0.0;
                let mut x: f64 = 0.0;
                let mut y: f64 = 0.0;
                let mut Hy1: f64 = 0.0;
                let mut Hy2: f64 = 0.0;
                let mut Hy3: f64 = 0.0;
                let mut Hy4: f64 = 0.0;
                let mut Hy5: f64 = 0.0;
                let mut d2Hy1: f64 = 0.0;
                let mut d2Hy2: f64 = 0.0;
                let mut d2Hy3: f64 = 0.0;
                let mut d2Hy4: f64 = 0.0;
                let mut d2Hy5: f64 = 0.0;
                let mut Ex1: f64 = 0.0;
                let mut Ex2: f64 = 0.0;
                let mut Ex3: f64 = 0.0;
                let mut Ex4: f64 = 0.0;
                let mut Ex5: f64 = 0.0;
                let p: f64 = solutions.E_x[0 + index * 13];
                let q: f64 = solutions.E_x[1 + index * 13];
                let a: f64 = solutions.E_x[2 + index * 13];
                let d: f64 = solutions.E_x[3 + index * 13];
                let n_0: f64 = solutions.E_x[4 + index * 13];
                let n_1: f64 = solutions.E_x[5 + index * 13];
                let k_x: f64 = solutions.E_x[7 + index * 13];
                let gamma_x: f64 = solutions.E_x[8 + index * 13];
                let k_y: f64 = solutions.E_x[9 + index * 13];
                let gamma_y: f64 = solutions.E_x[10 + index * 13];
                let beta: f64 = solutions.E_x[11 + index * 13];
                let phi: f64 = E_x::phi(p, PI);
                let psi: f64 = E_x::psi(q, PI);
                let mut range_x: Vec<f64> = vec![0.0;npoints];
                let mut range_y: Vec<f64> = vec![0.0;npoints];
                let mut current_step_a: f64 = 0.0;
                let mut current_step_d: f64 = 0.0;
                let mut E_x_data: Vec<f64> = vec![0.0;(npoints * npoints)];
                let mut H_y_data: Vec<f64> = vec![0.0;(npoints * npoints)];
                let mut data_counter: usize = 0;

                for i in 0..npoints {
                    current_step_a = -2.0 * a + a * 4.0 * (i as f64) / (npoints as f64);
                    current_step_d = -2.0 * d + d * 4.0 * (i as f64) / (npoints as f64);
                    range_x[i] = current_step_a;
                    range_y[i] = current_step_d;
                }

                let name: &str = &format!("mode_plot_data/E_x_data_P{}Q{}.csv", p, q);
                let mut plot_file = match File::create(name) {
                    Ok(contents) => contents,
                    Err(why)     => panic!("Unable to create file: {}", why.description()),
                };
                write!(plot_file, "E_y data\n");
                write!(plot_file, "p, q, x (width), y (depth), z (Electric field intensity)\n");
                for x in range_x.clone().into_iter() { //Changes made in this whole block
                    for y in range_y.clone().into_iter() {
                        if (x >= -a && x <= a) && (y >= -d && y <= d){
                            H_y_data[data_counter] = E_x::Hy1(1.0, k_x, k_y, phi, psi, x, y);
                            Hy1 = E_x_data[data_counter];
                            d2Hy1 = E_x::d2Hy1(1.0 , k_x, k_y, phi, psi, x, y);
                            Ex1 = E_x::Ex1(Hy1, d2Hy1, w_ang, mu_0, beta, epsilon_0, n_1);
                            E_x_data[data_counter] = Ex1;
                            write!(plot_file, "{}, {}, {}, {}, {}\n",
                            p, q, x, y, E_x_data[data_counter]);
                            data_counter += 1;
                        } else if (x > a) && (y >= -d && y <= d) {
                            H_y_data[data_counter] = E_x::Hy2(1.0, k_x, k_y, gamma_x,
                            phi, psi, x, y, a);
                            Hy2 = E_x_data[data_counter];
                            d2Hy2 = E_x::d2Hy2(1.0 , k_x, k_y, gamma_x, phi, psi, x, y, a);
                            Ex2 = E_x::Ex2(Hy2, d2Hy2, w_ang, mu_0, beta, epsilon_0, n_0);
                            E_x_data[data_counter] = Ex2;
                            write!(plot_file, "{}, {}, {}, {}, {}\n",
                            p, q, x, y, E_x_data[data_counter]);
                            data_counter += 1;
                        } else if (x >= -a && x <= a) && (y > d) {
                            H_y_data[data_counter] = E_x::Hy3(1.0, k_x, k_y, gamma_x,
                            phi, psi, x, y, d);
                            Hy3 = E_x_data[data_counter];
                            d2Hy3 = E_x::d2Hy3(1.0, k_x, k_y, gamma_y, phi, psi, x, y, d);
                            Ex3 = E_x::Ex3(Hy3, d2Hy3, w_ang, mu_0, beta, epsilon_0, n_0);
                            E_x_data[data_counter] = Ex3;
                            write!(plot_file, "{}, {}, {}, {}, {}\n",
                            p, q, x, y, E_x_data[data_counter]);
                            data_counter += 1;
                        } else if (x < -a) && (y >= -d && y <= d) {
                            H_y_data[data_counter] = E_x::Hy4(1.0, k_x, k_y, gamma_x,
                            phi, psi, x, y, a);
                            Hy4 = E_x_data[data_counter];
                            d2Hy4 = E_x::d2Hy4(1.0, k_x, k_y, gamma_x, phi, psi, x, y, a);
                            Ex4 = E_x::Ex4(Hy4, d2Hy4, w_ang, mu_0, beta, epsilon_0, n_0);
                            E_x_data[data_counter] = Ex4;
                            write!(plot_file, "{}, {}, {}, {}, {}\n",
                            p, q, x, y, E_x_data[data_counter]);
                            data_counter += 1;
                        } else if (x >= -a && x <= a) && (y < -d) {
                            H_y_data[data_counter] = E_x::Hy5(1.0, k_x, k_y, gamma_x,
                            phi, psi, x, y, d);
                            Hy5 = E_x_data[data_counter];
                            d2Hy5 = E_x::d2Hy5(1.0, k_x, k_y, gamma_y, phi, psi, x, y, d);
                            Ex5 = E_x::Ex5(Hy5, d2Hy5, w_ang, mu_0, beta, epsilon_0, n_0);
                            E_x_data[data_counter] = Ex5;
                            write!(plot_file, "{}, {}, {}, {}, {}\n",
                            p, q, x, y, E_x_data[data_counter]);
                            data_counter += 1;
                        } else {
                            E_x_data[data_counter] = 0.0;
                            write!(plot_file, "{}, {}, {}, {}, {}\n",
                            p, q, x, y, E_x_data[data_counter]);
                            data_counter += 1;
                        }
                    }
                }
                drop(plot_file);

                if booleans[1] == true {
                //Plot for E_x
                let mut figure = Figure::new();
                let mut figure2 = Figure::new();
                if booleans[3] == true {
                    figure.axes3d()
                    .set_title(&format!("Electric field E_x for P={}, Q ={}",p, q), &[])
                    .surface(E_x_data.iter(), npoints, npoints,
                    Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
                    .set_z_label("E_x", &[])
                    .set_x_label("Width", &[])
                    .set_y_label("Depth",&[])
                    .set_view(45.0, 45.0); //3D display
                    if booleans[2] == false {
                        figure.set_terminal("pngcairo", &format!("mode_graphs/E_x_field_P_{}_Q_{}.png", p, q));
                    }
                    figure.show();

                    figure2.axes3d()
                    .set_title(&format!("Magnetic field H_y for P={}, Q ={}",p, q), &[])
                    .surface(H_y_data.iter(), npoints, npoints,
                    Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
                    .set_z_label("H_y", &[])
                    .set_x_label("Width", &[])
                    .set_y_label("Depth",&[])
                    .set_view(45.0, 45.0); //3D display
                    if booleans[2] == false {
                        figure2.set_terminal("pngcairo", &format!("mode_graphs/H_y_field_P_{}_Q_{}.png", p, q));
                    }
                    figure2.show();

                } else {
                    figure.axes3d()
                    .set_title(&format!("Electric field E_x for P={}, Q ={}",p, q), &[])
                    .surface(E_x_data.iter(), npoints, npoints,
                    Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
                    .set_z_label("E_x", &[])
                    .set_x_label("Width", &[])
                    .set_y_label("Depth",&[])
                    .set_view_map();
                    if booleans[2] == false {
                        figure.set_terminal("pngcairo", &format!("mode_graphs/E_x_field_P_{}_Q_{}.png", p, q));
                    }
                    figure.show();

                    figure2.axes3d()
                    .set_title(&format!("Magnetic field H_y for P={}, Q ={}",p, q), &[])
                    .surface(H_y_data.iter(), npoints, npoints,
                    Some((-2.0 * a * 1e6, -2.0 * d * 1e6, 2.0 * a * 1e6, 2.0 * d * 1e6)), &[])
                    .set_z_label("H_y", &[])
                    .set_x_label("Width", &[])
                    .set_y_label("Depth",&[])
                    .set_view_map();
                    if booleans[2] == false {
                        figure2.set_terminal("pngcairo", &format!("mode_graphs/H_y_field_P_{}_Q_{}.png", p, q));
                    }
                    figure2.show();
                    }
                //End of if booleans[1] == true
                }
                //End of E_x plots
                }
            index +=1;
            //End of for q in 0..Q
            }
        //End of for p in 0..P
        }
        //End of if booleans[1] == true
        }
    //End of if (else) block
    }



//End of calculation block
}

//End of main
}
