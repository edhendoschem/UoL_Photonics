//Power containment page 34, Marcatilis starts in page 39
extern crate rand;
extern crate gnuplot;

use rand::distributions::{IndependentSample, Range};
use gnuplot::*;
use std::{f64, i64};
use std::f64::consts;
use std::error::Error;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

///////////////////////////////////CONSTANTS////////////////////////////////////
const STEPS: usize = 1000000;
//Calculation
const PI: f64 = consts::PI;
const mu_0: f64 = (4.0 * PI) * 1e-7;                                            // Permeability of free space in H/m
const epsilon_0: f64 = 8.854187817 * 1e-12;                                     // Permittivity of free space in F/m
const h: f64 = 6.626070e-34;                                                    // Planck's constant, J.s
const c: f64 =  299792458.0;                                                    // Speed of light, m/s

//////////////////////////////////STRUCTURES////////////////////////////////////
struct starting_params {
    n_0: f64,                                                                   // Refractive index of the cladding
    n_1: f64,                                                                   // Refractive index of the core
    wavelength: f64,                                                            // Wavelength in m
    frequency: f64,                                                             // Frequency in Hz
    w_ang: f64,                                                                // Angular frequency in rad/s
    k: f64,                                                                     // Wavenumber in Vacuum
    d: f64,                                                                     // depth
    a: f64,                                                                     // width
    p: f64,                                                                     // Modes in the x direction
    q: f64,                                                                     // Modes in the y direction
}

struct vectors {
    k_x: [f64;4],
    k_y: [f64;4],
    gamma_x: [f64; 4],
    gamma_y: [f64; 4],
    beta: [[f64;4];4],
    sols_Ex: [[f64; 13];16], //0=p, 1=q, 2=n_0, 3=n_1, 4=kx, 5=ky, 6=gamma_x, 7=gamma_y, 8=beta, 9=n_eff, 10=cont_x,
    sols_Ey: [[f64; 13];16], //0=p, 1=q, 2=n_0, 3=n_1, 4=kx, 5=ky, 6=gamma_x, 7=gamma_y, 8=beta, 9=n_eff, 10=cont_y,
}
//////////////////////////////////FUNCTIONS/////////////////////////////////////
fn find_freq(wavelength: f64) -> f64 {
    c / wavelength
}

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
//////////////////////////////////Mods//////////////////////////////////////////
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
    /////////////////////////////////PARAMETER SET UP///////////////////////////
    let mut params: starting_params = starting_params {
        n_0: 1.45,                                                              // Refractive index of the cladding
        n_1: 1.50,                                                             // Refractive index of the core
        wavelength: 1536.0 * 1e-9,                                              // Wavelength in m
        frequency: 0.0,                                                         // Frequency in Hz
        w_ang: 0.0,                                                             // Angular frequency in rad/s
        k: 0.0,                                                                 // Wavenumber in Vacuum
        d: 5.0e-6,                                                                     // depth
        a: 5.0e-6,                                                                     // width
        p: 4.0,                                                                     // Modes in the x direction
        q: 4.0,                                                                     // Modes in the y direction
    };

    let mut vectors: vectors = vectors {
        k_x: [0.0; 4],
        k_y: [0.0; 4],
        gamma_x: [0.0; 4],
        gamma_y: [0.0; 4],
        beta:[[0.0; 4]; 4],
        sols_Ex: [[0.0; 13]; 16], //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
        sols_Ey: [[0.0; 13]; 16], //0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
    };

    params.frequency = find_freq(params.wavelength);
    params.w_ang = find_ang_freq(params.frequency);
    params.k = params.w_ang / c;

    let mut a: f64 = params.a;
    let mut d: f64 = params.d;
    let mut k: f64 = params.k;
    let mut n_0: f64 = params.n_0;
    let mut n_1: f64 = params.n_1;
    let mut stepsize: f64 = k * n_1 / (STEPS as f64);
    let mut P: usize = params.p as usize;
    let mut Q: usize = params.q as usize;
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
        for i in 0..STEPS {
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
        for i in 0..STEPS {
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
//start of block
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
                            sum_x += E_x::Hy1(1.0, k_x, k_y, phi, psi, x, y)

                        }
                        sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                        sum_x = 0.0;
                    }
                }
                let P1: f64 = c1 * (2.0 * d) * (sum_y).powi(2)/((N as f64) * (cycles as f64));
                sum_x = 0.0;
                sum_y = 0.0;
                let range_x  = Range::new(a, 3.0 * a);
                let range_y  = Range::new(-1.0 * d, d);
                let mut rng = rand::thread_rng();
                //Region 2, -a<=x<=a, -d<=y<=d
                for z in 0..cycles {
                    for i in 0..N {
                        y = range_y.ind_sample(&mut rng);
                        for j in 0..N {
                            x = range_x.ind_sample(&mut rng);
                            sum_x += E_x::Hy2(1.0, k_x, k_y, gamma_x, phi, psi, x, y, a)

                        }
                        sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                        sum_x = 0.0;
                    }
                }

                let P2: f64 = c2 * (2.0 * d) * (sum_y).powi(2)/((N as f64) * (cycles as f64));
                let range_x  = Range::new(-1.0 * a, a);
                let range_y  = Range::new(d, 3.0 * d);
                let mut rng = rand::thread_rng();
                sum_x = 0.0;
                sum_y = 0.0;
                //Region 3, -a<=x<=a, -d<=y<=d
                for z in 0..cycles {
                    for i in 0..N {
                        y = range_y.ind_sample(&mut rng);
                        for j in 0..N {
                            x = range_x.ind_sample(&mut rng);
                            sum_x += E_x::Hy3(1.0, k_x, k_y, gamma_y, phi, psi, x, y, d)

                        }
                        sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                        sum_x = 0.0;
                    }
                }
                let P3: f64 = c3 * (2.0 * d) * (sum_y).powi(2)/((N as f64) * (cycles as f64));
                let containment = P1 / (P1+P2+P3);
                vectors.sols_Ex[t][12] = containment;
            }
        }
    /*
    println!("vectors.sols_Ex[0] = {:?}", vectors.sols_Ex);
    println!("0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x");
    */

    }
    println!("Number of modes E_x = {}", counter);
    println!("k = {}", params.k);
    println!("k * n_0 = {}", params.k * params.n_0);
    println!("k * n_1 = {}", params.k * params.n_1);
    for elements in vectors.sols_Ex.iter() {
        if elements[0] > 0.0 {
            println!("---------------------------------");
            println!("p = {}, q = {}, a (micron) = {}, d (micron) = {},
            n_1 = {}, n_0 = {}", elements[0], elements[1], elements[2] * 1e6,
            elements[3] * 1e6, elements[5], elements[4]);
            println!("kx = {}", elements[6]);
            println!("ky = {}", elements[7]);
            println!("gamma_x = {}", elements[8]);
            println!("gamma_y = {}", elements[9]);
            println!("beta = {}", elements[10]);
            println!("n_eff = {}", n_eff);
            println!("cont_x = {}", elements[12]);
        }
    }
//end of block
//0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
}

////////////////////////////////RESET////////////////////////////////////////
vectors.k_x = [0.0; 4];
vectors.k_y = [0.0; 4];
vectors.gamma_x = [0.0; 4];
vectors.gamma_y = [0.0; 4];
vectors.beta = [[0.0; 4]; 4];

a = params.a;
d = params.d;
k = params.k;
n_0 = params.n_0;
n_1 = params.n_1;
stepsize = k * n_1 / (STEPS as f64);
P = params.p as usize;
Q = params.q as usize;
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
    for i in 0..STEPS {
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
    for i in 0..STEPS {
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
                        sum_x += E_y::Hx1(1.0, k_x, k_y, phi, psi, x, y)

                    }
                    sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                    sum_x = 0.0;
                }
            }
            let P1: f64 = c1 * (2.0 * d) * (sum_y).powi(2)/((N as f64) * (cycles as f64));
            sum_x = 0.0;
            sum_y = 0.0;
            let range_x  = Range::new(a, 3.0 * a);
            let range_y  = Range::new(-1.0 * d, d);
            let mut rng = rand::thread_rng();
            //Region 2, -a<=x<=a, -d<=y<=d
            for z in 0..cycles {
                for i in 0..N {
                    y = range_y.ind_sample(&mut rng);
                    for j in 0..N {
                        x = range_x.ind_sample(&mut rng);
                        sum_x += E_y::Hx2(1.0, k_x, k_y, gamma_x, phi, psi, x, y, a)

                    }
                    sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                    sum_x = 0.0;
                }
            }

            let P2: f64 = c2 * (2.0 * d) * (sum_y).powi(2)/((N as f64) * (cycles as f64));
            let range_x  = Range::new(-1.0 * a, a);
            let range_y  = Range::new(d, 3.0 * d);
            let mut rng = rand::thread_rng();
            sum_x = 0.0;
            sum_y = 0.0;
            //Region 3, -a<=x<=a, -d<=y<=d
            for z in 0..cycles {
                for i in 0..N {
                    y = range_y.ind_sample(&mut rng);
                    for j in 0..N {
                        x = range_x.ind_sample(&mut rng);
                        sum_x += E_y::Hx3(1.0, k_x, k_y, gamma_y, phi, psi, x, y, d)

                    }
                    sum_y += ((2.0 * a)/(N as f64)) * sum_x;
                    sum_x = 0.0;
                }
            }
            let P3: f64 = c3 * (2.0 * d) * (sum_y).powi(2)/((N as f64) * (cycles as f64));
            let containment = P1 / (P1+P2+P3);
            vectors.sols_Ey[t][12] = containment;
        }
    }
/*
println!("vectors.sols_Ex[0] = {:?}", vectors.sols_Ex);
println!("0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x");
*/

}
println!("==================================================");
println!("Number of modes E_y = {}", counter);
println!("k = {}", params.k);
println!("k * n_0 = {}", params.k * params.n_0);
println!("k * n_1 = {}", params.k * params.n_1);
for elements in vectors.sols_Ey.iter() {
    if elements[0] > 0.0 {
        println!("---------------------------------");
        println!("p = {}, q = {}, a (micron) = {}, d (micron) = {},
        n_1 = {}, n_0 = {}", elements[0], elements[1], elements[2] * 1e6,
        elements[3] * 1e6, elements[5], elements[4]);
        println!("kx = {}", elements[6]);
        println!("ky = {}", elements[7]);
        println!("gamma_x = {}", elements[8]);
        println!("gamma_y = {}", elements[9]);
        println!("beta = {}", elements[10]);
        println!("n_eff = {}", n_eff);
        println!("cont_x = {}", elements[12]);
    }
}
//end of block
//0=p, 1=q, 2=a, 3=d, 4=n_0, 5=n_1, 6=kx, 7=ky, 8=gamma_x, 9=gamma_y, 10=beta, 11=n_eff, 12=cont_x,
}
}
