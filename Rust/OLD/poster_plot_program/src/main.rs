
extern crate gnuplot;

use gnuplot::*;
use std::{f64, i64};
use std::f64::consts;
use std::error::Error;
use std::str::FromStr;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;


fn main() {
    let mut z_data: Vec<f64> = vec![];
    let a: f64 = 5.0; //Width
    let d: f64 = 5.0; //Depth
    let npoints: usize = 5;
    let mut x: f64 = 0.0;
    let mut y: f64 = 0.0;

    z_data = vec![
    0.8939560708,
    0.8094055566,
    0.7204913537,
    0.6450471613,
    0.5807280139,
    0.7964471515,
    0.8742907304,
    0.8359691175,
    0.7909914637,
    0.7486411865,
    0.7084386308,
    0.8333166779,
    0.8704103903,
    0.8459685674,
    0.8167191899,
    0.6341647003,
    0.7880459672,
    0.8449732758,
    0.8676959586,
    0.8491461847,
    0.5738801681,
    0.7435197876,
    0.8154417988,
    0.8486742435,
    0.8673842536,
    ];

    for i in 0..5 {
        x = 1.0 + (i as f64);
        for j in 0..5 {
            y = 1.0 + (j as f64);

        }
    }

    //Plot for E_x
    let mut figure = Figure::new();
    figure.axes3d()
    .set_title("Containment", &[])
    .surface(z_data.iter(), npoints, npoints, Some((1.0, 1.0, a, d)), &[])
    .set_z_label("Containment", &[])
    .set_x_label("Width", &[])
    .set_y_label("Depth",&[])
    .set_view_map();
    //.set_view(45.0, 45.0); //3D display
    //figure.set_terminal("pngcairo", "mode_field.png");    let mut fg3 = Figure::new();
    figure.show();
}
