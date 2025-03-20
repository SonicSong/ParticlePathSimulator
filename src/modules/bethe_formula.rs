use std::any::TypeId;
use std::f64::consts;
use mendeleev::{Element, GramPerCubicCentimeter};
use std::convert::TryFrom;
use std::convert;
use std::iter::Iterator;
use modules::periodic_lookup;
use num::range_step;
use num::traits::real::Real;
// https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf
//TODO: Rewrite the entire formula to match for the universal version found on (https://pdg.lbl.gov/rpp/encoders/pdg_note_1901.pdf)

use crate::modules;
//TODO: Verify data and make sure calculation is correct

const PI_NUMBER: f64 = consts::PI;
const LIGHT_SPEED: f64 = 2.998e8; // m/s speed of light
const ELECTRON_MASS: f64 = 9.1093837e-31; // MeV/c² or kg - electron mass
const ELECTRON_CHARGE: f64 = 1.602176634e-19; // C - elementary charge
const VACUUM_PERMITTIVITY:f64 = 8.8541878188e-12; // F/m - electric constant
const MOLAR_MASS: f64 = 0.18384;
const I: f64 = 727.0;
const PLANCK_CONST: f64 = 6.62607015e-34;


const AVOGADRO_CONST: f64 = 6.02214076e23; // mol ^ -1

// const BETA: f64 = 0.7; // VELOCITY * LIGHT_SPEED (Velocity is stored as a vector) (Beta = V/c Velocity by light speed) v/c of incident particle
// const VELOCITY: Vec<f64> = range_step(0.1 * LIGHT_SPEED, LIGHT_SPEED, 100000.0).collect();

// const GAMMA: f64 = 1.0/(1.0 - BETA.powi(2)).sqrt();
//v <- seq(0.1*c,c,100000)
// GAMMA is going to be changed to work in function as a let and not as a const

pub fn low_energies_calc(name_of_element: &str) {
    println!("\n Element: {}", name_of_element);

    let energy: f64 = m_e_cpowit();
    let mut de_dx_array: Vec<f64> = vec![];
    let mut velocity: Vec<f64> = vec![];
    // TODO: Loop for going through various values for BETA and GAMMA. For example Beta = (0.1*C/C) 0.1*C can be considered V
    for i in 1i32..100i32 {
        let i = f64::from(i) * 0.1;
        let mut beta: f64 = (i * LIGHT_SPEED) / LIGHT_SPEED;
        let mut gamma: f64 = 1.0/(1.0 - beta.powi(2)).sqrt();
        let de_dx: f64 = k_z_two_z_a_1_b_two(name_of_element, 1.0) * twom_e_ctwo_btwo_dtwo_w(name_of_element, beta, gamma, energy);
        de_dx_array.push(de_dx);
    }
    println!("{:?}", de_dx_array);
    // println!("dE/dx:  J/m");
}

// const m_e_cpowit: f64 = ELECTRON_MASS * LIGHT_SPEED.powi(2);
fn m_e_cpowit() -> f64 {
    let result = ELECTRON_MASS * LIGHT_SPEED.powi(2);
    result
}

fn k_z_two_z_a_1_b_two(name_of_element: &str, beta: f64) -> f64 {
    //K = 4π * N * A* r^2_e * m_e * c^2     0.307075 MeV mol^−1 cm^2
    //z = -1 for electron
    //z = +1 for proton
    let k_z_two: f64 = 0.307075 * (1.0_f64).powi(2);
    if let Some((atom_density, atom_number, mass_number)) = periodic_lookup::look_up_element(name_of_element) {
        let z_by_a: f64 = atom_number / mass_number;
        let one_by_beta: f64 = 1.0 / beta.powi(2);
        let result: f64 = k_z_two * z_by_a * one_by_beta;
        result
    } else {
        eprintln!("Element {} not found or density is unavailable.", name_of_element);
        0.0
    }
}

fn twom_e_ctwo_btwo_dtwo_w(name_of_element: &str, beta: f64, gamma: f64, m_e_cpowit: f64) -> f64 {
    let twom_e_ctwo: f64 = 2.0 * m_e_cpowit * beta.powi(2) * gamma.powi(2) * wmax(beta, gamma, m_e_cpowit);
    let i_pow_two: f64 = 2.0; // Missing M^2 - mean excitation energy
    let result: f64 = twom_e_ctwo / i_pow_two;
    result
}

fn beta() -> f64 {

    let result: f64 = 1.0;
    result
}

fn wmax(beta: f64, gamma: f64, m_e_cpowit: f64) -> f64 {
    let two_m_e_ctwo: f64 = 2.0 * m_e_cpowit * beta.powi(2) * gamma.powi(2);
    let one_two_gamma_m_e: f64 = 1.0 + 2.0 * gamma * (ELECTRON_MASS / 1.5) + (ELECTRON_MASS / 1.5).powi(2);
    //^ Missing M - incident particle mass (Temporarly as 1.5)

    let result: f64 = two_m_e_ctwo/one_two_gamma_m_e;
    result
}