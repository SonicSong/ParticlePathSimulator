use std::any::TypeId;
use std::f64::consts;
use mendeleev::{Element, GramPerCubicCentimeter};
use std::convert::TryFrom;
use std::convert;
use modules::periodic_lookup;

use crate::modules;
//TODO: Verify data and make sure calculation is correct

const PI_NUMBER: f64 = consts::PI;
const AVOGADRO_CONST: f64 = 6.02214076e23; //N↓A
const LIGHT_SPEED: f64 = 2.998e8; // m/s speed of light
const ELECTRON_MASS: f64 = 9.1093837e-31; // MeV/c² or kg - electron mass
const ELECTRON_CHARGE: f64 = 1.602176634e-19; // C - elementary charge
const VACUUM_PERMITTIVITY:f64 = 8.8541878188e-12; // F/m - electric constant
const MOLAR_MASS: f64 = 0.18384;
const BETA: f64 = 0.7;
const I: f64 = 727.0;
const PLANCK_CONST: f64 = 6.62607015e-34; // kg * m^2 * s^-1 - planck constant

// Zdefiniowanie zmiennych i stałych

// Materiał (krzem):
// Liczba atomowa Z=14
// Masa atomowa A=28.0855g/mol
// Gęstość ρ=2.33g/cm3
// Średni potencjał wzbudzenia I=173eV
// Prędkość elektronu:
// Załóżmy, że elektron ma energię kinetyczną T=10MeV
// Stałe:
// Masa elektronu me=0.511MeV/c2
// Prędkość światła c=3×108m/s
// Ładunek elementarny e=1.602×10−19
// Liczba Avogadro NA=6.022×1023mol−1

pub fn low_energies_calc(name_of_element: &str) {
    println!("\n Element: {}", name_of_element);
    println!("\n fpi_na_zp: {:.20}", fpi_na_zp(name_of_element));
    println!("\n etwo_by_fpi: {:.20}", etwo_by_fpi());
    println!("\n ztwo_by_betatwo: {:.20}", ztwo_by_betatwo());
    println!("\n twom_e_ctwo: {:.20}", twom_e_ctwo_betatwo_tmax());
    let de_dx: f64 = fpi_na_zp(name_of_element) * etwo_by_fpi() * ztwo_by_betatwo() * twom_e_ctwo_betatwo_tmax();
    println!("dE/dx: {:.20} J/m", de_dx);
}


// 4pi N_A Z_p / A m_m m_e c²
fn fpi_na_zp(name_of_element: &str) -> f64 {
    if let Some((atom_density, atom_number, mass_number)) = periodic_lookup::look_up_element(name_of_element) {
        let top_calc: f64 = 4.0 * PI_NUMBER * AVOGADRO_CONST * atom_number * atom_density;
        let bottom_calc: f64 = mass_number * MOLAR_MASS * ELECTRON_MASS * (0.7 * LIGHT_SPEED).powi(2);
        top_calc / bottom_calc
    } else {
        eprintln!("Element {} not found or density is unavailable.", name_of_element);
        0.0
    }
}

// e² / 4pi e_0
fn etwo_by_fpi() -> f64 {
    let etwo: f64 = (ELECTRON_CHARGE.powi(2));
    let bot: f64 = 4.0 * PI_NUMBER * VACUUM_PERMITTIVITY;
    let top_by_bot: f64 = (etwo/bot)*(etwo/bot);
    top_by_bot
}

// z² / Beta²
fn ztwo_by_betatwo() -> f64 {
    let top: f64 = 1.0 * 1.0; // z^2 dla protonu
    let bottom: f64 = BETA.powi(2); // beta^2
    let result: f64 = top / bottom; // z^2 / beta^2
    result
}

fn twom_e_ctwo_betatwo_tmax() -> f64 {
    let top: f64 = 2.0 * ELECTRON_MASS * (LIGHT_SPEED.powi(2)) * (BETA.powi(2)) * tmax_calculation();
    let bottom: f64 = (1.0 - (BETA.powi(2))) * (I.powi(2));

    if bottom == 0.0 {
        panic!("Division by zero in twom_e_ctwo_betatwo_tmax");
    }

    let ratio: f64 = top / bottom;
    if ratio <= 0.0 {
        panic!("Attempting to take log of non-positive value in twom_e_ctwo_betatwo_tmax");
    }

    let result: f64 = 0.5 * f64::ln(ratio) - (BETA.powi(2)) - (delta() / 2.0);
    // println!("Top: {}, Bottom: {}, Ratio: {}, Result: {}", top, bottom, ratio, result);
    result
}


fn gamma() -> f64 {
    let gamma: f64 = 1.0 / (1.0 - (BETA.powi(2))).sqrt();
    gamma
}

fn delta() -> f64 {
    let hwp_i: f64 = f64::ln((PLANCK_CONST * 1.0) / I);
    let beta_gamma:f64 = f64::ln(BETA * gamma());
    let result: f64 = hwp_i + beta_gamma - 0.5;
    result
}

fn wp_plasma() -> f64 {
    let top_by_bottom: f64 = 1.0;
    let result: f64 = ELECTRON_CHARGE * top_by_bottom.sqrt();
    result
}

fn tmax_calculation() -> f64 {
    let m : f64 = 1.672e-27;
    let top: f64 = 2.0 * ELECTRON_MASS * (LIGHT_SPEED.powi(2)) * (BETA.powi(2)) * (gamma().powi(2));
    let bot: f64 = 1.0 + (2.0 * gamma() * ELECTRON_MASS) / m;
    let result: f64 = top / bot;
    result
}