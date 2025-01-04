use std::f64::consts;

const PI_NUMBER: f64 = consts::PI;
const AVOGDRO_CONST: f64 = 6.02214076e23; //N↓A
const LIGHT_SPEED: f64 = 2.998e8; // m/s speed of light
const ELECTRON_MASS: f64 = 9.1093837e-31; // MeV/c² or kg - electron mass
const ELECTRON_CHARGE: f64 = 1.602176634e-19; // C - elementary charge
const VACUUM_PERMITTIVITY:f64 = 8.8541878188e-12; // F/m - electric constant
const MOLAR_MASS: f64 = 0.18384;
const BETA: f64 = 0.7;

pub fn low_energies_calc() {
    println!("{}", fpi_na_zp());
    println!("{}", etwo_by_fpi());
    println!("{}", ztwo_by_betatwo());
    let de_dx: f64 = fpi_na_zp() * etwo_by_fpi() * ztwo_by_betatwo();
    println!("dE/dx: {}", de_dx);
}

// 4pi N_A Z_p / A m_m m_e c²
fn fpi_na_zp() -> f64 {
    const TUNGSTEN : (f64, f64, f64) = (19.3e-3, 74.0, 183.84);
    let (atom_density, atom_number, mass_number) = TUNGSTEN;
    let top_calc: f64 = 4.0 * PI_NUMBER * AVOGDRO_CONST * atom_number * atom_density;
    let bottom_calc : f64 = mass_number * MOLAR_MASS * ELECTRON_MASS * ((0.7 * LIGHT_SPEED) * (0.7 * LIGHT_SPEED));
    let top_by_bottom: f64 = top_calc / bottom_calc;
    top_by_bottom
}

// e² / 4pi e_0
fn etwo_by_fpi() -> f64 {
    let etwo: f64 = (ELECTRON_CHARGE * ELECTRON_CHARGE);
    let bot: f64 = 4.0 * PI_NUMBER * VACUUM_PERMITTIVITY;
    let top_by_bot: f64 = (etwo/bot)*(etwo/bot);
    top_by_bot
}

// z² / Beta²
fn ztwo_by_betatwo() -> f64 {
    let top: f64 = 1.0 * 1.0; // z^2 dla protonu
    let bottom: f64 = 0.7 * 0.7; // beta^2
    let result: f64 = top / bottom;
    result
}

fn twom_e_ctwo_betatwo_tmax() -> f64 {
    let top: f64 = 2.0 * ELECTRON_MASS * (LIGHT_SPEED * LIGHT_SPEED) * (BETA * BETA) * tmax_calculation();
    let bottom: f64 = 0.7 * 0.7;
    let result: f64 = top / bottom;
    result
}

fn tmax_calculation() -> f64 {
    let yotta: f64 = 1.0 / (1.0 - (BETA * BETA).sqrt());
    let top: f64 = 2.0 * ELECTRON_MASS * (LIGHT_SPEED * LIGHT_SPEED) * (BETA * BETA) * yotta;
}