// extern crate uom;

// use uom::si::{mole, molar_mass};
// use uom::num::One;
// use uom::si::f64::*;
// use uom::si::amount_of_substance;
// use uom::si::amount_of_substance::mole;

use std::f64::consts;

// const MOLENUM : MoleNumber = AmountOfSubstance::new::<mole>;
pub const PI_NUMBER: f64 = consts::PI;
pub const AVOGDRO_CONST: f64 = 6.02214076e23; //N↓A
pub const LIGHT_SPEED: f64 = 1.08e9;
pub const ELECTRON_MASS: f64 = 9.1093837e-31;
pub const ELECTRON_CHARGE: f64 = 1.602176634e-19;
pub const VACUUM_PERMITTIVITY:f64 = 8.8541878188e-12;


fn main() {
    // println!("{}", PI_NUMBER);
    low_energies_calc();
}

fn low_energies_calc() {
    println!("{}", fpi_na_zp());
    println!("{}", etwo_by_fpi());
    let de_dx: f64 = fpi_na_zp() * etwo_by_fpi();
    println!("dE/dx: {}", de_dx);
}

fn fpi_na_zp() -> f64 {
    const MOLAR_MASS: f64 = 79.2967;
    const TUNGSTEN : (f64, f64, f64) = (19.3e3, 74.0, 183.84);
    let (atom_density, atom_number, mass_number) = TUNGSTEN;
    let top_calc : f64 = 4.0 * PI_NUMBER * AVOGDRO_CONST * atom_number * atom_density;
    let bottom_calc : f64 = mass_number * MOLAR_MASS * ELECTRON_MASS * (LIGHT_SPEED * LIGHT_SPEED);
    let top_by_bottom: f64 = top_calc / bottom_calc;
    top_by_bottom
}

fn etwo_by_fpi() -> f64 {
    let top: f64 = (ELECTRON_CHARGE * ELECTRON_CHARGE);
    let bot: f64 = 4.0 * PI_NUMBER * VACUUM_PERMITTIVITY;
    let top_by_bot: f64 = (top/bot)*(top/bot);
    top_by_bot
}