// extern crate uom;

// use uom::si::{mole, molar_mass};
// use uom::num::One;
// use uom::si::f64::*;
// use uom::si::amount_of_substance;
// use uom::si::amount_of_substance::mole;
// use uom::si:
use std::f64::consts;

// const MOLENUM : MoleNumber = AmountOfSubstance::new::<mole>;
pub const PI_NUMBER: f64 = consts::PI;
pub const AVOGDRO_CONST: f64 = 6.02214076e23; //N↓A
pub const LIGHT_SPEED: f64 = 1.08e9;


fn main() {
    fpi_na_zp();
}

fn low_energies_calc() {

}

fn fpi_na_zp() {
    const TUNGSTEN : (f64, f64, f64) = (19.3e3, 74.0, 183.84);
    let (atom_density, atom_number, mass_number) = TUNGSTEN;
    let top_calc : f64 = 4.0 * PI_NUMBER * AVOGDRO_CONST * atom_number * atom_density;
    let bottom_calc : f64 = mass_number * 1.0 * (LIGHT_SPEED * LIGHT_SPEED);
}