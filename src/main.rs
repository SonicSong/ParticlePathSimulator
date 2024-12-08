extern crate uom;

// use uom::si::{mole, molar_mass};
use uom::num::One;
use uom::si::f64::*;
use uom::si::amount_of_substance;
use uom::si::amount_of_substance::mole;

// const MOLENUM : MoleNumber = AmountOfSubstance::new::<mole>;


fn main() {
    const AVOGDRO_CONST: f64 = 6.02214076e23; //N↓A
    const ATOM_NUM : f64 = 12.21;

    println!("{}", AVOGDRO_CONST);
}
