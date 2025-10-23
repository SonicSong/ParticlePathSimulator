mod modules;

use std::sync::atomic::Ordering;
// use modules::bethe_formula;
// use std::io::stdin;
// use mendeleev::Element;
use modules::atomic_vars;
use rug::{Assign, Float};
use crate::modules::atomic_vars::PRECISION;
use std::any::type_name;
use rug::float::Constant;

fn print_type_of<T>(_: &T) {
    println!("{} <- TYPE OF DATA ", type_name::<T>());
}

fn main() {

    //OLD: Change everything from f64 to "something" that will delay IEEE754 rounding error for the last possible moment
    //Switching to rug::float "fixes" problem of IEEE754 by having better precision while working on float values

    let PRECISION_BITS: u32 = PRECISION.load(Ordering::Relaxed) as u32;

    let TEST_PRECISION = Float::with_val(256, 0.001);

    let mut PI_NUMB = Float::new(PRECISION_BITS);
    PI_NUMB.assign(Constant::Pi);
    println!("PI = {}", PI_NUMB);
    print_type_of(&PI_NUMB);

    let mut TEST_PRECISION_S = Float::new(PRECISION_BITS);
    let mut TEST_PRECISION_S_t = Float::new(PRECISION_BITS);

    let mut SUM = Float::new(PRECISION_BITS);
    TEST_PRECISION_S.assign(Float::parse("0.001").expect("Invalid float string"));
    TEST_PRECISION_S_t.assign(Float::parse("2.000").expect("Invalid float string"));
    SUM.assign(&TEST_PRECISION_S + &TEST_PRECISION_S_t);
    println!("TEST ADDITION {:?}", SUM);
    print_type_of(&SUM);
    println!("PRECISION: {} <- RUG NO STRING", TEST_PRECISION);
    print_type_of(&TEST_PRECISION);
    println!("Precision str: {} <- RUG W/STRING", TEST_PRECISION_S);
    print_type_of(&TEST_PRECISION_S);
    println!("PRECISION: {:.64} <- f64", 0.001f64);
    print_type_of(&0.001f64);


    // print!("Write the symbol of element: ");
    // let mut name_of_element = String::new();
    // stdin().read_line(&mut name_of_element).expect("Element symbol");
    // name_of_element.pop();
    // let name_of_element = "Si";
    // let name_of_incident_particle = "Al";
    // modules::bethe_formula::stopping_power_intermediate_energies(&*name_of_element, name_of_incident_particle);
}

// Add function to calculate radiation damage