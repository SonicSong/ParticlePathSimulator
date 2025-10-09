mod modules;

use std::sync::atomic::Ordering;
use modules::bethe_formula;
// use std::io::stdin;
// use mendeleev::Element;
use modules::atomic_vars;
use rug::Float;

fn main() {

    //TODO: Change everything from f64 to "something" that will delay IEEE754 rounding error for the last possible moment

    let PRECISION_BITS: u32 = crate::modules::atomic_vars::PRECISION.load(Ordering::Relaxed) as u32;

    let TEST_PRECISION = Float::with_val(200, 0.001);
    println!("PRECISION: 0.0{}", TEST_PRECISION);
    println!("PRECISION: {:.64} <- f64", 0.001f64);


    // print!("Write the symbol of element: ");
    // let mut name_of_element = String::new();
    // stdin().read_line(&mut name_of_element).expect("Element symbol");
    // name_of_element.pop();
    // let name_of_element = "Si";
    // let name_of_incident_particle = "Proto";
    // modules::bethe_formula::low_energies_calc(&*name_of_element, name_of_incident_particle);
}

// Add function to calculate radiation damage