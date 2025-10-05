mod modules;
use modules::bethe_formula::*;
// use std::io::stdin;
// use mendeleev::Element;

fn main() {
    // print!("Write the symbol of element: ");
    // let mut name_of_element = String::new();
    // stdin().read_line(&mut name_of_element).expect("Element symbol");
    // name_of_element.pop();
    let name_of_element = "Si";
    let name_of_incident_particle = "Proto";
    modules::bethe_formula::low_energies_calc(&*name_of_element, name_of_incident_particle);
}

// Add function to calculate radiation damage