mod modules;
use std::io::stdin;
use mendeleev::Element;

fn main() {
    // print!("Write the symbol of element: ");
    // let mut name_of_element = String::new();
    // stdin().read_line(&mut name_of_element).expect("Element symbol");
    // name_of_element.pop();
    let name_of_element = "Si";
    // if let Some(energy) = name_of_element.custom_mean_excitation_energy() {
    //     println!("Mean excitation energy: {:.2} eV", energy);
    // }


    modules::bethe_formula::low_energies_calc(&*name_of_element);
}

// Add function to calculate radiation damage