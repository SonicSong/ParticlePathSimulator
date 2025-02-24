use std::fs::OpenOptions;
use std::ptr::read;
mod modules;
use mendeleev::Element;
use std::io::stdin;

// fn main() {
//     print!("Write the symbol of element: ");
//     let mut name_of_element = String::new();
//     stdin().read_line(&mut name_of_element).expect("Element symbol");
//
//     modules::bethe_formula::low_energies_calc(&*name_of_element);
// }
// Add function to calculate radiation damage

fn main() {
    // print!("{}", Element::H.symbol());
    // assert_eq!(Element::H.symbol(), "H");
    //
    // let mut name_of_chemical = String::new();
    // stdin().read_line(&mut name_of_chemical).expect("TODO: panic message");
    // print!("\n{}", name_of_chemical);

    let Some((atom_density, atom_number, mass_number)) = modules::periodic::look_up_element("W") else { todo!()};
    println!("{:?}", atom_density);
    println!("{}", atom_number);
    println!("{}", mass_number);


}