use rug::Float;
use modules::bethe_formula;
use crate::modules;
pub fn mean_excitation_energy(symbol_name: &str) -> Option<Float> {
    match symbol_name {
        "B" => Some(
            bethe_formula::precise("76.0")
        ),
        "Si" => Some(
            bethe_formula::precise("173.0")
        ),
        "C" => Some(
            bethe_formula::precise("78.0")
        ),
        "Al" => Some(
            bethe_formula::precise("166.0")
        ),
        "Au" => Some(
            bethe_formula::precise("790.0")
        ),
        "Cu" => Some(
            bethe_formula::precise("322.0")
        ),
        _ => None,
    }
}