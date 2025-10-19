use std::f64::consts;
use std::option;
use std::convert;
use mendeleev::{Element, GramPerCubicCentimeter};
use modules::bethe_formula;
use crate::modules;
use rug::Float;

//TODO: Phase out mendeleev for ANPM.
//Reason for phasing out mendeleev crate. It lacks a lot of information and is mainly for general use. Not advanced or custom properties.

// pub fn look_up_element(symbol_name: &str) -> Option<(f64, f64, f64)> {
//     Element::iter()
//         .find(|e| e.symbol().eq_ignore_ascii_case(symbol_name))
//         .map(|element| {
//             (
//                 element.density().map(|gpcc| gpcc.0).unwrap_or(0.0),
//                 element.atomic_number() as f64,
//                 element.atomic_weight().into()
//             )
//         })
// }

// For now using only few elements hardcoded for testing purposes. And later will expand the database or use external crate.

pub fn look_up_element(symbol_name: &str) -> Option<(Float, Float, Float)> {
    match symbol_name {
        "Si" => Some((
            bethe_formula::precise("2.329"),   // density
            bethe_formula::precise("14"),      // Z
            bethe_formula::precise("28.0855"), // A
        )),
        "C" => Some((
            bethe_formula::precise("2.267"),
            bethe_formula::precise("6"),
            bethe_formula::precise("12.0107"),
        )),
        "Al" => Some((
            bethe_formula::precise("2.70"),
            bethe_formula::precise("13"),
            bethe_formula::precise("26.9815"),
        )),
        "Au" => Some((
            bethe_formula::precise("19.32"),
            bethe_formula::precise("79"),
            bethe_formula::precise("196.966569"),
        )),
        "Cu" => Some((
            bethe_formula::precise("8.96"),
            bethe_formula::precise("29"),
            bethe_formula::precise("63.546"),
        )),
        _ => None, // fallback
    }
}

pub fn look_up_element_weight(symbol_name: &str) -> Option<Float> {
    match symbol_name {
        "Si" => Some(bethe_formula::precise("28.0855")), // A
        "C" => Some(bethe_formula::precise("12.0107")),
        _ => None, // fallback
    }
}

// Add function to calculate density of the material

pub fn multiple_material_density() {
    // p = m/V equation for the material density

}

pub fn periodic_table() {}