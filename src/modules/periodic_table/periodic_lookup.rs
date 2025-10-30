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
        "B" => Some((
            bethe_formula::precise("2.370"),     // density
            bethe_formula::precise("5"),       // Z
            bethe_formula::precise("10.81"),   // A
        )),
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

pub fn look_up_element_mass(symbol_name: &str) -> Option<Float> {
    match symbol_name {
        "B" => Some(bethe_formula::precise("10.817")),
        "Si" => Some(bethe_formula::precise("28.0855")), // A
        "C" => Some(bethe_formula::precise("12.0107")),
        "Al" => Some(bethe_formula::precise("26.9815")),
        "Au" => Some(bethe_formula::precise("196.966569")),
        "Cu" => Some(bethe_formula::precise("63.546")),
        _ => None, // fallback
    }
}

pub fn look_up_element_k_and_a(symbol_name: &str) -> Option<(Float, Float)> {
    match symbol_name {
        "B" => Some((
            bethe_formula::precise("0.56224"),  // a
            bethe_formula::precise("2.4512"),   // k =m_s
        )),
        "Si" => Some((
            bethe_formula::precise("0.14921"),
            bethe_formula::precise("3.2546"),
        )),
        "C" => Some((
            bethe_formula::precise("0.20240"),
            bethe_formula::precise("3.0036"),
        )),
        "Al" => Some((
            bethe_formula::precise("0.08024"),
            bethe_formula::precise("3.6345"),
        )),
        "Au" => Some((
            bethe_formula::precise("0.09756"),
            bethe_formula::precise("3.1101"),
        )),
        "Cu" => Some((
            bethe_formula::precise("0.14339"),
            bethe_formula::precise("2.9044"),
        )),
        _ => None
    }
}

pub fn look_up_element_muon(symbol_name: &str) -> Option<(Float, Float, Float, Float)> {
    match symbol_name {
        "B" => Some((

            bethe_formula::precise("0.0305"),   // x_0
            bethe_formula::precise("1.9688"),   // x_1
            bethe_formula::precise("2.8477"),   // C
            bethe_formula::precise("0.14")      // δ_0
        )),
        "Si" => Some((
            bethe_formula::precise("0.2015"),
            bethe_formula::precise("2.8716"),
            bethe_formula::precise("4.4355"),
            bethe_formula::precise("0.14")
        )),
        "C" => Some((
            bethe_formula::precise("-0.0351"),
            bethe_formula::precise("2.4860"),
            bethe_formula::precise("2.9925"),
            bethe_formula::precise("0.10")
        )),
        "Al" => Some((
            bethe_formula::precise("0.1708"),
            bethe_formula::precise("3.0127"),
            bethe_formula::precise("4.2395"),
            bethe_formula::precise("0.12")
        )),
        "Au" => Some((
            bethe_formula::precise("0.2021"),
            bethe_formula::precise("3.6979"),
            bethe_formula::precise("5.5747"),
            bethe_formula::precise("0.14")
        )),
        "Cu" => Some((
            bethe_formula::precise("-0.0254"),
            bethe_formula::precise("3.2492"),
            bethe_formula::precise("4.4190"),
            bethe_formula::precise("0.08")
        )),
        _ => None
    }
}