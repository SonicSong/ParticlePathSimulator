use std::f64::consts;
use std::option;
use std::convert;
use mendeleev::{Element, GramPerCubicCentimeter};

pub fn look_up_element(symbol_name: &str) -> Option<(f64, f64, f64)> {
    Element::iter()
        .find(|e| e.symbol().eq_ignore_ascii_case(symbol_name))
        .map(|element| {
            (
                element.density().map(|gpcc| gpcc.0).unwrap_or(0.0),
                element.atomic_number() as f64,
                element.atomic_weight().into()
            )
        })
}

pub fn look_up_element_weight(symbol_name: &str) -> Option<f64> {
    Element::iter()
        .find(|e|e.symbol().eq_ignore_ascii_case(symbol_name))
        .map(|element| {
            element.atomic_weight().into()
        })
}

// Add function to calculate density of the material

pub fn multiple_material_density() {
    // p = m/V equation for the material density

}

pub fn periodic_table() {}