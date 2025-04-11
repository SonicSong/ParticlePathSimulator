use std::f64::consts;
use mendeleev::{Element};
use std::iter::Iterator;
use modules::periodic_lookup;
use num::traits::real::Real;
use crate::modules::CustomProperties;

// https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf
//TODO: Rewrite the entire formula to match for the universal version found on (https://pdg.lbl.gov/rpp/encoders/pdg_note_1901.pdf)

use crate::modules;
//TODO: Verify data and make sure calculation is correct
//TODO: Make corrections for the data formats that use either CSG or SI or Both

//Constants https://pdg.lbl.gov/2024/reviews/rpp2024-rev-phys-constants.pdf

const PI_NUMBER: f64 = consts::PI;
const LIGHT_SPEED: f64 = 299792458.0; // m/s speed of light
const ELECTRON_MASS: f64 = 9.1093837e-31; // kg - electron mass (m_e)
const PLANCK_CONST: f64 = 6.62607015e-34;
const AVOGADRO_CONST: f64 = 6.02214076e23; // mol ^ -1

pub fn low_energies_calc(name_of_element: &str, name_of_incident_particle: &str) {
    println!("\n Element: {}", name_of_element);

    let element_exci_energy = match Element::iter().find(|e|
        e.symbol().eq_ignore_ascii_case(name_of_element) ||
            e.name().eq_ignore_ascii_case(name_of_element)) {
        Some(e) => e,
        None => panic!(),
    };

    let energy: f64 = m_e_cpowit();
    let mut de_dx_array: Vec<f64> = vec![];
    let mut velocity: Vec<f64> = vec![];
    // TODO: Loop for going through various values for BETA and GAMMA. For example Beta = (0.1*C/C) 0.1*C can be considered V

    for i in 1i32..1000i32 {
        let i_t = f64::from(i) * 0.001;
        let mut beta: f64 = (i_t * LIGHT_SPEED) / LIGHT_SPEED;

        let denominator = 1.0 - ((i_t*LIGHT_SPEED)/LIGHT_SPEED).powi(2);
        let gamma: f64;
        const VELOCITY_CUTOFF: f64 = 0.999;
        if denominator > 0.0  && i_t <= VELOCITY_CUTOFF {
            gamma = 1.0 / denominator.sqrt();
            velocity.push(gamma);
        } else {
            gamma = 0.0;
            // panic!();
        }
        let de_dx: f64 = k_z_two_z_a_1_b_two(name_of_element, beta, name_of_incident_particle) *
                (0.5 * (twom_e_ctwo_btwo_dtwo_w(beta, gamma, energy, element_exci_energy, name_of_incident_particle).ln() - beta.powi(2) -
                (1.0 / 2.0))); // Missing for now value of δ(βγ) which is currently 1.0
        de_dx_array.push(de_dx);
    }
    println!("{:?}", de_dx_array);
    // println!("dE/dx:  J/m");
}

// const m_e_cpowit: f64 = ELECTRON_MASS * LIGHT_SPEED.powi(2);
fn m_e_cpowit() -> f64 {
    let result = ELECTRON_MASS * LIGHT_SPEED.powi(2);
    // let result: f64 = 0.51099895000;
    result
}

fn density_effect_correction(beta: f64, gamma: f64) { // -> f64 {
    // δ(βγ)/2 → ln(ℏωp/I) + ln βγ − 1/2

}

fn k_z_two_z_a_1_b_two(name_of_element: &str, beta: f64, name_of_incident_particle: &str) -> f64 {
    //K = 4π * N * A* r^2_e * m_e * c^2     0.307075 MeV mol^−1 cm^2
    //z = -1 for electron
    //z = +1 for proton
    //z = +5 for Boron
    let mut z_inci;

    if let Some((atom_density, atom_number, mass_number)) = periodic_lookup::look_up_element(name_of_element) {
        if (name_of_incident_particle == "Ele") {
            z_inci = -1.0_f64;
        } else if (name_of_incident_particle == "Proto") {
            z_inci = 1.0_f64;
        } else if (name_of_element != "Ele" || name_of_element != "Proto") {
            z_inci = atom_number;
        } else {
            // Default electron
            z_inci = -1.0_f64
        }
        let k_z_two: f64 = 0.307075 * z_inci.powi(2);
        let z_by_a: f64 = atom_number / mass_number;
        let one_by_beta: f64 = 1.0 / beta.powi(2);
        let result: f64 = k_z_two * z_by_a * one_by_beta;
        result
    } else {
        eprintln!("Element {} not found or density is unavailable.", name_of_element);
        0.0
    }
}

fn twom_e_ctwo_btwo_dtwo_w(beta: f64, gamma: f64, m_e_cpowit: f64, element_exci_energy: Element, name_of_incident_particle: &str) -> f64 {
    let twom_e_ctwo: f64 = 2.0 * m_e_cpowit * beta.powi(2) * gamma.powi(2) * wmax(beta, gamma, m_e_cpowit, name_of_incident_particle);

    // Get mean excitation energy
    match element_exci_energy.custom_mean_excitation_energy() {
        Some(i_powi_two) => twom_e_ctwo / i_powi_two.powi(2),
        None => 0.0,  // Return 0.0 for elements without mean excitation energy
    }
}

fn calculate_incident_particle_mass(name_of_incident_particle: &str) -> f64{

    let uni_amu: f64 = 931.4941024228; // MeV
    if (name_of_incident_particle == "Ele" || name_of_incident_particle == "Proto") {
        if (name_of_incident_particle == "Ele") {
            let result: f64 = 0.51099895000; // MeV/c2
            result
        } else {
            let result: f64 = 938.27208816; // MeV/c2
            result
        }
    } else {
        if let Some((mass_number)) = periodic_lookup::look_up_element_weight(name_of_incident_particle) {
            let result: f64 = mass_number * uni_amu;
            result
        } else {
            eprintln!("Element {} not found or density is unavailable.", name_of_incident_particle);
            0.0
        }
    }
}

fn wmax(beta: f64, gamma: f64, m_e_cpowit: f64, name_of_incident_particle: &str) -> f64 {
    //TODO: Not compatible with electrons. Need to do find the calculation for Wmax that doesn't divide electron mass by electron mass.
    if (name_of_incident_particle == "Ele") {
        let result: f64 = 1.0;
        result
    } else {
        let two_m_e_ctwo: f64 = 2.0 * m_e_cpowit * beta.powi(2) * gamma.powi(2);
        let one_two_gamma_m_e: f64 = 1.0 + ((2.0 * gamma * 0.51099895000) / calculate_incident_particle_mass(name_of_incident_particle))
            + (0.51099895000 / calculate_incident_particle_mass(name_of_incident_particle)).powi(2);
        //TODO: Verify if the incident particle mass is calculated correctly.
        let result: f64 = two_m_e_ctwo/one_two_gamma_m_e;
        result
    }
}