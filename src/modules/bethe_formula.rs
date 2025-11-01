use std::f64::consts;
use std::iter::Iterator;
use modules::periodic_table::periodic_lookup;
use num::traits::real::Real;
use std::sync::atomic::Ordering;
use modules::atomic_vars;

//DONE: Replace all f64 to rug::Float for better precision during calculations
//TODO: Set global variable that can be changed in settings on runtime that defines the precision of Floats for "low-end machines"
// (Not sure how memory usage will spike so will start from small values)
use rug::{Assign, Float};
use rug::float::Constant;
use rug::ops::Pow;
use crate::modules::periodic_table::mean_excitation_energies::*;

// https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf
//TODO: Rewrite the entire formula to match for the universal version found on (https://pdg.lbl.gov/rpp/encoders/pdg_note_1901.pdf)


// For the ease of mine I simply listed here all the constants and variables used in the Bethe formula calculations
// Taken from https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf
/*
    Constants:
    m_e*c^2 = 0.51099895000 MeV electron mass × c^2
    r_e = 2.8179403227 fm classical electron radius (fm unit is femtometer = 10^-15 meter)
    α = 1/137.035999139 fine-structure constant
    N_A = 6.02214076×10^23 mol^−1 Avogadro constant

    Values:
    ρ (rho) = density of the absorber in g cm^−3
    x mass per unit area in g cm^−2
    M = incident particle mass in MeV/c^2
    E = incident particle energy γ*M*c^2 in MeV
    T = kinetic energy (γ -1)*M*c^2 in MeV
    W_max = maximum possible energy transfer to an electron in single collision in MeV
    k = bremsstrahlung photon energy in MeV
    z = charge number of incident particle
    Z = atomic number of absorber
    A = atomic mass of absorber in g mol^−1
    K = 0.307075 MeV mol^−1 cm^2
    I = mean excitation energy in eV
    δ(βγ) = density effect correction to ionization energy loss
    ℏω_p = plasma energy sqrt(ρ*〈Z/A〉) × 28.816 eV
    N_e = electron density in (units of r_e)^-3
    w_j = weight fraction of the jth element in a compound or mixture
    n_j = proportionality number of jth kind of atoms in a compound or mixture
    X_0 = radiation length of the absorber in g cm^−2
    E_c = critical energy for electrons in MeV
    E_μc = critical energy for muons in GeV
    E_s = scale energy sqrt(4π/α) * m_e * c^2 = 21.2052 MeV
    R_m = Molière radius in g cm^−2
*/
// Based on the above constants and values the Bethe formula can be constructed
// But also additionally for my ease of use I will simply separate those variables to be "separate" in a way, so I can easily identify which part does what
// And simply to know which part is responsible for what calculation
// Also I need to separate which variable is for incident particle and which one is for absorber... It gets confusing okay?


/* Absorber:
    ρ (rho) = density of the absorber in g cm^−3
    Z = atomic number of absorber
    A = atomic mass of absorber
*/

/* Incident Particle:
    M = incident particle mass in MeV/c^2
    E = incident particle energy γ*M*c^2 in MeV
    T = kinetic energy (γ -1)*M*c^2 in MeV
    W_max = maximum possible energy transfer to an electron in single collision in MeV
    z = charge number of incident particle
*/
use crate::modules;
use crate::modules::atomic_vars::PRECISION;
//TODO: Verify data and make sure calculation is correct
//TODO: Make corrections for the data formats that use either CSG or SI or Both

//Constants https://pdg.lbl.gov/2024/reviews/rpp2024-rev-phys-constants.pdf

// IDE reports that this function is being used over 99 times... I think I could cut down with that usage by reworking how periodic lookup works but either way it will be replaced by ANPM crate when it's done
// And then values returned by ANPM crate will be converted to Float so that rug::Float is used everywhere instead of f64
pub fn precise(val: &str) -> Float {
    let precision_bits: u32 = PRECISION.load(Ordering::Relaxed) as u32;
    let mut result = Float::new(precision_bits);
    result.assign(Float::parse(val).expect("Invalid float string"));
    result
}

//TODO: Do something akin of lazy loading or just not having to calculate or "something" to the values that don't change but still require to be used

//TODO: After moving to Float from f64 remove old unused const so that they won't waste space

fn pi_number_ret() -> Float {
    let precision_bits: u32 = PRECISION.load(Ordering::Relaxed) as u32;
    let mut pi_numb = Float::new(precision_bits);
    pi_numb.assign(Constant::Pi);
    pi_numb
}

fn light_speed_ret() -> Float {
    let light_spd = precise("299792458.0");
    light_spd
}

fn electron_mass_constant() -> Float {
    let electron_mass = precise("9.1093837e-31");
    electron_mass
}

fn planck_constant() -> Float {
    let planck = precise("6.62607015e-34");
    planck
}

fn avogadro_constant() -> Float {
    let avogadro = precise("6.02214076e23"); // mol^-1
    avogadro
}

pub fn stopping_power_intermediate_energies(name_of_incident_particle: &str, name_of_absorber: &str) {
    // https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf
    // Equation 34.5

    println!("Name of incident particle: {}", name_of_incident_particle);
    println!("Name of absorber {}", name_of_absorber);

    let element_exci_energy = mean_excitation_energy(name_of_absorber).unwrap_or_else(|| { // This is returned in eV, but it should be converted to MeV to not make unit mismatch
        eprintln!("Element {} not found or density is unavailable.", name_of_absorber);
        precise("0.0")
    });

    let energy: Float = m_e_cpowit();
    let mut de_dx_array: Vec<Float> = vec![];
    let mut velocity: Vec<Float> = vec![];
    println!("Mean excitation energy of absorber {} : {}", name_of_absorber.clone(), element_exci_energy.clone());

    // TODO: Loop for going through various values for BETA and GAMMA. For example Beta = (0.1*C/C) 0.1*C can be considered V

    for i in 1i32..1000i32 {
        // TODO: Verify if the loop values are correct and make sense physically.
        // TODO: Verify if the step of 0.001 makes sense or should it be smaller/larger
        // TODO: Verify if the gamma and beta calculations are correct and make sense physically.
        let i_t: Float = precise(&format!("{}", i)) * precise("0.001");
        let beta: Float = (i_t.clone() * light_speed_ret()) / light_speed_ret();

        let denominator: Float = precise("1.0") - ((i_t.clone() * light_speed_ret()) / light_speed_ret()).pow(2);
        let gamma: Float;
        let velocity_cutoff: Float = precise("0.999");
        if denominator > 0.0  && i_t.clone() <= velocity_cutoff {
            gamma = precise("1.0") / denominator.sqrt();
            velocity.push(gamma.clone());
        } else {
            gamma = precise("0.0");
            // panic!();
        }

        // Equation taken from https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf at 34.2.3
        // This formula is for Stopping power at intermediate energies
        let de_dx: Float = k_z_two_z_a_1_b_two(name_of_incident_particle, beta.clone(), name_of_absorber) *
                (0.5 * (twom_e_ctwo_btwo_dtwo_w(beta.clone(), gamma.clone(), energy.clone(), element_exci_energy.clone(), name_of_absorber).ln() - beta.clone().pow(2) -
                density_effect_correction(beta.clone(), gamma.clone(), plasma_energy(name_of_incident_particle), name_of_incident_particle))); // Missing for now value of δ(βγ) which is currently just
        de_dx_array.push(de_dx);
    }
    println!("<-dE/dx>: {:?}", de_dx_array);
    // println!("Velocity: {:?}", velocity);

    // println!("dE/dx:  J/m");
}

fn m_e_cpowit() -> Float {
    // https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf#table.caption.1
    // Table 34.1
    // m_e * c^2
    // Electron mass x c^2
    let result: Float = precise("0.51099895000"); // MeV
    // 0.51099895000(15) MeV
    // No idea what 15 in brackets means in the PDG document. The best thing I can think of is repeated decimal but not sure.
    // let result = electron_mass_constant() * light_speed_ret().pow(2);
    result
}

fn density_effect_correction(beta: Float, gamma: Float, plasma: Float, mean_exci_energy: &str, ) -> Float {
    // Important to calculate Density effect correction
    // NOTE: From what I understand this is basically density effect
    // https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf
    // Equation 34.6
    // δ(βγ)/2 → ln(ℏωp/I) + ln βγ − 1/2
    // No idea if I should use this calculation or just focus on Sternheimer's parametrization... OR HOW IT EVEN WORKS.
    // let mut logarithm_plasma_energy: Float = precise("1.0");
    // if let Some(mean_exci_energy_for_plasma) = mean_excitation_energy(mean_exci_energy) {
    //     logarithm_plasma_energy = (plasma / mean_exci_energy_for_plasma).ln();
    // }
    //
    // let logarithm_beta_gamma: Float = (beta.clone() * gamma.clone()).ln();
    // let res_beta_gamma_plasma: Float = logarithm_plasma_energy.clone() + logarithm_beta_gamma.clone() - precise("0.5");

    // 34.2.5 Density effect
    // PDG: https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf
    // NOTE: And this should be correction for density effect...
    // Use Sternheimer parameterization (Equation 34.7)

    // δ(βγ)/2 = { 2 (ln 10)x - C^-,                x >= x1
    //           { 2 (ln 10)x - C^- + a(x1 - x)^k , x0 <= x < x1
    //           { 0,                               x < x0 (nonconductors)
    //           { δ_0 10^2(x-x_0),                 x < x0 (conductors)

    // x = log_10 βγ = log_10 (p/Mc)

    let log_beta_gamma: Float = (beta.clone() * gamma.clone()).log10(); // x
    // x1 and x0 are maybe derived from table of muons dE/dx and Range. (example for Al: https://pdg.lbl.gov/2024/AtomicNuclearProperties/MUE/muE_aluminum_Al.pdf)
    let mut res_beta_gamma_sternheimer: Float;

    if let Some((x_0, x_1, c_minus, delta_0)) = periodic_lookup::look_up_element_x_c_delta(mean_exci_energy) {
        if log_beta_gamma.clone() >= x_1.clone() {
            println!("Duh. Normal Beta Gamma bigger than x_1 WAS USED");
            res_beta_gamma_sternheimer = precise("2.0") * precise("10").ln() * log_beta_gamma.clone() - c_minus.clone();
        } else if x_0.clone() <= log_beta_gamma.clone() && log_beta_gamma.clone() < x_1.clone() {
            if let Some((a_param, k_param)) = periodic_lookup::look_up_element_k_and_a(mean_exci_energy) {
                println!("K AND A WERE USED");
                res_beta_gamma_sternheimer = precise("2.0") * precise("10").ln() * log_beta_gamma.clone() - c_minus.clone() +
                a_param.clone() * (x_1.clone() - log_beta_gamma.clone()).pow(k_param.clone());
            } else {
                println!("Look up element x_c_delta failed...");
                res_beta_gamma_sternheimer = precise("2.0") * precise("10").ln() * log_beta_gamma.clone() - c_minus.clone();
            }
        } else {
            // TODO: Find a way to differentiate between conductors and nonconductors
            if delta_0.clone() == precise("0") {
                println!("delta 0 clone WAS USED");
                res_beta_gamma_sternheimer = precise("0");
            } else {
                println!("delta 0 clone ELSE WAS USED");
                res_beta_gamma_sternheimer = delta_0.clone() * precise("10").pow(precise("2") * log_beta_gamma.clone() - x_0.clone());
            }
        }
    } else {
        res_beta_gamma_sternheimer = precise("0");
    }

    // println!("Parametr sternh: {}", res_beta_gamma_sternheimer.clone()/precise("2"));
    // I know I need to use Sternheimer parameterization for this but for now it's better than nothing.
    let result = res_beta_gamma_sternheimer / precise("2");
    result
}

fn plasma_energy(name_of_absorber: &str) -> Float {
    // Important to calculate Density effect correction
    // https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf
    // ℏωp = √ρ*〈Z/A〉 × 28.816 eV
    // ρ in g cm^-3

    // 27.10.2025
    // In theory, I could just use ANPM (WHEN IT'S DONE FOR THE LOVE OF GOD I AM STILL NOT DONE WITH THAT CRATE AND I NEED IT.) and get values for ℏωp from there as for example
    // plasma energy of Si is just 31.05 eV (Taken from here https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/silicon_Si.html).
    // So I might just skip calculating it as I might be not sure if it's being calculated correctly (But even so the calculations seem close enough as I don't round up the result)

    // TODO: Verify if the calculation is correct as for now it might be simply wrong because of my lack of knowledge

    // TODO: Figure out a way to calculate for compounds
    // Based on the data found for example on https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/calcium_fluoride.html it seems like it does provide <Z/A> value

    let calc_const: Float = precise("28.816");
    if let Some((atom_density, atom_number, mass_number)) = periodic_lookup::look_up_element(name_of_absorber) {
        let mut result: Float = atom_density.clone() * (atom_number.clone() / mass_number.clone());
        result = result.sqrt() * calc_const;
        // println!("Plasma energy: {}", result.clone());
        result
    } else {
        eprintln!("Element {} not found or density is unavailable.", name_of_absorber);
        precise("0.0")
    }
    // let result: Float = precise("32.86");
    // result
}

fn shell_correction() -> Float{
    //TODO: Implement shell correction formula

    precise("1.0")
}

fn k_z_two_z_a_1_b_two(name_of_absorber: &str, beta: Float, name_of_incident_particle: &str) -> Float {
    //K = 4π * N * A* r^2_e * m_e * c^2     0.307075 MeV mol^−1 cm^2
    //z = -1 for electron
    //z = +1 for proton
    //Example for particle
    //z = +5 for Boron

    // From what I learned z in practice should accept custom inputs because of how charge number of incident particle works.
    // So more correct for Boron would be -3 and +3. Not +5 because it shouldn't be element atom_number

    if let Some((atom_density, atom_number, mass_number)) = periodic_lookup::look_up_element(name_of_incident_particle) {
        let z_inci = if name_of_incident_particle == "Ele" {
            precise("-1.0")
        } else if name_of_incident_particle == "Proto" {
            precise("1.0")
        } else if name_of_incident_particle != "Ele" && name_of_incident_particle != "Proto"  && name_of_incident_particle != "" {
            atom_number.clone()
        } else {
            // Default electron
            precise("1.0")
        };

        let k_z_two: Float = precise("0.307075") * z_inci.pow(2);
        let z_by_a: Float = atom_number.clone() / mass_number.clone();
        let one_by_beta: Float = precise("1.0") / beta.pow(2);
        let result: Float = k_z_two * z_by_a * one_by_beta;
        result
    } else {
        eprintln!("Element {} not found or density is unavailable.", name_of_absorber);
        precise("0.0")
    }
}

fn twom_e_ctwo_btwo_dtwo_w(beta: Float, gamma: Float, m_e_cpowit: Float, element_exci_energy: Float, name_of_incident_particle: &str) -> Float {
    // 2 * m_e * c^2 * β^2 * γ^2 * W_max
    // / I^2

    // Clone the values before using them in wmax
    let wmax_result: Float = wmax(beta.clone(), gamma.clone(), m_e_cpowit.clone(), name_of_incident_particle);
    let twom_e_ctwo: Float = precise("2.0") * m_e_cpowit.clone() * beta.clone().pow(2) * gamma.clone().pow(2) * wmax_result;

    // What the hell... Why did I call mean_excitation_energy when I already passed down element_exci_energy.
    let result: Float = twom_e_ctwo.clone() / element_exci_energy.clone();
    result
}

fn calculate_incident_particle_mass(name_of_incident_particle: &str) -> Float{
    // M    incident particle mass    MeV/c^2

    let uni_amu: Float = precise("931.4941024228");
    if (name_of_incident_particle == "Ele" || name_of_incident_particle == "Proto") {
        if (name_of_incident_particle == "Ele") {
            let result = precise("0.51099895000"); // MeV/c2
            result
        } else {
            let result = precise("938.27208816"); // MeV/c2
            result
        }
    } else {
        if let Some((mass_number)) = periodic_lookup::look_up_element_mass(name_of_incident_particle) {
            let result = mass_number * uni_amu;
            result
        } else {
            eprintln!("Element {} not found or density is unavailable. (Incident particle mass)", name_of_incident_particle);
            precise("0.0")
        }
    }
}

fn wmax(beta: Float, gamma: Float, m_e_cpowit: Float, name_of_incident_particle: &str) -> Float {
    // https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf
    // 34.2.2 Maximum energy transfer to an electron in a single collision

    // W_max = (2 * m_e * c^2 * β^2 * γ^2)
    // /(1 + 2 * γ * m_e / M + (m_e / M)^2)


    // 26.10.2025
    // Electron mass is required, but it simply confuses me because in the paper it uses both m_e * c^2 and m_e alone without c^2. In constants, it only states for m_e without c^2.
    // So it doesn't make any sense to me and how to differentiate when to use which one.
    // In this paper https://pdg.lbl.gov/2024/reviews/rpp2024-rev-phys-constants.pdf it states electron mass as m_e = 0.51099895000 MeV/c^2.
    // But in the review for passage of particles through matter it uses m_e * c^2 = 0.51099895000 MeV.
    // So which one is it supposed to be? Because in the formula for W_max it uses both m_e * c^2 and m_e alone without c^2.
    // In notation https://pdg.lbl.gov/2024/reviews/rpp2024-rev-passage-particles-matter.pdf#table.caption.1 it states that it's m_e * c^2 = 0.51099895000 MeV. and it's definition is electron mass x c^2.

    // Worst possible source but let's go with it for now. Based on what I found on wikipedia for electron mass (https://en.wikipedia.org/wiki/Electron_mass).
    // Both m_e and m_e*c^2 values are the same. Only difference between those two are units. One is with MeV/c^2 and other one is with just MeV

    //TODO: Not compatible with electrons. Need to do find the calculation for Wmax that doesn't divide electron mass by electron mass.
    if (name_of_incident_particle == "Ele") {
        let result = precise("1.0");
        result
    } else {
        let two_m_e_ctwo: Float = precise("2.0") * m_e_cpowit.clone() * beta.clone().pow(2) * gamma.clone().pow(2);
        let one_two_gamma_m_e: Float = precise("1.0") + precise("2.0") * gamma * m_e_cpowit.clone() //precise("0.51099895000"))
            / calculate_incident_particle_mass(name_of_incident_particle).clone()
            + (m_e_cpowit.clone() //(precise("0.51099895000")
            / calculate_incident_particle_mass(name_of_incident_particle).clone()).pow(2);
        //TODO: Verify if the incident particle mass is calculated correctly.
        let result = two_m_e_ctwo/one_two_gamma_m_e;
        result
    }
}