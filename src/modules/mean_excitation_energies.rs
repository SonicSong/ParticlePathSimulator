use mendeleev::Element;

pub trait CustomProperties {
    fn custom_mean_excitation_energy(&self) -> Option<f64>;
}

impl CustomProperties for Element {
    fn custom_mean_excitation_energy(&self) -> Option<f64> {
        // Return None for elements without data
        // Skipping gasses and liquids
        match self {
            Element::Li => Some(40.0),
            Element::Be => Some(63.7),
            Element::B => Some(76.0),
            Element::Na => Some(149.0),
            Element::Mg => Some(156.0),
            Element::Al => Some(166.0),
            Element::Si => Some(173.0),
            Element::P => Some(173.0),
            Element::S => Some(180.0),
            Element::K => Some(190.0),
            Element::Ca => Some(191.0),
            Element::Sc => Some(216.0),
            Element::Ti => Some(233.0),
            Element::V => Some(245.0),
            Element::Cr => Some(257.0),
            Element::Mn => Some(272.0),
            Element::Fe => Some(286.0),
            Element::Co => Some(297.0),
            Element::Ni => Some(311.0),
            Element::Cu => Some(322.0),
            Element::Zn => Some(330.0),
            Element::Ga => Some(334.0),
            Element::Ge => Some(350.0),
            Element::As => Some(347.0),
            Element::Se => Some(348.0),
            Element::Rb => Some(363.0),
            Element::Sr => Some(366.0),
            Element::Y => Some(379.0),
            Element::Zr => Some(393.0),
            Element::Nb => Some(417.0),
            Element::Mo => Some(424.0),
            Element::Tc => Some(428.0),
            Element::Ru => Some(441.0),
            Element::Rh => Some(449.0),
            Element::Pd => Some(470.0),
            Element::Ag => Some(470.0),
            Element::Cd => Some(469.0),
            Element::In => Some(488.0),
            Element::Sn => Some(488.0),
            Element::Sb => Some(487.0),
            Element::Te => Some(485.0),
            Element::I => Some(491.0),
            Element::Cs => Some(488.0),
            Element::Ba => Some(491.0),
            Element::La => Some(501.0),
            Element::Hf => Some(705.0),
            Element::Ta => Some(718.0),
            Element::W => Some(727.0),
            Element::Re => Some(736.0),
            _ => None,
        }
    }
}