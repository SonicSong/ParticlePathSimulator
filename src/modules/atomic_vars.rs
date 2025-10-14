use std::sync::atomic::{AtomicI32, Ordering};

pub static PRECISION: AtomicI32 = AtomicI32::new(256);