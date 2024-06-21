//! Types

mod ownership;
pub use ownership::Ownership;

use num::Float;
use rlst::{LinAlg, RlstScalar};

pub trait RealScalar: Float + LinAlg + RlstScalar<Real = Self> {}
