//! Types

mod ownership;
pub use ownership::Ownership;

use num::Float;
use rlst::{LinAlg, RlstScalar};

/// A real scalar
pub trait RealScalar: Float + LinAlg + RlstScalar<Real = Self> {}
