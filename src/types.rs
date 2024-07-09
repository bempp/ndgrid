//! Types

mod ownership;
pub use ownership::Ownership;
mod connectivity;
pub use connectivity::CellLocalIndexPair;

use num::Float;
use rlst::{Array, BaseArray, LinAlg, RlstScalar, VectorContainer};

/// A real scalar
pub trait RealScalar: Float + LinAlg + RlstScalar<Real = Self> {}

/// A 2-dimensional array
pub type Array2D<T> = Array<T, BaseArray<T, VectorContainer<T>, 2>, 2>;
