//! Types
use num::Float;
use rlst::{Array, BaseArray, LinAlg, RlstScalar, VectorContainer};

/// A real scalar
pub trait RealScalar: Float + LinAlg + RlstScalar<Real = Self> {}

impl RealScalar for f32 {}
impl RealScalar for f64 {}

/// An N-dimensional array
pub type ArrayND<const N: usize, T> = Array<T, BaseArray<T, VectorContainer<T>, N>, N>;
/// A 2-dimensional array
pub type Array2D<T> = ArrayND<2, T>;

/// A (cell, local index) pair
///
/// The local index is the index of a subentity (eg vertex, edge) within the cell as it is numbered in the reference cell
#[derive(Debug, Clone)]
pub struct CellLocalIndexPair<IndexType: std::fmt::Debug + Eq + Copy> {
    /// The cell's index
    pub cell: IndexType,
    /// The local index of the subentity
    pub local_index: usize,
}

impl<IndexType: std::fmt::Debug + Eq + Copy> CellLocalIndexPair<IndexType> {
    /// Create a (cell, local index) pair
    pub fn new(cell: IndexType, local_index: usize) -> Self {
        Self { cell, local_index }
    }
}

/// Ownership
#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
pub enum Ownership {
    /// Owned by the current process
    Owned,
    /// Not owned on the current process. The two values are the process that owns it
    /// and its local index on that process
    Ghost(usize, usize),
}
