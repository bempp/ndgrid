//! Types
use mpi::traits::Equivalence;
use num::Float;
use rlst::{Array, BaseArray, LinAlg, RlstScalar, SliceContainer, VectorContainer};

/// A real scalar
pub trait RealScalar: Float + LinAlg + RlstScalar<Real = Self> {}

impl RealScalar for f32 {}
impl RealScalar for f64 {}

/// An N-dimensional array
pub type ArrayND<const N: usize, T> = Array<T, BaseArray<T, VectorContainer<T>, N>, N>;
/// A 2-dimensional array
pub type Array2D<T> = ArrayND<2, T>;

/// An N-dimensional array with borrowed data
pub type ArrayNDBorrowed<'a, const N: usize, T> =
    Array<T, BaseArray<T, SliceContainer<'a, T>, N>, N>;
/// A 2-dimensional array with borrowed data
pub type Array2DBorrowed<'a, T> = ArrayNDBorrowed<'a, 2, T>;

/// A (cell, local index) pair
///
/// The local index is the index of a subentity (eg vertex, edge) within the cell as it is numbered in the reference cell
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, PartialEq, Eq, Clone, Copy, Hash)]
#[repr(u8)]
pub enum Ownership {
    /// Undefined ownership. This should only be used as a placeholder value
    Undefined,
    /// Owned by the current process
    Owned,
    /// Not owned on the current process. The two values are the process that owns it
    /// and its local index on that process
    Ghost(usize, usize),
}

unsafe impl Equivalence for Ownership {
    type Out = <u8 as Equivalence>::Out;
    fn equivalent_datatype() -> <u8 as Equivalence>::Out {
        <u8 as Equivalence>::equivalent_datatype()
    }
}

/// Graph partitioner
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Debug, PartialEq, Eq, Clone, Hash)]
#[repr(u8)]
pub enum GraphPartitioner {
    /// No partioner: cells are split into approximately equal sized lists by index
    None,
    /// A manual partition
    Manual(Vec<usize>),
    #[cfg(feature = "coupe")]
    /// Use Coupe's KMeans implementation to generate a partition
    Coupe,
    #[cfg(feature = "scotch")]
    /// Use Scotch to generate a partition
    Scotch,
}
