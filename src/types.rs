//! Types
#[cfg(feature = "mpi")]
use mpi::traits::Equivalence;
use rlst::{
    RlstScalar,
    dense::linalg::lapack::interface::{getrf::Getrf, getri::Getri},
};

/// A scalar.
pub trait Scalar: RlstScalar + Getrf + Getri {}

impl<T: RlstScalar + Getrf + Getri> Scalar for T {}

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

#[cfg(feature = "mpi")]
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
