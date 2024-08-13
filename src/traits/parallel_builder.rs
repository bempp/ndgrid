//! Builder for MPI parallel grids
use super::Builder;

#[cfg(feature = "mpi")]
pub trait ParallelBuilder: Builder {
    //! Trait to build a MPI parellelised mesh from a grid builder.
}
