//! Traits

mod builder;
mod entity;
mod geometry;
mod geometry_map;
mod grid;
mod io;
mod topology;

pub use builder::Builder;
#[cfg(feature = "mpi")]
pub use builder::ParallelBuilder;
pub(crate) use builder::{GeometryBuilder, GridBuilder, TopologyBuilder};
pub use entity::Entity;
pub use geometry::{Geometry, Point};
pub use geometry_map::GeometryMap;
pub use grid::Grid;
#[cfg(feature = "mpi")]
pub use grid::{DistributableGrid, ParallelGrid};
#[cfg(feature = "serde")]
pub(crate) use io::ConvertToSerializable;
pub use io::{GmshExport, GmshImport};
#[cfg(feature = "serde")]
pub use io::{RONExport, RONImport};
#[cfg(feature = "mpi")]
#[cfg(feature = "serde")]
pub use io::{RONExportParallel, RONImportParallel};
pub use topology::Topology;
