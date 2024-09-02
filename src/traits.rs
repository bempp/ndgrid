//! Traits

mod builder;
mod entity;
mod geometry;
mod geometry_map;
mod grid;
mod io;
#[cfg(feature = "mpi")]
mod parallel;
mod topology;

pub use builder::{Builder, GeometryBuilder, GridBuilder, TopologyBuilder};
pub use entity::Entity;
pub use geometry::{Geometry, Point};
pub use geometry_map::GeometryMap;
pub use grid::Grid;
#[cfg(feature = "serde")]
pub(crate) use io::ConvertToSerializable;
pub use io::GmshExport;
#[cfg(feature = "serde")]
pub use io::{RONExport, RONImport};
#[cfg(all(feature = "serde", feature = "mpi"))]
pub use io::{RONExportParallel, RONImportParallel};
#[cfg(feature = "mpi")]
pub use parallel::{ParallelBuilder, ParallelGrid};
pub use topology::Topology;
