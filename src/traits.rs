//! Traits

mod builder;
mod entity;
mod geometry;
mod geometry_map;
mod grid;
mod io;
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
#[cfg(feature = "serde")]
pub use io::{RONExportParallel, RONImportParallel};
pub use parallel::{ParallelBuilder, ParallelGrid};
pub use topology::Topology;
