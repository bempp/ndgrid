//! Traits

mod builder;
mod entity;
mod geometry;
mod geometry_map;
mod grid;
mod io;
mod topology;

pub use builder::{Builder, GeometryBuilder, GridBuilder, ParallelBuilder, TopologyBuilder};
pub use entity::Entity;
pub use geometry::{Geometry, Point};
pub use geometry_map::GeometryMap;
pub use grid::{Grid, ParallelGrid};
#[cfg(feature = "serde")]
pub(crate) use io::ConvertToSerializable;
pub use io::{GmshExport, GmshImport};
#[cfg(feature = "serde")]
pub use io::{RONExport, RONExportParallel, RONImport, RONImportParallel};
pub use topology::Topology;
