//! Traits

mod builder;
mod entity;
mod geometry;
mod geometry_map;
mod grid;
#[cfg(feature = "mpi")]
mod parallel;
mod topology;

pub use builder::{Builder, GeometryBuilder, GridBuilder, TopologyBuilder};
pub use entity::{Entity, EntityId};
pub use geometry::{Geometry, Point};
pub use geometry_map::GeometryMap;
pub use grid::Grid;
#[cfg(feature = "mpi")]
pub use parallel::{ParallelBuilder, ParallelGrid};
pub use topology::Topology;
