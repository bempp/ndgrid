//! Traits

mod builder;
mod entity;
mod geometry;
mod geometry_map;
mod grid;
#[cfg(features = "mpi")]
mod parallel;
mod topology;

pub use builder::Builder;
pub use entity::{Entity, EntityId};
pub use geometry::{Geometry, Point};
pub use geometry_map::GeometryMap;
pub use grid::Grid;
#[cfg(features = "mpi")]
pub use parallel::{ParallelBuilder, ParallelGrid};
pub use topology::Topology;
