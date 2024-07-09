//! Traits

mod entity;
mod geometry;
mod grid;
mod topology;

pub use entity::{Entity, EntityId};
pub use geometry::{Geometry, Point};
pub use grid::Grid;
pub use topology::Topology;
