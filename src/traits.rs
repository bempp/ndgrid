//! Traits

mod entity;
mod geometry;
mod grid;
mod point;
mod topology;

pub use entity::{Entity, EntityGeometry, EntityId, EntityTopology};
pub use geometry::Geometry;
pub use grid::Grid;
pub use point::Point;
pub use topology::Topology;
