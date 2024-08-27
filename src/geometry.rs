//! Grid geometry
mod geometry_map;
mod point;
pub(crate) mod single_element;
pub use geometry_map::GeometryMap;
pub use point::{Point, PointIter};
pub use single_element::{SingleElementEntityGeometry, SingleElementGeometry};
