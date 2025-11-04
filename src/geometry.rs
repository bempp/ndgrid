//! Grid geometry
mod geometry_map;
pub(crate) mod mixed;
mod point;
pub(crate) mod single_element;
pub use geometry_map::GeometryMap;
pub use mixed::{MixedEntityGeometry, MixedGeometry};
pub use point::{Point, PointIter};
pub use single_element::{SingleElementEntityGeometry, SingleElementGeometry};
