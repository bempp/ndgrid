//! Grid geometry
mod geometry_map;
mod point;
pub(crate) mod single_element;
pub(crate) mod mixed;
pub use geometry_map::GeometryMap;
pub use point::{Point, PointIter};
pub use single_element::{
    SingleElementEntityGeometry, SingleElementEntityGeometryBorrowed, SingleElementGeometry,
    SingleElementGeometryBorrowed,
};
pub use mixed::{
    MixedEntityGeometry, MixedEntityGeometryBorrowed, MixedGeometry,
    MixedGeometryBorrowed,
};
