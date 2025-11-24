//! Data structures to work with n-dimensional grids.
//!
//! `ndgrid` builds upen [ndelement] to provide data structures for grids on either a single node or distributed via MPI.
//!
//! ## Creating a grid with `ndgrid`
//!
//! To explain the library we use the following example of a grid consisting of two triangles that together form the unit rectangle.
//! To that effect we introduce the following boundary vertices.
//! - Point 0: (0, 0)
//! - Point 1: (1, 0)
//! - Point 2: (0, 1)
//! - Point 3: (1, 1)
//!
//! To make matters more interesting we will define a grid of second order elements.
//! This means that each edge of the triangle also has a middle point. The corresponding points are
//! given as follows:
//!
//! - Point 4: (0.5, 0.5)
//! - Point 5: (0.0, 0.5)
//! - Point 6: (0.5, 0.0)
//! - Point 7: (0.5, 1.0)
//! - Point 8: (1.0, 0.5)
//!
//! The order of points for each element is the same as the one on [defelement.org](https://defelement.org).
//! For second order triangles we use
//! [this ordering](https://defelement.org/elements/examples/triangle-lagrange-equispaced-2.html).
//!
//! `ndgrid` does an important distinction between points and topological vertices. The boundary points
//! of the triangle are called vertices and determine the connectivity relationship of this triangle with other
//! triangles. Topologically, the two triangles are defined through the points 0, 1, 2 for the first triangle
//! and 1, 3, 2 for the second triangle. However, the full cell definition is
//! - Cell 1: 0, 1, 2, 4, 5, 6
//! - Cell 2: 1, 3, 2, 7, 4, 8
//!
//! in terms of the points.
//!
//! Let us generate the corresponding data structures.
//!
//! ```
//! use ndgrid::traits::Builder;
//! use ndelement::types::ReferenceCellType;
//! let mut builder = ndgrid::SingleElementGridBuilder::<f64>::new_with_capacity(
//!     2,
//!     9,
//!     2,
//!     (ReferenceCellType::Triangle, 2),
//! );
//! ```
//! The [SingleElementGridBuilder] is for grids that
//! only use a single element type. The first parameter is the geometric dimension. Here we choose 2,
//! meaning the grid lives in two-dimensionals pace. The next parameter is the number of points, 9 in this case,
//! and the third parameter is the number of cells, which is also 2 here.
//! The element type is [ReferenceCellType::Triangle](ndelement::types::ReferenceCellType::Triangle).
//! The order of the triangles is 2.
//!
//! We now add the definitions of the points and cells.
//! ```
//! # use ndgrid::traits::Builder;
//! # use ndelement::types::ReferenceCellType;
//! # let mut builder = ndgrid::SingleElementGridBuilder::<f64>::new_with_capacity(
//! #     2,
//! #     9,
//! #     2,
//! #     (ReferenceCellType::Triangle, 2),
//! # );
//! builder.add_point(0, &[0.0, 0.0]);
//! builder.add_point(1, &[1.0, 0.0]);
//! builder.add_point(2, &[0.0, 1.0]);
//! builder.add_point(3, &[1.0, 1.0]);
//! builder.add_point(4, &[0.5, 0.5]);
//! builder.add_point(5, &[0.0, 0.5]);
//! builder.add_point(6, &[0.5, 0.0]);
//! builder.add_point(7, &[0.5, 1.0]);
//! builder.add_point(8, &[1.0, 0.5]);
//!
//! builder.add_cell(0, &[0, 1, 2, 4, 5, 6]);
//! builder.add_cell(1, &[1, 3, 2, 7, 4, 8]);
//! ```
//! Finally, we generate the grid.
//! ```
//! # use ndgrid::traits::Builder;
//! # use ndelement::types::ReferenceCellType;
//! # let mut builder = ndgrid::SingleElementGridBuilder::<f64>::new_with_capacity(
//! #     2,
//! #     9,
//! #     2,
//! #     (ReferenceCellType::Triangle, 2),
//! # );
//! # builder.add_point(0, &[0.0, 0.0]);
//! # builder.add_point(1, &[1.0, 0.0]);
//! # builder.add_point(2, &[0.0, 1.0]);
//! # builder.add_point(3, &[1.0, 1.0]);
//! # builder.add_point(4, &[0.5, 0.5]);
//! # builder.add_point(5, &[0.0, 0.5]);
//! # builder.add_point(6, &[0.5, 0.0]);
//! # builder.add_point(7, &[0.5, 1.0]);
//! # builder.add_point(8, &[1.0, 0.5]);
//!
//! # builder.add_cell(0, &[0, 1, 2, 4, 5, 6]);
//! # builder.add_cell(1, &[1, 3, 2, 7, 4, 8]);
//! let grid = builder.create_grid();
//! ```
//! ## Querying the grid
//!
//! A grid is a hierarchy of entities. The highest dimension entities are the cells. Each cell consists of subentities,
//! which are faces, edges, etc. For each entity there are two types of information, the topology information and the
//! geometry information. The topology describes how entities are connected. The geometry describes how entities related
//! to their associated physical points. Each entity is associated with an `index`. Indices are unique within the class
//! of entities, that is there is a point with index 0 and a cell with index 0 but no two points with index 0. Points and
//! cells also have an associated `id`. `ids` are the indices provided by the user with the `add_point` or `add_cell`
//! methods in the grid builder. These ids will usually be different from the internal indices of entities.
//!
//! The following code extracts all topological vertices for each cell and prints the corresponding physical coordinates.
//! ```
//! # use ndgrid::traits::Builder;
//! # use ndelement::types::ReferenceCellType;
//! # let mut builder = ndgrid::SingleElementGridBuilder::<f64>::new_with_capacity(
//! #     2,
//! #     9,
//! #     2,
//! #     (ReferenceCellType::Triangle, 2),
//! # );
//! # builder.add_point(0, &[0.0, 0.0]);
//! # builder.add_point(1, &[1.0, 0.0]);
//! # builder.add_point(2, &[0.0, 1.0]);
//! # builder.add_point(3, &[1.0, 1.0]);
//! # builder.add_point(4, &[0.5, 0.5]);
//! # builder.add_point(5, &[0.0, 0.5]);
//! # builder.add_point(6, &[0.5, 0.0]);
//! # builder.add_point(7, &[0.5, 1.0]);
//! # builder.add_point(8, &[1.0, 0.5]);
//!
//! # builder.add_cell(0, &[0, 1, 2, 4, 5, 6]);
//! # builder.add_cell(1, &[1, 3, 2, 7, 4, 8]);
//! # let grid = builder.create_grid();
//! use ndgrid::traits::{Grid, Entity, Topology, Geometry, Point};
//! for cell in grid.entity_iter(ReferenceCellType::Triangle) {
//!     for vertex in cell.topology().sub_entity_iter(ReferenceCellType::Point) {
//!         let vertex = grid.entity(ReferenceCellType::Point, vertex).unwrap();
//!         let mut coords = [0.0; 2];
//!         vertex
//!             .geometry()
//!             .points()
//!             .next()
//!             .unwrap()
//!             .coords(&mut coords);
//!         println!(
//!             "Cell {} has vertex {} with coordinate [{}, {}]",
//!             cell.id().unwrap(),
//!             vertex.id().unwrap(),
//!             coords[0],
//!             coords[1]
//!         )
//!     }
//! }
//! ```
//! Let us dissect what is going on here. First, we iteratre through the cells of the grid.
//! For this we use the [Grid::entity_iter](crate::traits::Grid::entity_iter) function.
//! For each cell we then access the topology information via [Entity::topology](crate::traits::Entity::topology)
//! and iterate through the point subentities via [Topology::sub_entity_iter](crate::traits::Topology::sub_entity_iter).
//! This gives us the vertices of the triangles. The topology information only considers the points that define the topology.
//! So the middle points on each edge which are necessary for the order of the triangle, are not returned. Also, the iterator
//! returns integer indices of entities. To convert an entity index to an actual entity use the
//! [Grid::entity](crate::traits::Grid::entity) function. We now want to get the actual physical coordinate of a vertex.
//! Since the geometric dimension is 2 we instantiate an array `[f64; 2]` for this. We now call on the vertex the
//! [Entity::geometry](crate::traits::Entity::geometry) function to obtain its geometry information. We then
//! call [Geometry::points](crate::traits::Geometry::points) to get an iterator to all physical points
//! associated with the vertex. Since a vertex only has one associated physical point (namely the vertex itself) we just
//! call `next` once on this iterator to get the actual point. Finally, we call [Point::coords](crate::traits::Point::coords)
//! to get the values of the physical coordinate.
//!
//!
//!
//!
//!
//!

#![cfg_attr(feature = "strict", deny(warnings), deny(unused_crate_dependencies))]
#![warn(missing_docs)]

pub mod geometry;
pub mod grid;
mod io;
pub mod shapes;
pub mod topology;
pub mod traits;
pub mod types;

#[cfg(feature = "mpi")]
pub use grid::ParallelGridImpl;
pub use grid::{MixedGrid, MixedGridBuilder, SingleElementGrid, SingleElementGridBuilder};
pub use ndelement;

// Hack to avoid unused dependency errors if partitioner features are used without the mpi feature
#[cfg(not(feature = "mpi"))]
#[cfg(feature = "coupe")]
use coupe as _;
#[cfg(not(feature = "mpi"))]
#[cfg(feature = "scotch")]
use coupe as _;
