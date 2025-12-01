//! Data structures to work with n-dimensional grids.
//!
//! `ndgrid` builds upen [ndelement] to provide data structures for grids on either a single node or distributed via MPI.
//!
//! ## Creating a grid with `ndgrid`
//!
//! To demonstrate the library, we use an example grid consisting of two triangles that together form the unit square.
//! We introduce the following points. As we will make a grid of second oder elements, we include points at the midpoint
//! of each edge as well as at the vertices of the square
//! - Point 0: (0, 0)
//! - Point 1: (1, 0)
//! - Point 2: (0, 1)
//! - Point 3: (1, 1)
//! - Point 4: (0.5, 0.5)
//! - Point 5: (0.0, 0.5)
//! - Point 6: (0.5, 0.0)
//! - Point 7: (0.5, 1.0)
//! - Point 8: (1.0, 0.5)
//!
//! The order of points for each cell follows the point ordering of the
//! [Lagrange element on DefElement](https://defelement.org/elements/lagrange.html).
//! For second order triangles we use
//! [this ordering](https://defelement.org/elements/examples/triangle-lagrange-equispaced-2.html).
//! Following this ordering, the two cells of our example grid are:
//!
//! - Cell 1: 0, 1, 2, 4, 5, 6
//! - Cell 2: 1, 3, 2, 7, 4, 8
//!
//! In order to create our grid using ndgrid, we first create a grid builder.
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
//!
//! The [SingleElementGridBuilder] is for grids that only use a single element type. The parameters passed when
//! initialising the build are:
//!
//! - The geometric dimension: our example grid lives in two-dimensional space, so we use 2.
//! - The number of points: 9 for our example grid.
//! - The number of cells: 2 for our example grid.
//! - The cell type and element degree: for our example, this is ([ReferenceCellType::Triangle](ndelement::types::ReferenceCellType::Triangle), 2)
//!   as our geometry cells are triangles and we use quadratic geometry for each triangle.
//!
//! If we did not know the number of points and cells that we will include in out grid when creating ths builder,
//! we could instead use the function [SingleElementGridBuilder::new] when initialising the grid.
//!
//! Now that we have created a grid builder, we can add the points and cells:
//!
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
//! builder.add_cell(1, &[0, 1, 2, 4, 5, 6]);
//! builder.add_cell(2, &[1, 3, 2, 7, 4, 8]);
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
//! #
//! # builder.add_cell(1, &[0, 1, 2, 4, 5, 6]);
//! # builder.add_cell(2, &[1, 3, 2, 7, 4, 8]);
//! let grid = builder.create_grid();
//! ```
//!
//! ## Querying the grid
//!
//! A grid is a hierarchy of entities. We follow the standard name conventions for entities of a given topological dimension:
//! 0-, 1-, 2- and 3-dimensional entities and called vertices, edges, faces and volumes (respectively).
//! The highest dimensional entities are called cells. If $d$ the (topological) dimension of the cells,
//! then $d-1$-, $d-2$- and $d-3$-dimensional entities are called facets, ridges and peaks (respectively).
//!
//! For each entity there are two types of information: the topology and the geometry.
//! The topology describes how entities are connected. The geometry describes how entities are positioned in physical space.
//! As the topology is only concerned with the connectivity between entities, it only includes the cell's points that are at
//! the vertices of a cell (eg for the triangle cells in our grid, the topology onle includes the first three points for each cell).
//! In the geometry, all the points that define the cell are stored.
//!
//! Each entity has an associated `index`. Indices are unique within entities of a given type:
//! there is a vertex with index 0 and a cell with index 0 but there cannot be two vertices with index 0. Points and
//! cells may also have an associated `id`: these are the values provided by the user when using the `add_point` or `add_cell`
//! methods in the grid builder. These ids are not guaranteed to be equal to the indices of the entities.
//!
//! The following code extracts the vertices of each cell and prints their corresponding physical coordinates.
//! ```
//! # use ndgrid::traits::Builder;
//! use ndgrid::traits::{Grid, Entity, Topology, Geometry, Point};
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
//! #
//! # builder.add_cell(1, &[0, 1, 2, 4, 5, 6]);
//! # builder.add_cell(2, &[1, 3, 2, 7, 4, 8]);
//! # let grid = builder.create_grid();
//!
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
//!
//! This snippets starts by using [Grid::entity_iter](crate::traits::Grid::entity_iter) to iterate through each
//! cell (ie each entity that is a triangle).
//! For each cell, we then access the topology information via [Entity::topology](crate::traits::Entity::topology)
//! and iterate through the vertices (ie the subentities that are points) using [Topology::sub_entity_iter](crate::traits::Topology::sub_entity_iter).
//! This iterators gives us the index of each vertex: to convert an entity index to an entity, we use [Grid::entity](crate::traits::Grid::entity).
//! We now want to get the actual physical coordinate of a vertex.
//! Since the geometric dimension is 2 we instantiate an array `[f64; 2]` for this. We use
//! [Entity::geometry](crate::traits::Entity::geometry) to obtain the geometry for the vertex, then use
//! [Geometry::points](crate::traits::Geometry::points) to get an iterator over the physical points
//! associated with the vertex. Since a vertex has only one associated physical point, we
//! call `next` once on this iterator to get the point. Finally, we call [Point::coords](crate::traits::Point::coords)
//! to get the values of the physical coordinate.

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
use scotch as _;
