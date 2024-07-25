//! Map from reference to physical space.

use crate::types::RealScalar;

pub trait GeometryMap {
    //! Reference to physical geometry  map

    /// Scalar type
    type T: RealScalar;

    /// The topoloical dimension of the entity being mapped
    fn entity_topology_dimension(&self) -> usize;

    /// The geometric dimension of the physical space
    fn geometry_dimension(&self) -> usize;

    /// The number of reference points that this map uses
    fn point_count(&self) -> usize;

    /// Write the physical points for the entity with index `entity_index` into `value`
    ///
    /// `value` should have shape [geometry_dimension, npts] and use column-major ordering
    fn points(&self, entity_index: usize, value: &mut [Self::T]);

    /// Write the jacobians at the physical points for the entity with index `entity_index` into `value`
    ///
    /// `value` should have shape [geometry_dimension, entity_topology_dimension, npts] and use column-major ordering
    fn jacobians(&self, entity_index: usize, value: &mut [Self::T]);

    /// Write the normals at the physical points for the entity with index `entity_index` into `value`
    ///
    /// `value` should have shape [geometry_dimension, npts] and use column-major ordering
    fn normals(&self, entity_index: usize, value: &mut [Self::T]);
}
