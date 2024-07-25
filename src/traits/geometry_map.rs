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

    /// Write the physical points for the entity with index `entity_index` into `points`
    ///
    /// `points` should have shape [geometry_dimension, npts] and use column-major ordering
    fn points(&self, entity_index: usize, points: &mut [Self::T]);

    /// Write the jacobians at the physical points for the entity with index `entity_index` into `jacobians`
    ///
    /// `jacobians` should have shape [geometry_dimension, entity_topology_dimension, npts] and use column-major ordering
    fn jacobians(&self, entity_index: usize, jacobians: &mut [Self::T]);

    /// Write the jacobians, their determinants, and the normals at the physical points for the entity with
    /// index `entity_index` into `jacobians`, `jdets` and `normals`
    ///
    /// `jacobians` should have shape [geometry_dimension, entity_topology_dimension, npts] and use column-major ordering;
    /// `jdets` should have shape \[npts\];
    /// `normals` should have shape [geometry_dimension, npts] and use column-major ordering
    fn jacobians_dets_normals(
        &self,
        entity_index: usize,
        jacobians: &mut [Self::T],
        jdets: &mut [Self::T],
        normals: &mut [Self::T],
    );
}
