//! Map from reference to physical space.

use crate::types::Scalar;

/// A geometry map allows the computation of maps from reference to physical space and their derivatives.
///
/// A geometry map is typically initialised with a number of points on a reference entity.
/// We can then for each physical entity compute
/// - The associated physical points as maps from the reference points.
/// - The jacobian of the map at the physical points.
/// - The jacobians, transformation determinants and the normals of the physical entity.
pub trait GeometryMap {
    /// Scalar type
    type T: Scalar;

    /// The topoloical dimension of the entity being mapped.
    ///
    /// The topological dimension is e.g. two for a triangle, independent
    /// of whether it is embedded in two or three dimensional space.
    fn entity_topology_dimension(&self) -> usize;

    /// The geometric dimension of the physical space.
    fn geometry_dimension(&self) -> usize;

    /// The number of reference points that this map uses.
    fn point_count(&self) -> usize;

    /// Write the physical points for the entity with index `entity_index` into `points`
    ///
    /// `points` should have shape [geometry_dimension, npts] and use column-major ordering.
    fn physical_points(&self, entity_index: usize, points: &mut [Self::T]);

    /// Write the jacobians at the physical points for the entity with index `entity_index` into `jacobians`
    ///
    /// `jacobians` should have shape [geometry_dimension, entity_topology_dimension, npts] and use column-major ordering
    fn jacobians(&self, entity_index: usize, jacobians: &mut [Self::T]);

    /// Write the jacobians, their determinants, and the normals at the physical points for the entity with
    /// index `entity_index` into `jacobians`, `jdets` and `normals`.
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
