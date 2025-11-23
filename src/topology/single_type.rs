//! Topology for grids where entities of each tdim are the same type

#[cfg(feature = "serde")]
use crate::traits::ConvertToSerializable;
use crate::traits::Topology;
use itertools::izip;
use ndelement::{orientation::compute_orientation, reference_cell, types::ReferenceCellType};
use rlst::{DynArray, rlst_dynamic_array};
use std::collections::HashMap;
use std::iter::Copied;

/// Topology of a single element grid
#[derive(Debug)]
pub struct SingleTypeTopology {
    dim: usize,
    pub(crate) ids: Vec<Option<Vec<usize>>>,
    pub(crate) ids_to_indices: Vec<HashMap<usize, usize>>,
    entity_types: Vec<ReferenceCellType>,
    entity_counts: Vec<usize>,
    pub(crate) downward_connectivity: Vec<Vec<DynArray<usize, 2>>>,
    pub(crate) upward_connectivity: Vec<Vec<Vec<Vec<usize>>>>,
    pub(crate) orientation: Vec<Vec<i32>>,
}

#[cfg(feature = "serde")]
#[derive(serde::Serialize, Debug, serde::Deserialize)]
/// Serde serializable topology
pub struct SerializableTopology {
    dim: usize,
    ids: Vec<Option<Vec<usize>>>,
    ids_to_indices: Vec<HashMap<usize, usize>>,
    entity_types: Vec<ReferenceCellType>,
    entity_counts: Vec<usize>,
    downward_connectivity: Vec<Vec<(Vec<usize>, [usize; 2])>>,
    upward_connectivity: Vec<Vec<Vec<Vec<usize>>>>,
    orientation: Vec<Vec<i32>>,
}

#[cfg(feature = "serde")]
impl ConvertToSerializable for SingleTypeTopology {
    type SerializableType = SerializableTopology;
    fn to_serializable(&self) -> SerializableTopology {
        SerializableTopology {
            dim: self.dim,
            ids: self.ids.clone(),
            ids_to_indices: self.ids_to_indices.clone(),
            entity_types: self.entity_types.clone(),
            entity_counts: self.entity_counts.clone(),
            downward_connectivity: self
                .downward_connectivity
                .iter()
                .map(|a| {
                    a.iter()
                        .map(|b| (b.data().unwrap().to_vec(), b.shape()))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>(),
            upward_connectivity: self.upward_connectivity.clone(),
            orientation: self.orientation.clone(),
        }
    }
    fn from_serializable(s: SerializableTopology) -> Self {
        Self {
            dim: s.dim,
            ids: s.ids,
            ids_to_indices: s.ids_to_indices,
            entity_types: s.entity_types,
            entity_counts: s.entity_counts,
            downward_connectivity: s
                .downward_connectivity
                .iter()
                .map(|a| {
                    a.iter()
                        .map(|(data, shape)| {
                            let mut c = DynArray::<usize, 2>::from_shape(*shape);
                            c.data_mut().unwrap().copy_from_slice(data);
                            c
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>(),
            upward_connectivity: s.upward_connectivity,
            orientation: s.orientation,
        }
    }
}

impl SingleTypeTopology {
    /// Create a topology
    pub fn new(
        cells: &[usize],
        cell_type: ReferenceCellType,
        vertex_ids: Option<Vec<usize>>,
        cell_ids: Option<Vec<usize>>,
    ) -> Self {
        let size = reference_cell::entity_counts(cell_type)[0];
        let ncells = cells.len() / size;
        let dim = reference_cell::dim(cell_type);

        // The reference connectivity is a 4-dimensional array.
        // The first dimension is the topological dimension of the entities for which we want the connectivity.
        // The second dimension is the index of the entity for which we want to enquire connectivity.
        // The third dimension is the dimension of the sub-entities for which we want to enquire connectivity.
        // Hence, for a triangle ref_connectivity[0][1][1] gives us the indices of the two lines that point one
        // is connected two.
        let ref_conn = reference_cell::connectivity(cell_type);

        // This gives us all the subentity-types of the cell type. etypes[0] is a list of
        // all the entity types of dimension 0 (i.e. vertices), etypes[1] is a list of all
        // the entity types of dimension 1, etc.
        let etypes = reference_cell::entity_types(cell_type);

        // The following iterates over all the entity types of each dimension and makes
        // sure that within a dimension all entity types are identical.
        // Cells where faces are mixture of triangles and quads not supported
        for t in &etypes {
            for i in t {
                if *i != t[0] {
                    panic!("Unsupported cell type for SingleTypeTopology: {cell_type:?}");
                }
            }
        }

        // This simply collects all the different entity types that are possible within an element.
        let entity_types = etypes
            .iter()
            .filter(|i| !i.is_empty())
            .map(|i| i[0])
            .collect::<Vec<_>>();

        // List of entities by dimension
        // Entity are stored by dimension. We are not storing vertices. Hence, dim - 1.
        // entities[0] are a HashMap of all zero dimensional entities to their respective entity index.
        // The key of this HashMap is the ordered set of associated vertex indices for this entity.
        let mut entities = Vec::<HashMap<Vec<usize>, usize>>::new();
        if dim > 0 {
            for _ in 0..dim - 1 {
                entities.push(HashMap::new());
            }
        }
        // We setup entity counters for each dimension.

        let mut entity_counter = vec![0; entities.len()];

        // Now iterate through all the cells
        for cell_index in 0..ncells {
            // Get all the vertices of the cell.
            let cell = &cells[cell_index * size..(cell_index + 1) * size];
            // Iterate over topological dimension.
            // The variable e_i is a mutable reference to the list of entities of dimension i.
            // The variable rc_i is a reference to the connectivity of the cell type for entities of dimension i.
            // Hence, rc_i[0] is the connectivity information about the entity with reference index 0 and
            // dimension i.
            for (entity_dim, (e_i, rc_i)) in izip!(
                entities.iter_mut(),
                // Skip the vertices in the connectivity
                ref_conn.iter().take(dim).skip(1),
            )
            .enumerate()
            {
                // We iterate through the concrete entities of dimension i and the corresponding
                // entity types. c_ij is a reference to the connectivity information of the jth entity
                // with dimension i.
                for c_ij in rc_i {
                    // c_ij[0] is the list of reference cell vertex indices that are associated with the jth entity of dimension i.
                    // cell[*i] below maps the local reference cell vertex index to the actual id of the vertex.
                    // Hence, the following command gives us all the vertex indicies of the entity.
                    let mut entity = c_ij[0].iter().map(|i| cell[*i]).collect::<Vec<_>>();
                    // We sort entities so that entities with same vertices but different order of vertices
                    // are treated as the same entity.
                    entity.sort();
                    // An entity only gets added if it has not been added before.
                    e_i.entry(entity).or_insert_with(|| {
                        let old = entity_counter[entity_dim];
                        entity_counter[entity_dim] += 1;
                        old
                    });
                }
            }
        }
        // Number of entities by dimension
        let mut entity_counts = vec![0; dim + 1];
        // Vertices are indexed consecutively starting from 0.
        // So to get the number of vertices take the highest index
        // and add 1.
        entity_counts[0] = cells.iter().max().unwrap() + 1;
        entity_counts[dim] = ncells;
        for d in 1..dim {
            // In entities the dimension 0 is skipped.
            // So entity_counts[d] is the number of entities in
            // entities[d - 1].
            entity_counts[d] = entities[d - 1].len();
        }

        // Downward connectivity: The entities of dimension dim1 that are subentities of
        // entities  of dimension dim0 (with dim0>=dim1) (eg edges of a triangle, vertices
        // of a tetrahedron, etc)
        // downward_connectivity[dim0][dim1][[.., dim0_entity_index]] = [dim1_entity_index]
        // Hence, downward_connectivity[2][1][0, 5] for a triangle grid gives the index of the
        // first edge of the 6th triangle.
        // The following loop only creates space for the downward connectivity. It does not fill
        // the actual entries.
        let mut downward_connectivity = entity_counts
            .iter()
            .enumerate()
            // i is the topological dimension of the entities for which we want the connectivity.
            // j is the reference cell index of the entity of dimension j.
            .map(|(i, j)| {
                // ref_conn[i][0] gives us the first entity of dimension i. We don't need to look
                // at other entities of dimension j since the information about number of
                // subentities is the same.
                ref_conn[i][0]
                    // We iterate through all the dimensions it can be connected to.
                    .iter()
                    // We take i + 1 dimensions of the connection.
                    // If i = 2 (triangle) then we take vertices, edges, and triangles.
                    // That's why take(i + 1).
                    .take(i + 1)
                    // We now iterate through the connected dimensions, e.g. in the case
                    // of i = 2, the dimensions 0, 1, 2. These are in r.
                    // r.len() is the number of these entities connected to our i-dim entity.
                    // and j is the actual number of i-dimensional entities.
                    .map(|r| rlst_dynamic_array!(usize, [r.len(), *j]))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        // Upward connectivity: The entities of dimension dim1 that are superentities of
        // entities of dimension dim0 (with dim0<dim1) (eg triangles connected to an edge,
        // tetrahedra connected to a vertex, etc)
        // upward_connectivity[dim0][dim1 - dim0 - 1][dim0_entity_index][..] = [dim1_entity_index]
        let mut upward_connectivity = entity_counts
            .iter()
            .take(dim)
            .enumerate()
            .map(|(i, j)| vec![vec![vec![]; *j]; dim - i])
            .collect::<Vec<Vec<Vec<Vec<usize>>>>>();

        // downward_connectivity[d][d][i] = [i] (ie each entity is a sub-entity of itself)
        for (d, dc) in downward_connectivity.iter_mut().enumerate() {
            for (i, mut j) in dc[d].col_iter_mut().enumerate() {
                // Each entity has exactly one sub-entity of the same dimension, namely itself.
                // Hence, the matrix is a 1x1 matrix that refers back to itself.
                j[[0]] = i;
            }
        }

        // downward_connectivity[dim][0] = vertices of each cell
        // This is just the vertices that each cell is made of.
        for (i, mut dc_d0i) in downward_connectivity[dim][0].col_iter_mut().enumerate() {
            for (dc_d0ij, c_j) in izip!(dc_d0i.iter_mut(), &cells[i * size..(i + 1) * size]) {
                *dc_d0ij = *c_j;
                // We can also fill out the upward connectivity if dim > 0.
                if dim > 0 {
                    upward_connectivity[0][dim - 1][*c_j].push(i);
                }
            }
        }

        // downward_connectivity[i][0] = vertices of entity
        for (i, (es_i, dc_i)) in izip!(
            entities.iter(),
            downward_connectivity.iter_mut().take(dim).skip(1)
        )
        .enumerate()
        {
            for es_ij in es_i.iter() {
                let entity_vertices = es_ij.0;
                let entity_index = es_ij.1;
                let mut dc_i0j = dc_i[0].r_mut().slice::<1>(1, *entity_index);
                for (vertex, dc_i0jk) in izip!(entity_vertices.iter(), dc_i0j.iter_mut()) {
                    *dc_i0jk = *vertex;
                    if !upward_connectivity[0][i][*vertex].contains(entity_index) {
                        upward_connectivity[0][i][*vertex].push(*entity_index);
                    }
                }
            }
        }

        // We now have to fill the upward and downward connectivity for the other cases. This can only happen if dim > 0.
        if dim > 0 {
            // We collect the cell entities in an array. We do not need the connectivity of the vertices as that
            // has already been dealt with. Hence, the skip(1). The index i specifies the subentity of the cell.
            // and i.len() is the number of subentities of the given dimension.
            // Hence, cell_entities[0] is a vector of length of all subentities of dimension 1 in the reference cell.
            let mut cell_entities = ref_conn
                .iter()
                .skip(1)
                .map(|i| vec![0; i.len()])
                .collect::<Vec<_>>();
            // Now we iterate through each individual cell.
            for cell_index in 0..ncells {
                // Collect indices of each subentity of the cell.
                // The subentities of dimension the cell itself are just the cell indicies. Remember that a dimension index of
                // dim - 1 denotes actually the dim dimensional subentity as we skipped vertices.
                cell_entities[dim - 1][0] = cell_index;
                // We now get all the vertices of the current cell.
                let cell = &cells[cell_index * size..(cell_index + 1) * size];
                // We now iterate through the connectivity dimensions. e_i is all entities of dimension i + 1.
                // ce_i is the list of all cell entities of dimension i + 1. rc_i is the connectivity information
                // of dimension i + 1 (because we do skip(1))
                for (e_i, ce_i, rc_i) in izip!(
                    entities.iter(),
                    cell_entities.iter_mut(),
                    ref_conn.iter().skip(1),
                ) {
                    // This iterates over the actual sub entities.
                    // - ce_ij is the index of the jth subentity of dimension i + 1 (that's what we want to get in this loop).
                    // - rc_ij is the reference connectivity information of the ith subentity of dimension i + 1.
                    for (ce_ij, rc_ij) in izip!(ce_i.iter_mut(), rc_i) {
                        // We get the actual entity by mapping reference cell vertices to actual vertices
                        let mut entity = rc_ij[0].iter().map(|i| cell[*i]).collect::<Vec<_>>();
                        entity.sort();
                        *ce_ij = *e_i.get(&entity).unwrap();
                    }
                }
                // The following is a double loop that enters all k-dimensional entities connected to all
                // i-dimensional entities into the downward connectivities.
                // Copy these indices into connectivity for each dim
                // Loop over dim0
                // We now fill out the downward connectivity.
                // dc_i is downward connectivity for dim i.
                // rc_i is the reference connectivity for dim i.
                // ce_i are the cell entitiies for dim i.
                // The loop ignores dim 0 and dim 1 entities as everything for dim 1 entities
                // has been processed at this point (entities with themselves have been done, and entities with vertices).
                // The index i can only be 0 or one. We have 4 dimensions in total (0, 1, 2, 3). If we skip 0 and 1 then
                // the possible dimensions are 2 and 3. Since i is a counter i=0 is associated with dimension 2 and i=1 is
                // associated with dimension 3.
                for (i, (dc_i, rc_i, ce_i)) in izip!(
                    downward_connectivity.iter_mut().skip(2),
                    ref_conn.iter().skip(2),
                    cell_entities.iter().skip(1)
                )
                .enumerate()
                {
                    // ce_ij is the jth subentity of dimension i.
                    // rc_ij is the corresponding reference connectivity of the jth subentity of dimension i.
                    for (ce_ij, rc_ij) in izip!(ce_i, rc_i) {
                        // k loops over the dimensions of the cell entities.
                        // rc_ijk is the kth dimensional reference connectivity of the jth subentity of dimension i.
                        // dc_ik is the kth dimensional downward connectivity of the subentities of dimension i.
                        // If i = 0 we consider 2-dimensional entities (since we have skipped vertices and edges). Hence,
                        // dc_i.iter_mut().take(i+2) takes the first two dimensions of the downward connectivity,
                        // that is vertices and edges. It then skips over the vertices. So k = 0 is only option.
                        // If i = 1 we consider 3-dimensional entities. dc_i.iter_mut().take(i+2) takes the first 3
                        // dimensions, that is vertices, edges, faces. It then skips over the vertices. So k can be 0 and 1.
                        for (k, (dc_ik, rc_ijk, ce_k)) in izip!(
                            dc_i.iter_mut().take(i + 2).skip(1),
                            rc_ij.iter().take(i + 2).skip(1),
                            &cell_entities
                        )
                        .enumerate()
                        {
                            // l is the concrete index of the kth dimensional connectivity of the jth subentity of dimension i.
                            // We copy this into the right position of the downward connectivity.
                            // Remember that ce_ij is the jth subentity of dimension i.
                            for (l, rc_ijkl) in rc_ijk.iter().enumerate() {
                                dc_ik[[l, *ce_ij]] = ce_k[*rc_ijkl];
                                // Now will in reverse the corresponding upward connectivity.
                                // If i = 0 then we consider 2-dimensional entities (faces). k=0 corresponds to edges.
                                // The upward connectivity between the two is upward_connectivity[1][0].
                                // The general setup is upward_connectivity[dim0][dim1 - dim0 - 1][dim0_entity_index][..] = [dim1_entity_index].
                                // We hae dim0 = k + 1 and dim1 = i + 2. Hence, dim1 - dim0 - 1 = i - k.
                                if !upward_connectivity[k + 1][i - k][ce_k[*rc_ijkl]]
                                    .contains(ce_ij)
                                {
                                    upward_connectivity[k + 1][i - k][ce_k[*rc_ijkl]].push(*ce_ij);
                                }
                            }
                        }
                    }
                }
            }
        }

        let mut ids = vec![vertex_ids];
        for _ in 1..dim {
            ids.push(None);
        }
        ids.push(cell_ids);

        let ids_to_indices = ids
            .iter()
            .map(|ids_option| {
                let mut ie = HashMap::new();
                if let Some(ids) = ids_option {
                    for (i, j) in ids.iter().enumerate() {
                        ie.insert(*j, i);
                    }
                }
                ie
            })
            .collect::<Vec<_>>();

        let mut orientation = vec![];

        for d in 1..dim + 1 {
            let dc = &downward_connectivity[d][0];
            orientation.push(
                (0..dc.shape()[1])
                    .map(|i| {
                        compute_orientation(
                            entity_types[d],
                            &dc.data().unwrap()[i * dc.shape()[0]..(i + 1) * dc.shape()[0]],
                        )
                    })
                    .collect::<Vec<_>>(),
            );
        }

        Self {
            dim,
            ids,
            ids_to_indices,
            entity_types,
            entity_counts,
            downward_connectivity,
            upward_connectivity,
            orientation,
        }
    }
    /// Topological dimension
    pub fn dim(&self) -> usize {
        self.dim
    }
    /// Entity types
    pub fn entity_types(&self) -> &[ReferenceCellType] {
        &self.entity_types
    }
    /// Entity count
    pub fn entity_counts(&self) -> &[usize] {
        &self.entity_counts
    }
    /// Entity count
    pub fn entity_count(&self, entity_type: ReferenceCellType) -> usize {
        if !self.entity_types.contains(&entity_type) {
            0
        } else {
            self.entity_counts[reference_cell::dim(entity_type)]
        }
    }
    /// Downward connectivity
    pub fn downward_connectivity(&self) -> &[Vec<DynArray<usize, 2>>] {
        &self.downward_connectivity
    }
    /// Upward connectivity
    pub fn upward_connectivity(&self) -> &[Vec<Vec<Vec<usize>>>] {
        &self.upward_connectivity
    }
    /// Cell sub-entity index
    pub fn cell_entity_index(
        &self,
        cell_index: usize,
        entity_dim: usize,
        entity_index: usize,
    ) -> usize {
        self.downward_connectivity[self.dim][entity_dim][[entity_index, cell_index]]
    }
    /// Entity id
    pub fn entity_id(&self, entity_dim: usize, entity_index: usize) -> Option<usize> {
        self.ids[entity_dim].as_ref().map(|a| a[entity_index])
    }
}

/// Topology of a cell
#[derive(Debug)]
pub struct SingleTypeEntityTopology<'a> {
    topology: &'a SingleTypeTopology,
    entity_index: usize,
    dim: usize,
}

impl<'t> SingleTypeEntityTopology<'t> {
    /// Create new
    pub fn new(
        topology: &'t SingleTypeTopology,
        entity_type: ReferenceCellType,
        entity_index: usize,
    ) -> Self {
        Self {
            topology,
            entity_index,
            dim: reference_cell::dim(entity_type),
        }
    }
}
impl Topology for SingleTypeEntityTopology<'_> {
    type EntityDescriptor = ReferenceCellType;
    type EntityIndexIter<'a>
        = Copied<std::slice::Iter<'a, usize>>
    where
        Self: 'a;

    type ConnectedEntityIndexIter<'a>
        = Copied<std::slice::Iter<'a, usize>>
    where
        Self: 'a;

    fn connected_entity_iter(
        &self,
        entity_type: ReferenceCellType,
    ) -> Copied<std::slice::Iter<'_, usize>> {
        let dim = reference_cell::dim(entity_type);
        if entity_type == self.topology.entity_types()[dim] {
            self.topology.upward_connectivity[self.dim][dim - self.dim - 1][self.entity_index]
                .iter()
                .copied()
        } else {
            [].iter().copied()
        }
    }

    fn sub_entity_iter(
        &self,
        entity_type: ReferenceCellType,
    ) -> Copied<std::slice::Iter<'_, usize>> {
        let dim = reference_cell::dim(entity_type);
        if entity_type == self.topology.entity_types()[dim] {
            let rows = self.topology.downward_connectivity[self.dim][dim].shape()[0];
            self.topology.downward_connectivity[self.dim][dim]
                .data()
                .unwrap()[rows * self.entity_index..rows * (self.entity_index + 1)]
                .iter()
                .copied()
        } else {
            [].iter().copied()
        }
    }

    fn sub_entity(&self, entity_type: ReferenceCellType, index: usize) -> usize {
        let dim = reference_cell::dim(entity_type);
        if entity_type != self.topology.entity_types()[dim] {
            panic!("Invalid entity type");
        }
        self.topology.downward_connectivity[self.dim][dim][[index, self.entity_index]]
    }

    fn orientation(&self) -> i32 {
        if self.dim == 0 {
            0
        } else {
            self.topology.orientation[self.dim - 1][self.entity_index]
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn example_topology_point() -> SingleTypeTopology {
        //! An example topology
        SingleTypeTopology::new(&[0, 1], ReferenceCellType::Point, None, None)
    }

    fn example_topology_interval() -> SingleTypeTopology {
        //! An example topology
        SingleTypeTopology::new(&[0, 1, 1, 2], ReferenceCellType::Interval, None, None)
    }

    fn example_topology_triangle() -> SingleTypeTopology {
        //! An example topology
        SingleTypeTopology::new(&[0, 1, 2, 2, 1, 3], ReferenceCellType::Triangle, None, None)
    }

    fn example_topology_tetrahedron() -> SingleTypeTopology {
        //! An example topology
        SingleTypeTopology::new(
            &[0, 1, 2, 3, 4, 0, 2, 3],
            ReferenceCellType::Tetrahedron,
            None,
            None,
        )
    }

    #[test]
    #[should_panic]
    fn test_prism() {
        let _ = SingleTypeTopology::new(&[0, 1, 2, 3, 4, 5], ReferenceCellType::Prism, None, None);
    }
    #[test]
    #[should_panic]
    fn test_pyramid() {
        let _ = SingleTypeTopology::new(&[0, 1, 2, 3, 4], ReferenceCellType::Pyramid, None, None);
    }

    #[test]
    fn test_sub_entities_point() {
        //! Test sub-entities of a point
        let t = example_topology_point();
        let cell0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 0), 0);
        let cell1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 1);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 0), 1);
    }

    #[test]
    fn test_sub_entities_interval() {
        //! Test sub-entities of an interval
        let t = example_topology_interval();
        let cell0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 1), 1);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 0), 0);
        let cell1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 1);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 0), 1);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 1), 2);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 0), 1);

        let vertex0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 0);
        assert_eq!(vertex0.sub_entity(ReferenceCellType::Point, 0), 0);
        let vertex1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 1);
        assert_eq!(vertex1.sub_entity(ReferenceCellType::Point, 0), 1);
        let vertex2 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 2);
        assert_eq!(vertex2.sub_entity(ReferenceCellType::Point, 0), 2);
    }

    #[test]
    fn test_sub_entities_triangle() {
        //! Test sub-entities of a triangle
        let t = example_topology_triangle();
        let cell0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 1), 1);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 2), 2);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 0), 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 1), 1);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 2), 2);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Triangle, 0), 0);
        let cell1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 1);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 0), 2);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 1), 1);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 2), 3);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 0), 3);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 1), 4);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 2), 0);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Triangle, 0), 1);

        let edge0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 0);
        assert_eq!(edge0.sub_entity(ReferenceCellType::Point, 0), 1);
        assert_eq!(edge0.sub_entity(ReferenceCellType::Point, 1), 2);
        assert_eq!(edge0.sub_entity(ReferenceCellType::Interval, 0), 0);
        let edge1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 1);
        assert_eq!(edge1.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(edge1.sub_entity(ReferenceCellType::Point, 1), 2);
        assert_eq!(edge1.sub_entity(ReferenceCellType::Interval, 0), 1);
        let edge2 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 2);
        assert_eq!(edge2.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(edge2.sub_entity(ReferenceCellType::Point, 1), 1);
        assert_eq!(edge2.sub_entity(ReferenceCellType::Interval, 0), 2);
        let edge3 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 3);
        assert_eq!(edge3.sub_entity(ReferenceCellType::Point, 0), 1);
        assert_eq!(edge3.sub_entity(ReferenceCellType::Point, 1), 3);
        assert_eq!(edge3.sub_entity(ReferenceCellType::Interval, 0), 3);
        let edge4 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 4);
        assert_eq!(edge4.sub_entity(ReferenceCellType::Point, 0), 2);
        assert_eq!(edge4.sub_entity(ReferenceCellType::Point, 1), 3);
        assert_eq!(edge4.sub_entity(ReferenceCellType::Interval, 0), 4);

        let vertex0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 0);
        assert_eq!(vertex0.sub_entity(ReferenceCellType::Point, 0), 0);
        let vertex1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 1);
        assert_eq!(vertex1.sub_entity(ReferenceCellType::Point, 0), 1);
        let vertex2 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 2);
        assert_eq!(vertex2.sub_entity(ReferenceCellType::Point, 0), 2);
        let vertex3 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 3);
        assert_eq!(vertex3.sub_entity(ReferenceCellType::Point, 0), 3);
    }

    #[test]
    fn test_sub_entities_tetrahedron() {
        //! Test sub-entities of a tetrahedron
        let t = example_topology_tetrahedron();
        let cell0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Tetrahedron, 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 1), 1);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 2), 2);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 3), 3);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 0), 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 1), 1);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 2), 2);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 3), 3);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 4), 4);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 5), 5);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Triangle, 0), 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Triangle, 1), 1);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Triangle, 2), 2);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Triangle, 3), 3);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Tetrahedron, 0), 0);
        let cell1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Tetrahedron, 1);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 0), 4);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 1), 0);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 2), 2);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 3), 3);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 0), 0);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 1), 3);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 2), 4);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 3), 6);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 4), 7);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 5), 8);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Triangle, 0), 1);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Triangle, 1), 4);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Triangle, 2), 5);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Triangle, 3), 6);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Tetrahedron, 0), 1);

        let face0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 0);
        assert_eq!(face0.sub_entity(ReferenceCellType::Point, 0), 1);
        assert_eq!(face0.sub_entity(ReferenceCellType::Point, 1), 2);
        assert_eq!(face0.sub_entity(ReferenceCellType::Point, 2), 3);
        assert_eq!(face0.sub_entity(ReferenceCellType::Interval, 0), 0);
        assert_eq!(face0.sub_entity(ReferenceCellType::Interval, 1), 1);
        assert_eq!(face0.sub_entity(ReferenceCellType::Interval, 2), 2);
        assert_eq!(face0.sub_entity(ReferenceCellType::Triangle, 0), 0);
        let face1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 1);
        assert_eq!(face1.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(face1.sub_entity(ReferenceCellType::Point, 1), 2);
        assert_eq!(face1.sub_entity(ReferenceCellType::Point, 2), 3);
        assert_eq!(face1.sub_entity(ReferenceCellType::Interval, 0), 0);
        assert_eq!(face1.sub_entity(ReferenceCellType::Interval, 1), 3);
        assert_eq!(face1.sub_entity(ReferenceCellType::Interval, 2), 4);
        assert_eq!(face1.sub_entity(ReferenceCellType::Triangle, 0), 1);
        let face2 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 2);
        assert_eq!(face2.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(face2.sub_entity(ReferenceCellType::Point, 1), 1);
        assert_eq!(face2.sub_entity(ReferenceCellType::Point, 2), 3);
        assert_eq!(face2.sub_entity(ReferenceCellType::Interval, 0), 1);
        assert_eq!(face2.sub_entity(ReferenceCellType::Interval, 1), 3);
        assert_eq!(face2.sub_entity(ReferenceCellType::Interval, 2), 5);
        assert_eq!(face2.sub_entity(ReferenceCellType::Triangle, 0), 2);
        let face3 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 3);
        assert_eq!(face3.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(face3.sub_entity(ReferenceCellType::Point, 1), 1);
        assert_eq!(face3.sub_entity(ReferenceCellType::Point, 2), 2);
        assert_eq!(face3.sub_entity(ReferenceCellType::Interval, 0), 2);
        assert_eq!(face3.sub_entity(ReferenceCellType::Interval, 1), 4);
        assert_eq!(face3.sub_entity(ReferenceCellType::Interval, 2), 5);
        assert_eq!(face3.sub_entity(ReferenceCellType::Triangle, 0), 3);
        let face4 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 4);
        assert_eq!(face4.sub_entity(ReferenceCellType::Point, 0), 2);
        assert_eq!(face4.sub_entity(ReferenceCellType::Point, 1), 3);
        assert_eq!(face4.sub_entity(ReferenceCellType::Point, 2), 4);
        assert_eq!(face4.sub_entity(ReferenceCellType::Interval, 0), 0);
        assert_eq!(face4.sub_entity(ReferenceCellType::Interval, 1), 6);
        assert_eq!(face4.sub_entity(ReferenceCellType::Interval, 2), 7);
        assert_eq!(face4.sub_entity(ReferenceCellType::Triangle, 0), 4);
        let face5 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 5);
        assert_eq!(face5.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(face5.sub_entity(ReferenceCellType::Point, 1), 3);
        assert_eq!(face5.sub_entity(ReferenceCellType::Point, 2), 4);
        assert_eq!(face5.sub_entity(ReferenceCellType::Interval, 0), 3);
        assert_eq!(face5.sub_entity(ReferenceCellType::Interval, 1), 6);
        assert_eq!(face5.sub_entity(ReferenceCellType::Interval, 2), 8);
        assert_eq!(face5.sub_entity(ReferenceCellType::Triangle, 0), 5);
        let face6 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 6);
        assert_eq!(face6.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(face6.sub_entity(ReferenceCellType::Point, 1), 2);
        assert_eq!(face6.sub_entity(ReferenceCellType::Point, 2), 4);
        assert_eq!(face6.sub_entity(ReferenceCellType::Interval, 0), 4);
        assert_eq!(face6.sub_entity(ReferenceCellType::Interval, 1), 7);
        assert_eq!(face6.sub_entity(ReferenceCellType::Interval, 2), 8);
        assert_eq!(face6.sub_entity(ReferenceCellType::Triangle, 0), 6);

        let edge0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 0);
        assert_eq!(edge0.sub_entity(ReferenceCellType::Point, 0), 2);
        assert_eq!(edge0.sub_entity(ReferenceCellType::Point, 1), 3);
        assert_eq!(edge0.sub_entity(ReferenceCellType::Interval, 0), 0);
        let edge1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 1);
        assert_eq!(edge1.sub_entity(ReferenceCellType::Point, 0), 1);
        assert_eq!(edge1.sub_entity(ReferenceCellType::Point, 1), 3);
        assert_eq!(edge1.sub_entity(ReferenceCellType::Interval, 0), 1);
        let edge2 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 2);
        assert_eq!(edge2.sub_entity(ReferenceCellType::Point, 0), 1);
        assert_eq!(edge2.sub_entity(ReferenceCellType::Point, 1), 2);
        assert_eq!(edge2.sub_entity(ReferenceCellType::Interval, 0), 2);
        let edge3 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 3);
        assert_eq!(edge3.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(edge3.sub_entity(ReferenceCellType::Point, 1), 3);
        assert_eq!(edge3.sub_entity(ReferenceCellType::Interval, 0), 3);
        let edge4 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 4);
        assert_eq!(edge4.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(edge4.sub_entity(ReferenceCellType::Point, 1), 2);
        assert_eq!(edge4.sub_entity(ReferenceCellType::Interval, 0), 4);
        let edge5 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 5);
        assert_eq!(edge5.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(edge5.sub_entity(ReferenceCellType::Point, 1), 1);
        assert_eq!(edge5.sub_entity(ReferenceCellType::Interval, 0), 5);
        let edge6 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 6);
        assert_eq!(edge6.sub_entity(ReferenceCellType::Point, 0), 3);
        assert_eq!(edge6.sub_entity(ReferenceCellType::Point, 1), 4);
        assert_eq!(edge6.sub_entity(ReferenceCellType::Interval, 0), 6);
        let edge7 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 7);
        assert_eq!(edge7.sub_entity(ReferenceCellType::Point, 0), 2);
        assert_eq!(edge7.sub_entity(ReferenceCellType::Point, 1), 4);
        assert_eq!(edge7.sub_entity(ReferenceCellType::Interval, 0), 7);
        let edge8 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 8);
        assert_eq!(edge8.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(edge8.sub_entity(ReferenceCellType::Point, 1), 4);
        assert_eq!(edge8.sub_entity(ReferenceCellType::Interval, 0), 8);

        let vertex0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 0);
        assert_eq!(vertex0.sub_entity(ReferenceCellType::Point, 0), 0);
        let vertex1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 1);
        assert_eq!(vertex1.sub_entity(ReferenceCellType::Point, 0), 1);
        let vertex2 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 2);
        assert_eq!(vertex2.sub_entity(ReferenceCellType::Point, 0), 2);
        let vertex3 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 3);
        assert_eq!(vertex3.sub_entity(ReferenceCellType::Point, 0), 3);
        let vertex4 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 4);
        assert_eq!(vertex4.sub_entity(ReferenceCellType::Point, 0), 4);
    }

    macro_rules! make_tests {
        ($cellname:ident) => {
            paste::item! {
                #[test]
                fn [< test_up_and_down_ $cellname >]() {
                    //! Test that upward and downward connectivities agree
                    let t = [< example_topology_ $cellname >]();

                    for (i, dc_i) in t.downward_connectivity.iter().enumerate() {
                        for (j, dc_ij) in dc_i.iter().enumerate() {
                            if i != j {
                                let uc_ji = &t.upward_connectivity[j][i - j - 1];
                                for (c, col) in dc_ij.col_iter().enumerate() {
                                    for value in col.iter_ref() {
                                        assert!(uc_ji[*value].contains(&c));
                                    }
                                }
                                for (k, uc_jik) in uc_ji.iter().enumerate() {
                                    for value in uc_jik {
                                        assert!(dc_ij.r().slice::<1>(1, *value).iter_ref()
                                            .position(|&i| i == k).is_some());
                                    }
                                }
                            }
                        }
                    }
                }
                #[test]
                fn [< test_sub_entity_iter_ $cellname >]() {
                    //! Test sub-entity iterators
                    let t = [< example_topology_ $cellname >]();
                    for cell_type in t.entity_types() {
                        for index in 0..t.entity_count(*cell_type) {
                            let cell = SingleTypeEntityTopology::new(&t, *cell_type, index);
                            for dim in 0..reference_cell::dim(*cell_type) + 1 {
                                let sub_entity_type = t.entity_types()[dim];
                                for (i, j) in cell.sub_entity_iter(sub_entity_type).enumerate() {
                                    assert_eq!(j, cell.sub_entity(sub_entity_type, i));
                                }
                            }
                        }
                    }
                }
            }
        };
    }

    make_tests!(point);
    make_tests!(interval);
    make_tests!(triangle);
    make_tests!(tetrahedron);

    #[test]
    fn test_orientation_triangle() {
        let t =
            SingleTypeTopology::new(&[0, 1, 2, 0, 2, 1], ReferenceCellType::Triangle, None, None);
        assert_ne!(t.orientation[1][0], t.orientation[1][1]);
    }

    #[test]
    fn test_orientation_quadrilateral() {
        let t = SingleTypeTopology::new(
            &[0, 1, 2, 3, 0, 3, 1, 2],
            ReferenceCellType::Quadrilateral,
            None,
            None,
        );
        assert_ne!(t.orientation[1][0], t.orientation[1][1]);
    }
}
