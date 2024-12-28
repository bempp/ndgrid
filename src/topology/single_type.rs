//! Topology for grids where entities of each tdim are the same type

#[cfg(feature = "serde")]
use crate::traits::ConvertToSerializable;
use crate::traits::Topology;
use crate::types::{Array2D, Array2DBorrowed};
use itertools::izip;
use ndelement::reference_cell;
use ndelement::types::ReferenceCellType;
#[cfg(feature = "serde")]
use rlst::RawAccessMut;
use rlst::{rlst_dynamic_array2, DefaultIteratorMut, RawAccess, Shape};
use std::collections::{HashMap, HashSet};
use std::iter::Copied;
use std::sync::Arc;
use std::time::Instant;

/// Topology of a single element grid
#[derive(Debug)]
pub struct SingleTypeTopology {
    dim: usize,
    pub(crate) ids: Vec<Option<Vec<usize>>>,
    entity_types: Vec<ReferenceCellType>,
    entity_counts: Vec<usize>,
    pub(crate) downward_connectivity: Vec<Vec<Array2D<usize>>>,
    pub(crate) upward_connectivity: Vec<Vec<Vec<Vec<usize>>>>,
}

#[cfg(feature = "serde")]
#[derive(serde::Serialize, Debug, serde::Deserialize)]
/// Serde serializable topology
pub struct SerializableTopology {
    dim: usize,
    ids: Vec<Option<Vec<usize>>>,
    entity_types: Vec<ReferenceCellType>,
    entity_counts: Vec<usize>,
    downward_connectivity: Vec<Vec<(Vec<usize>, [usize; 2])>>,
    upward_connectivity: Vec<Vec<Vec<Vec<usize>>>>,
}

#[cfg(feature = "serde")]
impl ConvertToSerializable for SingleTypeTopology {
    type SerializableType = SerializableTopology;
    fn to_serializable(&self) -> SerializableTopology {
        SerializableTopology {
            dim: self.dim,
            ids: self.ids.clone(),
            entity_types: self.entity_types.clone(),
            entity_counts: self.entity_counts.clone(),
            downward_connectivity: self
                .downward_connectivity
                .iter()
                .map(|a| {
                    a.iter()
                        .map(|b| (b.data().to_vec(), b.shape()))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>(),
            upward_connectivity: self.upward_connectivity.clone(),
        }
    }
    fn from_serializable(s: SerializableTopology) -> Self {
        Self {
            dim: s.dim,
            ids: s.ids,
            entity_types: s.entity_types,
            entity_counts: s.entity_counts,
            downward_connectivity: s
                .downward_connectivity
                .iter()
                .map(|a| {
                    a.iter()
                        .map(|(data, shape)| {
                            let mut c = rlst_dynamic_array2!(usize, *shape);
                            c.data_mut().copy_from_slice(data);
                            c
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>(),
            upward_connectivity: s.upward_connectivity,
        }
    }
}

// TODO: delete following line
// unsafe impl Sync for SingleTypeTopology {}
// Entities are reoriented so that the smallest vertex index comes first.
fn orient_entity(entity_type: ReferenceCellType, vertices: &mut [usize]) {
    match entity_type {
        ReferenceCellType::Point => {}
        ReferenceCellType::Interval => {
            if vertices[0] > vertices[1] {
                vertices.swap(0, 1);
            }
        }
        ReferenceCellType::Triangle => {
            if vertices[0] > vertices[1] {
                vertices.swap(0, 1);
            }
            if vertices[1] > vertices[2] {
                vertices.swap(1, 2);
            }
            if vertices[0] > vertices[1] {
                vertices.swap(0, 1);
            }
        }
        ReferenceCellType::Quadrilateral => {
            let minimum = *vertices.iter().min().unwrap();
            if vertices[1] == minimum {
                vertices.swap(0, 1);
                vertices.swap(2, 3);
            } else if vertices[2] == minimum {
                vertices.swap(0, 2);
                vertices.swap(1, 3);
            } else if vertices[3] == minimum {
                vertices.swap(0, 3);
            }
            if vertices[1] > vertices[2] {
                vertices.swap(1, 2);
            }
        }
        _ => {
            unimplemented!();
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

        //The reference connectivity is a 4-dimensional array.
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
        for _ in 0..dim - 1 {
            entities.push(HashMap::new());
        }
        // Old inefficient entity definition.
        //let mut entities = vec![vec![]; if dim == 0 { 0 } else { dim - 1 }];

        // Now iterate through all the cells
        for cell_index in 0..ncells {
            // Get all the vertices of the cell.
            let cell = &cells[cell_index * size..(cell_index + 1) * size];
            // Iterate over topological dimension.
            // The variable e_i is a mutable reference to the list of entities of dimension i.
            // The variable rc_i is a reference to the connectivity of the cell type for entities of dimension i.
            // Hence, rc_i[0] is the connectivity information about the entity with reference index 0 and
            // dimension i and et_i[0] gives the corresponding entity type.
            // etypes is a reference to the entity types of the cell type at dimension i.
            for (e_i, rc_i, et_i) in izip!(
                entities.iter_mut(),
                // Skip the vertices in the connectivity
                ref_conn.iter().take(dim).skip(1),
                // Skip the vertices in the entity types.
                etypes.iter().take(dim).skip(1)
            ) {
                // A simple counter for an entity index of dimension i.
                let mut entity_counter = 0;
                // We iterate through the concrete entities of dimension i and the corresponding
                // entity types. c_ij is a reference to the connectivity information of the jth entity
                // with dimension i. et_ij is the corresponding entity type.
                for (c_ij, et_ij) in izip!(rc_i, et_i) {
                    // c_ij[0] is the list of reference cell vertex indices that are associated with the jth entity of dimension i.
                    // cell[*i] below maps the local reference cell vertex index to the actual id of the vertex.
                    // Hence, the following command gives us all the vertex indicies of the entity.
                    let mut entity = c_ij[0].iter().map(|i| cell[*i]).collect::<Vec<_>>();
                    // We reorient entities so that identities with same vertices but different order of vertices
                    // are treated as the same entity.
                    orient_entity(*et_ij, &mut entity);
                    // An entity only gets added if it has not been added before.
                    e_i.entry(entity).or_insert_with(|| {
                        let old = entity_counter;
                        entity_counter += 1;
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

        println!("In new topology 1");
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
                    .map(|r| rlst_dynamic_array2!(usize, [r.len(), *j]))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        println!("In new topology2");
        // Upward connectivity: The entities of dimension dim1 that are superentities of
        // entities of dimension dim0 (with dim0<dim1) (eg triangles connected to an edge,
        // tetrahedra connected to a vertex, etc)
        // upward_connectivity[dim0][dim1 - dim0 - 1][dim0_entity_index][..] = [dim1_entity_index]
        // TODO: The last part of upwards connectivity should be a set.
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

        println!("In new topology3");
        // downward_connectivity[dim][0] = vertices of each cell
        // This is just the vertices that each cell is made of.
        for (i, mut dc_d0i) in downward_connectivity[dim][0].col_iter_mut().enumerate() {
            for (dc_d0ij, c_j) in izip!(dc_d0i.iter_mut(), &cells[i * size..(i + 1) * size]) {
                *dc_d0ij = *c_j;
                // We can also fill out the upward connectivity if dim > 0.
                // The contains test should not be necessary as for eaach vertex a cell can only
                // be added once. However, it seems better to move upward connectivity to sets anyway.
                if dim > 0 && !upward_connectivity[0][dim - 1][*c_j].contains(&i) {
                    upward_connectivity[0][dim - 1][*c_j].push(i);
                }
            }
        }

        // We iterate through all the cells
        for cell_index in 0..ncells {
            // We get the cell vertices.
            let cell = &cells[cell_index * size..(cell_index + 1) * size];
            for (i, (e_i, rc_i, et_i, dc_i)) in izip!(
                entities.iter_mut(),
                // Skip the vertices in the connectivity
                ref_conn.iter().take(dim).skip(1),
                // Skip the vertices in the entity types.
                etypes.iter().take(dim).skip(1),
                downward_connectivity.iter_mut().take(dim).skip(1)
            )
            .enumerate()
            {
                for (c_ij, et_ij) in izip!(rc_i, et_i) {
                    // c_ij[0] is the list of reference cell vertex indices that are associated with the jth entity of dimension i.
                    // cell[*i] below maps the local reference cell vertex index to the actual id of the vertex.
                    // Hence, the following command gives us all the vertex indicies of the entity.
                    let mut entity = c_ij[0].iter().map(|i| cell[*i]).collect::<Vec<_>>();
                    // We reorient entities so that identities with same vertices but different order of vertices
                    // are treated as the same entity.
                    orient_entity(*et_ij, &mut entity);
                    let entity_index = *e_i.get(&entity).unwrap();

                    // We now have entity and entity index, so can fill up the downward and upward connectivity.
                    assert_eq!(dc_i[0].r().slice(1, entity_index).len(), entity.len());
                    for (dc_i_vertex, &entity_vertex) in izip!(
                        dc_i[0].r_mut().slice(1, entity_index).iter_mut(),
                        entity.iter()
                    ) {
                        *dc_i_vertex = entity_vertex;
                        if !upward_connectivity[0][i][entity_vertex].contains(&entity_index) {
                            upward_connectivity[0][i][entity_index].push(entity_index);
                        }
                    }
                }
            }
        }

        // println!("In new topology4");
        // // downward_connectivity[i][0] = vertices of entity
        // // This is similar to before. But now we actually add the vertices of each entity.
        // // We iterate through the entities and the downward connectivity. Since the entities
        // // don't have information about dim 0 we skip the first entry in the downward connectivity.
        // // i + 1 is our current dimension. es_i is the list of entities of dimension i in the
        // // reference cell, and dc_i is the corresponding iterator over downward connectivities.
        // for (i, (es_i, dc_i)) in izip!(
        //     entities.iter(),
        //     downward_connectivity.iter_mut().take(dim).skip(1)
        // )
        // .enumerate()
        // {
        //     // We are now iterating over the entities of dimension i + 1. es_ij is the jth entity of dimension i + 1.
        //     // dc_i0j is a column vector that contains the vertices that make up the jth entity of dimension i +  1.
        //     for (j, (es_ij, mut dc_i0j)) in izip!(es_i.iter(), dc_i[0].col_iter_mut()).enumerate() {
        //         // We are now copying the vertices from the entity es_ij into the column vector dc_i0j.
        //         for (es_ijk, dc_i0jk) in izip!(es_ij.iter(), dc_i0j.iter_mut()) {
        //             *dc_i0jk = *es_ijk;
        //             // In each iteration we test for the upward connectivity whether the entity j is already in the list
        //             // of upward connectivities for the vertex es_ijk. If not we add it.
        //             // Again, I don't see how this could require a set as this is a unique mapping.
        //             if !upward_connectivity[0][i][*es_ijk].contains(&j) {
        //                 upward_connectivity[0][i][*es_ijk].push(j);
        //             }
        //         }
        //     }
        // }

        println!("In new topology5");
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
                // of dimension i + 1 (because we do skip(1)) and etypes is the corresponding reference cell type
                // (again with skip(1)).
                for (e_i, ce_i, rc_i, et_i) in izip!(
                    entities.iter(),
                    cell_entities.iter_mut(),
                    ref_conn.iter().skip(1),
                    etypes.iter().skip(1)
                ) {
                    // This iterates over the actual sub entities.
                    // - et_ij is the jth subentity of dimension i + 1.
                    // - ce_ij is the index of the jth subentity of dimension i + 1 (that's what we want to get in this loop).
                    // - rc_ij is the reference connectivity information of the ith subentity of dimension i + 1.
                    for (ce_ij, rc_ij, et_ij) in izip!(ce_i.iter_mut(), rc_i, et_i) {
                        // We get the actual entity by mapping reference cell vertices to actual vertices
                        // and then reorienting.
                        let mut entity = rc_ij[0].iter().map(|i| cell[*i]).collect::<Vec<_>>();
                        orient_entity(*et_ij, &mut entity);
                        // AAARGGHH!!! This is expensive!!. We should use a hash map here.
                        // This finds the index of the ith entity by using `position` method on the
                        // entity list.
                        *ce_ij = *e_i.get(&entity).unwrap();
                        // Old inefficient version when entities was an array.
                        //*ce_ij = e_i.iter().position(|r| *r == entity).unwrap();
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
                        // If i = 0 we consider 1-dimensional entities (since we have skipped vertices). Hence,
                        // dc_i.iter_mut().take(i+2) takes the first two dimensions of the downward connectivity,
                        // that is vertices and edges. It then skips over the vertices.
                        // If i = 2 we consider 3-dimensional entities. dc_i.iter_mut().take(i+2) takes the first 4
                        // dimensions, that is vertices, edges, faces, and volumes. It then skips over the vertices.
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

        println!("Finished topology loop.");
        let mut ids = vec![vertex_ids];
        for _ in 1..dim {
            ids.push(None);
        }
        ids.push(cell_ids);
        Self {
            dim,
            ids,
            entity_types,
            entity_counts,
            downward_connectivity,
            upward_connectivity,
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
    pub fn downward_connectivity(&self) -> &[Vec<Array2D<usize>>] {
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
    //entity_type: ReferenceCellType,
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
            //entity_type,
            entity_index,
            dim: reference_cell::dim(entity_type),
        }
    }
}
impl Topology for SingleTypeEntityTopology<'_> {
    type EntityIndexIter<'a>
        = Copied<std::slice::Iter<'a, usize>>
    where
        Self: 'a;

    type ConnectedEntityIndexIter<'a>
        = Copied<std::slice::Iter<'a, usize>>
    where
        Self: 'a;

    fn connected_entity_iter(&self, dim: usize) -> Copied<std::slice::Iter<'_, usize>> {
        self.topology.upward_connectivity[self.dim][dim - self.dim - 1][self.entity_index]
            .iter()
            .copied()
    }

    fn sub_entity_iter(&self, dim: usize) -> Copied<std::slice::Iter<'_, usize>> {
        let rows = self.topology.downward_connectivity[self.dim][dim].shape()[0];
        self.topology.downward_connectivity[self.dim][dim].data()
            [rows * self.entity_index..rows * (self.entity_index + 1)]
            .iter()
            .copied()
    }

    fn sub_entity(&self, dim: usize, index: usize) -> usize {
        self.topology.downward_connectivity[self.dim][dim][[index, self.entity_index]]
    }
}

/// Topology of a single element grid with borrowed data
#[derive(Debug)]
pub struct SingleTypeTopologyBorrowed<'a> {
    dim: usize,
    pub(crate) ids: Vec<Option<&'a [usize]>>,
    entity_types: &'a [ReferenceCellType],
    entity_counts: &'a [usize],
    pub(crate) downward_connectivity: Vec<Vec<Array2DBorrowed<'a, usize>>>,
    pub(crate) upward_connectivity: Vec<Vec<Vec<&'a [usize]>>>,
}

impl<'a> SingleTypeTopologyBorrowed<'a> {
    /// Create new
    pub fn new(
        dim: usize,
        ids: Vec<Option<&'a [usize]>>,
        entity_types: &'a [ReferenceCellType],
        entity_counts: &'a [usize],
        downward_connectivity: Vec<Vec<Array2DBorrowed<'a, usize>>>,
        upward_connectivity: Vec<Vec<Vec<&'a [usize]>>>,
    ) -> Self {
        Self {
            dim,
            ids,
            entity_types,
            entity_counts,
            downward_connectivity,
            upward_connectivity,
        }
    }
    /// Topological dimension
    pub fn dim(&self) -> usize {
        self.dim
    }
    /// Entity types
    pub fn entity_types(&self) -> &[ReferenceCellType] {
        self.entity_types
    }
    /// Entity counts
    pub fn entity_count(&self, entity_type: ReferenceCellType) -> usize {
        if !self.entity_types.contains(&entity_type) {
            0
        } else {
            self.entity_counts[reference_cell::dim(entity_type)]
        }
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
pub struct SingleTypeEntityTopologyBorrowed<'a> {
    topology: &'a SingleTypeTopologyBorrowed<'a>,
    entity_index: usize,
    dim: usize,
}

impl<'t> SingleTypeEntityTopologyBorrowed<'t> {
    /// Create new
    pub fn new(
        topology: &'t SingleTypeTopologyBorrowed<'t>,
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
impl Topology for SingleTypeEntityTopologyBorrowed<'_> {
    type EntityIndexIter<'a>
        = Copied<std::slice::Iter<'a, usize>>
    where
        Self: 'a;

    type ConnectedEntityIndexIter<'a>
        = Copied<std::slice::Iter<'a, usize>>
    where
        Self: 'a;

    fn connected_entity_iter(&self, dim: usize) -> Copied<std::slice::Iter<'_, usize>> {
        self.topology.upward_connectivity[self.dim][dim - self.dim - 1][self.entity_index]
            .iter()
            .copied()
    }

    fn sub_entity_iter(&self, dim: usize) -> Copied<std::slice::Iter<'_, usize>> {
        let rows = self.topology.downward_connectivity[self.dim][dim].shape()[0];
        self.topology.downward_connectivity[self.dim][dim].data()
            [rows * self.entity_index..rows * (self.entity_index + 1)]
            .iter()
            .copied()
    }

    fn sub_entity(&self, dim: usize, index: usize) -> usize {
        self.topology.downward_connectivity[self.dim][dim][[index, self.entity_index]]
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rlst::DefaultIterator;

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
        let _ = SingleTypeTopology::new(&[0, 1, 2, 3, 4], ReferenceCellType::Prism, None, None);
    }

    #[test]
    fn test_sub_entities_point() {
        //! Test sub-entities of a point
        let t = example_topology_point();
        let cell0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 0);
        assert_eq!(cell0.sub_entity(0, 0), 0);
        let cell1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 1);
        assert_eq!(cell1.sub_entity(0, 0), 1);
    }

    #[test]
    fn test_sub_entities_interval() {
        //! Test sub-entities of an interval
        let t = example_topology_interval();
        let cell0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 0);
        assert_eq!(cell0.sub_entity(0, 0), 0);
        assert_eq!(cell0.sub_entity(0, 1), 1);
        assert_eq!(cell0.sub_entity(1, 0), 0);
        let cell1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 1);
        assert_eq!(cell1.sub_entity(0, 0), 1);
        assert_eq!(cell1.sub_entity(0, 1), 2);
        assert_eq!(cell1.sub_entity(1, 0), 1);

        let vertex0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 0);
        assert_eq!(vertex0.sub_entity(0, 0), 0);
        let vertex1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 1);
        assert_eq!(vertex1.sub_entity(0, 0), 1);
        let vertex2 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 2);
        assert_eq!(vertex2.sub_entity(0, 0), 2);
    }

    #[test]
    fn test_sub_entities_triangle() {
        //! Test sub-entities of a triangle
        let t = example_topology_triangle();
        let cell0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 0);
        assert_eq!(cell0.sub_entity(0, 0), 0);
        assert_eq!(cell0.sub_entity(0, 1), 1);
        assert_eq!(cell0.sub_entity(0, 2), 2);
        assert_eq!(cell0.sub_entity(1, 0), 0);
        assert_eq!(cell0.sub_entity(1, 1), 1);
        assert_eq!(cell0.sub_entity(1, 2), 2);
        assert_eq!(cell0.sub_entity(2, 0), 0);
        let cell1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 1);
        assert_eq!(cell1.sub_entity(0, 0), 2);
        assert_eq!(cell1.sub_entity(0, 1), 1);
        assert_eq!(cell1.sub_entity(0, 2), 3);
        assert_eq!(cell1.sub_entity(1, 0), 3);
        assert_eq!(cell1.sub_entity(1, 1), 4);
        assert_eq!(cell1.sub_entity(1, 2), 0);
        assert_eq!(cell1.sub_entity(2, 0), 1);

        let edge0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 0);
        assert_eq!(edge0.sub_entity(0, 0), 1);
        assert_eq!(edge0.sub_entity(0, 1), 2);
        assert_eq!(edge0.sub_entity(1, 0), 0);
        let edge1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 1);
        assert_eq!(edge1.sub_entity(0, 0), 0);
        assert_eq!(edge1.sub_entity(0, 1), 2);
        assert_eq!(edge1.sub_entity(1, 0), 1);
        let edge2 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 2);
        assert_eq!(edge2.sub_entity(0, 0), 0);
        assert_eq!(edge2.sub_entity(0, 1), 1);
        assert_eq!(edge2.sub_entity(1, 0), 2);
        let edge3 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 3);
        assert_eq!(edge3.sub_entity(0, 0), 1);
        assert_eq!(edge3.sub_entity(0, 1), 3);
        assert_eq!(edge3.sub_entity(1, 0), 3);
        let edge4 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 4);
        assert_eq!(edge4.sub_entity(0, 0), 2);
        assert_eq!(edge4.sub_entity(0, 1), 3);
        assert_eq!(edge4.sub_entity(1, 0), 4);

        let vertex0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 0);
        assert_eq!(vertex0.sub_entity(0, 0), 0);
        let vertex1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 1);
        assert_eq!(vertex1.sub_entity(0, 0), 1);
        let vertex2 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 2);
        assert_eq!(vertex2.sub_entity(0, 0), 2);
        let vertex3 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 3);
        assert_eq!(vertex3.sub_entity(0, 0), 3);
    }

    #[test]
    fn test_sub_entities_tetrahedron() {
        //! Test sub-entities of a tetrahedron
        let t = example_topology_tetrahedron();
        let cell0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Tetrahedron, 0);
        assert_eq!(cell0.sub_entity(0, 0), 0);
        assert_eq!(cell0.sub_entity(0, 1), 1);
        assert_eq!(cell0.sub_entity(0, 2), 2);
        assert_eq!(cell0.sub_entity(0, 3), 3);
        assert_eq!(cell0.sub_entity(1, 0), 0);
        assert_eq!(cell0.sub_entity(1, 1), 1);
        assert_eq!(cell0.sub_entity(1, 2), 2);
        assert_eq!(cell0.sub_entity(1, 3), 3);
        assert_eq!(cell0.sub_entity(1, 4), 4);
        assert_eq!(cell0.sub_entity(1, 5), 5);
        assert_eq!(cell0.sub_entity(2, 0), 0);
        assert_eq!(cell0.sub_entity(2, 1), 1);
        assert_eq!(cell0.sub_entity(2, 2), 2);
        assert_eq!(cell0.sub_entity(2, 3), 3);
        assert_eq!(cell0.sub_entity(3, 0), 0);
        let cell1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Tetrahedron, 1);
        assert_eq!(cell1.sub_entity(0, 0), 4);
        assert_eq!(cell1.sub_entity(0, 1), 0);
        assert_eq!(cell1.sub_entity(0, 2), 2);
        assert_eq!(cell1.sub_entity(0, 3), 3);
        assert_eq!(cell1.sub_entity(1, 0), 0);
        assert_eq!(cell1.sub_entity(1, 1), 3);
        assert_eq!(cell1.sub_entity(1, 2), 4);
        assert_eq!(cell1.sub_entity(1, 3), 6);
        assert_eq!(cell1.sub_entity(1, 4), 7);
        assert_eq!(cell1.sub_entity(1, 5), 8);
        assert_eq!(cell1.sub_entity(2, 0), 1);
        assert_eq!(cell1.sub_entity(2, 1), 4);
        assert_eq!(cell1.sub_entity(2, 2), 5);
        assert_eq!(cell1.sub_entity(2, 3), 6);
        assert_eq!(cell1.sub_entity(3, 0), 1);

        let face0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 0);
        assert_eq!(face0.sub_entity(0, 0), 1);
        assert_eq!(face0.sub_entity(0, 1), 2);
        assert_eq!(face0.sub_entity(0, 2), 3);
        assert_eq!(face0.sub_entity(1, 0), 0);
        assert_eq!(face0.sub_entity(1, 1), 1);
        assert_eq!(face0.sub_entity(1, 2), 2);
        assert_eq!(face0.sub_entity(2, 0), 0);
        let face1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 1);
        assert_eq!(face1.sub_entity(0, 0), 0);
        assert_eq!(face1.sub_entity(0, 1), 2);
        assert_eq!(face1.sub_entity(0, 2), 3);
        assert_eq!(face1.sub_entity(1, 0), 0);
        assert_eq!(face1.sub_entity(1, 1), 3);
        assert_eq!(face1.sub_entity(1, 2), 4);
        assert_eq!(face1.sub_entity(2, 0), 1);
        let face2 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 2);
        assert_eq!(face2.sub_entity(0, 0), 0);
        assert_eq!(face2.sub_entity(0, 1), 1);
        assert_eq!(face2.sub_entity(0, 2), 3);
        assert_eq!(face2.sub_entity(1, 0), 1);
        assert_eq!(face2.sub_entity(1, 1), 3);
        assert_eq!(face2.sub_entity(1, 2), 5);
        assert_eq!(face2.sub_entity(2, 0), 2);
        let face3 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 3);
        assert_eq!(face3.sub_entity(0, 0), 0);
        assert_eq!(face3.sub_entity(0, 1), 1);
        assert_eq!(face3.sub_entity(0, 2), 2);
        assert_eq!(face3.sub_entity(1, 0), 2);
        assert_eq!(face3.sub_entity(1, 1), 4);
        assert_eq!(face3.sub_entity(1, 2), 5);
        assert_eq!(face3.sub_entity(2, 0), 3);
        let face4 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 4);
        assert_eq!(face4.sub_entity(0, 0), 2);
        assert_eq!(face4.sub_entity(0, 1), 3);
        assert_eq!(face4.sub_entity(0, 2), 4);
        assert_eq!(face4.sub_entity(1, 0), 0);
        assert_eq!(face4.sub_entity(1, 1), 6);
        assert_eq!(face4.sub_entity(1, 2), 7);
        assert_eq!(face4.sub_entity(2, 0), 4);
        let face5 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 5);
        assert_eq!(face5.sub_entity(0, 0), 0);
        assert_eq!(face5.sub_entity(0, 1), 3);
        assert_eq!(face5.sub_entity(0, 2), 4);
        assert_eq!(face5.sub_entity(1, 0), 3);
        assert_eq!(face5.sub_entity(1, 1), 6);
        assert_eq!(face5.sub_entity(1, 2), 8);
        assert_eq!(face5.sub_entity(2, 0), 5);
        let face6 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Triangle, 6);
        assert_eq!(face6.sub_entity(0, 0), 0);
        assert_eq!(face6.sub_entity(0, 1), 2);
        assert_eq!(face6.sub_entity(0, 2), 4);
        assert_eq!(face6.sub_entity(1, 0), 4);
        assert_eq!(face6.sub_entity(1, 1), 7);
        assert_eq!(face6.sub_entity(1, 2), 8);
        assert_eq!(face6.sub_entity(2, 0), 6);

        let edge0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 0);
        assert_eq!(edge0.sub_entity(0, 0), 2);
        assert_eq!(edge0.sub_entity(0, 1), 3);
        assert_eq!(edge0.sub_entity(1, 0), 0);
        let edge1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 1);
        assert_eq!(edge1.sub_entity(0, 0), 1);
        assert_eq!(edge1.sub_entity(0, 1), 3);
        assert_eq!(edge1.sub_entity(1, 0), 1);
        let edge2 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 2);
        assert_eq!(edge2.sub_entity(0, 0), 1);
        assert_eq!(edge2.sub_entity(0, 1), 2);
        assert_eq!(edge2.sub_entity(1, 0), 2);
        let edge3 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 3);
        assert_eq!(edge3.sub_entity(0, 0), 0);
        assert_eq!(edge3.sub_entity(0, 1), 3);
        assert_eq!(edge3.sub_entity(1, 0), 3);
        let edge4 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 4);
        assert_eq!(edge4.sub_entity(0, 0), 0);
        assert_eq!(edge4.sub_entity(0, 1), 2);
        assert_eq!(edge4.sub_entity(1, 0), 4);
        let edge5 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 5);
        assert_eq!(edge5.sub_entity(0, 0), 0);
        assert_eq!(edge5.sub_entity(0, 1), 1);
        assert_eq!(edge5.sub_entity(1, 0), 5);
        let edge6 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 6);
        assert_eq!(edge6.sub_entity(0, 0), 3);
        assert_eq!(edge6.sub_entity(0, 1), 4);
        assert_eq!(edge6.sub_entity(1, 0), 6);
        let edge7 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 7);
        assert_eq!(edge7.sub_entity(0, 0), 2);
        assert_eq!(edge7.sub_entity(0, 1), 4);
        assert_eq!(edge7.sub_entity(1, 0), 7);
        let edge8 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Interval, 8);
        assert_eq!(edge8.sub_entity(0, 0), 0);
        assert_eq!(edge8.sub_entity(0, 1), 4);
        assert_eq!(edge8.sub_entity(1, 0), 8);

        let vertex0 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 0);
        assert_eq!(vertex0.sub_entity(0, 0), 0);
        let vertex1 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 1);
        assert_eq!(vertex1.sub_entity(0, 0), 1);
        let vertex2 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 2);
        assert_eq!(vertex2.sub_entity(0, 0), 2);
        let vertex3 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 3);
        assert_eq!(vertex3.sub_entity(0, 0), 3);
        let vertex4 = SingleTypeEntityTopology::new(&t, ReferenceCellType::Point, 4);
        assert_eq!(vertex4.sub_entity(0, 0), 4);
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
                                    for value in col.iter() {
                                        assert!(uc_ji[value].contains(&c));
                                    }
                                }
                                for (k, uc_jik) in uc_ji.iter().enumerate() {
                                    for value in uc_jik {
                                        assert!(dc_ij.r().slice(1, *value).data().contains(&k));
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
                                for (i, j) in cell.sub_entity_iter(dim).enumerate() {
                                    assert_eq!(j, cell.sub_entity(dim, i));
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
}
