//! Topology for grids where entities of each tdim may be a mixture of types

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
pub struct MixedTopology {
    dim: usize,
    pub(crate) ids: HashMap<ReferenceCellType, Vec<usize>>,
    pub(crate) ids_to_indices: HashMap<ReferenceCellType, HashMap<usize, usize>>,
    pub(crate) insertion_indices: HashMap<ReferenceCellType, Vec<usize>>,
    entity_types: Vec<Vec<ReferenceCellType>>,
    entity_counts: HashMap<ReferenceCellType, usize>,
    pub(crate) downward_connectivity:
        HashMap<ReferenceCellType, HashMap<ReferenceCellType, DynArray<usize, 2>>>,
    pub(crate) upward_connectivity:
        HashMap<ReferenceCellType, HashMap<ReferenceCellType, Vec<Vec<usize>>>>,
    pub(crate) orientation: HashMap<ReferenceCellType, Vec<i32>>,
}

#[cfg(feature = "serde")]
#[derive(serde::Serialize, Debug, serde::Deserialize)]
/// Serde serializable topology
pub struct SerializableTopology {
    dim: usize,
    ids: HashMap<ReferenceCellType, Vec<usize>>,
    ids_to_indices: HashMap<ReferenceCellType, HashMap<usize, usize>>,
    insertion_indices: HashMap<ReferenceCellType, Vec<usize>>,
    entity_types: Vec<Vec<ReferenceCellType>>,
    entity_counts: HashMap<ReferenceCellType, usize>,
    #[allow(clippy::type_complexity)]
    downward_connectivity:
        HashMap<ReferenceCellType, HashMap<ReferenceCellType, (Vec<usize>, [usize; 2])>>,
    upward_connectivity: HashMap<ReferenceCellType, HashMap<ReferenceCellType, Vec<Vec<usize>>>>,
    orientation: HashMap<ReferenceCellType, Vec<i32>>,
}

#[cfg(feature = "serde")]
impl ConvertToSerializable for MixedTopology {
    type SerializableType = SerializableTopology;
    fn to_serializable(&self) -> SerializableTopology {
        SerializableTopology {
            dim: self.dim,
            ids: self.ids.clone(),
            ids_to_indices: self.ids_to_indices.clone(),
            insertion_indices: self.insertion_indices.clone(),
            entity_types: self.entity_types.clone(),
            entity_counts: self.entity_counts.clone(),
            downward_connectivity: self
                .downward_connectivity
                .iter()
                .map(|(a, b)| {
                    (
                        *a,
                        b.iter()
                            .map(|(c, d)| (*c, (d.data().unwrap().to_vec(), d.shape())))
                            .collect::<HashMap<_, _>>(),
                    )
                })
                .collect::<HashMap<_, _>>(),
            upward_connectivity: self.upward_connectivity.clone(),
            orientation: self.orientation.clone(),
        }
    }
    fn from_serializable(s: SerializableTopology) -> Self {
        Self {
            dim: s.dim,
            ids: s.ids,
            ids_to_indices: s.ids_to_indices,
            insertion_indices: s.insertion_indices,
            entity_types: s.entity_types,
            entity_counts: s.entity_counts,
            downward_connectivity: s
                .downward_connectivity
                .iter()
                .map(|(a, b)| {
                    (
                        *a,
                        b.iter()
                            .map(|(c, (data, shape))| {
                                (*c, {
                                    let mut d = DynArray::<usize, 2>::from_shape(*shape);
                                    d.data_mut().unwrap().copy_from_slice(data);
                                    d
                                })
                            })
                            .collect::<HashMap<_, _>>(),
                    )
                })
                .collect::<HashMap<_, _>>(),
            upward_connectivity: s.upward_connectivity,
            orientation: s.orientation,
        }
    }
}

impl MixedTopology {
    /// Create a topology
    pub fn new(
        cells: &[usize],
        cell_types: &[ReferenceCellType],
        vertex_ids: Option<Vec<usize>>,
        cell_ids: Option<Vec<usize>>,
    ) -> Self {
        let mut start = 0;

        let dim = reference_cell::dim(cell_types[0]);
        // Check that all cells have the same topological dim
        for ct in cell_types {
            if reference_cell::dim(*ct) != dim {
                panic!(
                    "MixedTopology does not support cells of a mixture of topological dimensions"
                );
            }
        }

        let mut entities = HashMap::new();
        let mut entity_counts = HashMap::new();
        entity_counts.insert(ReferenceCellType::Point, cells.iter().max().unwrap() + 1);

        let mut ids = HashMap::new();
        let mut ids_to_indices = HashMap::new();
        if let Some(i) = vertex_ids {
            ids_to_indices.insert(ReferenceCellType::Point, {
                let mut vids = HashMap::new();
                for (j, k) in i.iter().enumerate() {
                    vids.insert(*k, j);
                }
                vids
            });
            ids.insert(ReferenceCellType::Point, i);
        }

        let mut insertion_indices = HashMap::new();
        for (cell_index, cell_type) in cell_types.iter().enumerate() {
            let size = reference_cell::entity_counts(*cell_type)[0];
            insertion_indices
                .entry(*cell_type)
                .or_insert(vec![])
                .push(cell_index);

            // Get vertices of cell
            let cell = &cells[start..start + size];
            if let Some(i) = &cell_ids {
                ids.entry(*cell_type).or_insert(vec![]).push(i[cell_index]);
                ids_to_indices
                    .entry(*cell_type)
                    .or_insert(HashMap::new())
                    .insert(
                        i[cell_index],
                        if let Some(j) = entity_counts.get(cell_type) {
                            *j
                        } else {
                            0
                        },
                    );
            }

            // Add each subentity of cell that's not already included in entities
            for (etypes, rc_i) in izip!(
                &reference_cell::entity_types(*cell_type),
                &reference_cell::connectivity(*cell_type)
            ) {
                for (etype, c_ij) in izip!(etypes, rc_i) {
                    let mut entity = c_ij[0].iter().map(|i| cell[*i]).collect::<Vec<_>>();
                    entity.sort();
                    entities
                        .entry(*etype)
                        .or_insert(HashMap::new())
                        .entry(entity)
                        .or_insert_with(|| {
                            let old = *entity_counts.entry(*etype).or_insert(0);
                            *entity_counts.get_mut(etype).unwrap() += 1;
                            old
                        });
                }
            }
            start += size;
        }

        let mut entity_types = vec![vec![]; dim + 1];
        for e in entity_counts.keys() {
            entity_types[reference_cell::dim(*e)].push(*e);
        }

        // Downward connectivity: The entities of type t1 that are subentities of
        // entities of type t0 (with dim(t0)>=dim(t1))
        // downward_connectivity[t0][t1][[.., t0_entity_index]] = dim1_entity_index
        let mut downward_connectivity = HashMap::new();
        // Assign memory for downward connectivity
        for (etype, ecount) in &entity_counts {
            downward_connectivity.insert(*etype, {
                let mut dc = HashMap::new();
                for sub_etypes in &reference_cell::entity_types(*etype) {
                    for e in sub_etypes {
                        dc.entry(*e).or_insert(rlst_dynamic_array!(
                            usize,
                            [sub_etypes.iter().filter(|&i| *i == *e).count(), *ecount]
                        ));
                    }
                }
                dc
            });
        }

        // Upward connectivity: The entities of type t1 that are superentities of
        // entities of type t0 (with dim(t0)<dim(t1))
        // upward_connectivity[t0][t1][t0_entity_index][..] = t1_entity_index
        let mut upward_connectivity = HashMap::new();
        // Add correct nu ber of items to upward connectivity
        for etype in entity_counts.keys() {
            for sub_etypes in reference_cell::entity_types(*etype)
                .iter()
                .take(reference_cell::dim(*etype))
            {
                for e in sub_etypes {
                    upward_connectivity
                        .entry(*e)
                        .or_insert(HashMap::new())
                        .entry(*etype)
                        .or_insert(vec![vec![]; entity_counts[e]]);
                }
            }
        }

        // downward_connectivity[t][t][i] = [i] (ie each entity is a sub-entity of itself)
        for (d, dc) in downward_connectivity.iter_mut() {
            for (i, mut j) in dc.get_mut(d).unwrap().col_iter_mut().enumerate() {
                j[[0]] = i;
            }
        }

        // downward_connectivity[t][Point][i] = vertices of entity i of type t
        for (etype, entities_vertices) in &entities {
            for (vs, entity_index) in entities_vertices {
                for (i, vertex_index) in vs.iter().enumerate() {
                    downward_connectivity
                        .get_mut(etype)
                        .unwrap()
                        .get_mut(&ReferenceCellType::Point)
                        .unwrap()[[i, *entity_index]] = *vertex_index;
                    // We can also fill out the upward connectivity if dim > 0.
                    if reference_cell::dim(*etype) > 0 {
                        let uc = upward_connectivity
                            .get_mut(&ReferenceCellType::Point)
                            .unwrap()
                            .get_mut(etype)
                            .unwrap();
                        if !uc[*vertex_index].contains(entity_index) {
                            uc[*vertex_index].push(*entity_index);
                        }
                    }
                }
            }
        }

        // Fill remaining connectivity
        if dim > 0 {
            for etype in &entity_types[dim] {
                let entity_types = reference_cell::entity_types(*etype);
                let ref_conn = reference_cell::connectivity(*etype);
                // Cell entities, skipping vertices
                let mut cell_entities = ref_conn
                    .iter()
                    .skip(1)
                    .map(|i| vec![0; i.len()])
                    .collect::<Vec<_>>();
                for (cell, cell_index) in &entities[etype] {
                    // dim - 1 denotes the dim dimensional subentity as we skipped vertices
                    cell_entities[dim - 1][0] = *cell_index;

                    // Fill in subentity indices
                    for (ce_i, rc_i, et_i) in izip!(
                        cell_entities.iter_mut(),
                        ref_conn.iter().skip(1),
                        entity_types.iter().skip(1),
                    ) {
                        for (ce_ij, rc_ij, et_ij) in izip!(ce_i.iter_mut(), rc_i, et_i) {
                            let mut entity = rc_ij[0].iter().map(|i| cell[*i]).collect::<Vec<_>>();
                            entity.sort();
                            *ce_ij = *entities[et_ij].get(&entity).unwrap();
                        }
                    }

                    for (i, (rc_i, ce_i, et_i)) in izip!(
                        ref_conn.iter().skip(2),
                        cell_entities.iter().skip(1),
                        entity_types.iter().skip(2),
                    )
                    .enumerate()
                    {
                        for (ce_ij, rc_ij, et_ij) in izip!(ce_i, rc_i, et_i) {
                            for (k, (rc_ijk, ce_k)) in
                                izip!(rc_ij.iter().take(i + 2).skip(1), &cell_entities).enumerate()
                            {
                                let mut sub_entity_indices = HashMap::new();
                                for (sub_entity_type, rc_ijkl) in
                                    izip!(&entity_types[k + 1], rc_ijk)
                                {
                                    let sub_entity_index = sub_entity_indices
                                        .entry(sub_entity_type)
                                        .and_modify(|e| *e += 1)
                                        .or_insert(0);
                                    downward_connectivity
                                        .get_mut(et_ij)
                                        .unwrap()
                                        .get_mut(sub_entity_type)
                                        .unwrap()[[*sub_entity_index, *ce_ij]] = ce_k[*rc_ijkl];
                                    let uc = upward_connectivity
                                        .get_mut(sub_entity_type)
                                        .unwrap()
                                        .get_mut(et_ij)
                                        .unwrap();
                                    if !uc[ce_k[*rc_ijkl]].contains(ce_ij) {
                                        uc[ce_k[*rc_ijkl]].push(*ce_ij);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        let mut orientation = HashMap::new();
        for types in entity_types.iter().skip(1) {
            for t in types {
                let dc = &downward_connectivity[t][&ReferenceCellType::Point];
                orientation.insert(
                    *t,
                    (0..dc.shape()[1])
                        .map(|i| {
                            compute_orientation(
                                *t,
                                &dc.data().unwrap()[i * dc.shape()[0]..(i + 1) * dc.shape()[0]],
                            )
                        })
                        .collect::<Vec<_>>(),
                );
            }
        }

        Self {
            dim,
            ids,
            ids_to_indices,
            insertion_indices,
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
    pub fn entity_types(&self) -> &[Vec<ReferenceCellType>] {
        &self.entity_types
    }
    /// Entity counts
    pub fn entity_counts(&self) -> &HashMap<ReferenceCellType, usize> {
        &self.entity_counts
    }
    /// Entity counts
    pub fn entity_count(&self, entity_type: ReferenceCellType) -> usize {
        if let Some(n) = self.entity_counts.get(&entity_type) {
            *n
        } else {
            0
        }
    }
    /// Cell sub-entity index
    pub fn cell_entity_index(
        &self,
        cell_type: ReferenceCellType,
        cell_index: usize,
        entity_type: ReferenceCellType,
        entity_index: usize,
    ) -> usize {
        self.downward_connectivity[&cell_type][&entity_type][[entity_index, cell_index]]
    }
    /// Entity id
    pub fn entity_id(&self, entity_type: ReferenceCellType, entity_index: usize) -> Option<usize> {
        self.ids.get(&entity_type).map(|a| a[entity_index])
    }
}

/// Topology of a cell
#[derive(Debug)]
pub struct MixedEntityTopology<'a> {
    topology: &'a MixedTopology,
    entity_type: ReferenceCellType,
    entity_index: usize,
}

impl<'t> MixedEntityTopology<'t> {
    /// Create new
    pub fn new(
        topology: &'t MixedTopology,
        entity_type: ReferenceCellType,
        entity_index: usize,
    ) -> Self {
        Self {
            topology,
            entity_type,
            entity_index,
        }
    }
}

impl Topology for MixedEntityTopology<'_> {
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
        self.topology.upward_connectivity[&self.entity_type][&entity_type][self.entity_index]
            .iter()
            .copied()
    }

    fn sub_entity_iter(
        &self,
        entity_type: ReferenceCellType,
    ) -> Copied<std::slice::Iter<'_, usize>> {
        let rows = self.topology.downward_connectivity[&self.entity_type][&entity_type].shape()[0];
        self.topology.downward_connectivity[&self.entity_type][&entity_type]
            .data()
            .unwrap()[rows * self.entity_index..rows * (self.entity_index + 1)]
            .iter()
            .copied()
    }

    fn sub_entity(&self, entity_type: ReferenceCellType, index: usize) -> usize {
        self.topology.downward_connectivity[&self.entity_type][&entity_type]
            [[index, self.entity_index]]
    }

    fn orientation(&self) -> i32 {
        if reference_cell::dim(self.entity_type) == 0 {
            0
        } else {
            self.topology.orientation[&self.entity_type][self.entity_index]
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn example_topology_triangle_and_quad() -> MixedTopology {
        //! An example topology
        MixedTopology::new(
            &[0, 1, 2, 3, 1, 3, 4],
            &[
                ReferenceCellType::Quadrilateral,
                ReferenceCellType::Triangle,
            ],
            None,
            None,
        )
    }

    #[test]
    fn test_prism() {
        let _ = MixedTopology::new(&[0, 1, 2, 3, 4, 5], &[ReferenceCellType::Prism], None, None);
    }
    #[test]
    fn test_pyramid() {
        let _ = MixedTopology::new(&[0, 1, 2, 3, 4], &[ReferenceCellType::Pyramid], None, None);
    }

    #[test]
    fn test_sub_entities_triangle_and_quad() {
        //! Test sub-entities of a triangle
        let t = example_topology_triangle_and_quad();
        let cell0 = MixedEntityTopology::new(&t, ReferenceCellType::Quadrilateral, 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 1), 1);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 2), 2);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Point, 3), 3);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 0), 0);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 1), 1);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 2), 2);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Interval, 3), 3);
        assert_eq!(cell0.sub_entity(ReferenceCellType::Quadrilateral, 0), 0);
        let cell1 = MixedEntityTopology::new(&t, ReferenceCellType::Triangle, 0);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 0), 1);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 1), 3);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Point, 2), 4);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 0), 4);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 1), 5);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Interval, 2), 2);
        assert_eq!(cell1.sub_entity(ReferenceCellType::Triangle, 0), 0);

        let edge0 = MixedEntityTopology::new(&t, ReferenceCellType::Interval, 0);
        assert_eq!(edge0.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(edge0.sub_entity(ReferenceCellType::Point, 1), 1);
        assert_eq!(edge0.sub_entity(ReferenceCellType::Interval, 0), 0);
        let edge1 = MixedEntityTopology::new(&t, ReferenceCellType::Interval, 1);
        assert_eq!(edge1.sub_entity(ReferenceCellType::Point, 0), 0);
        assert_eq!(edge1.sub_entity(ReferenceCellType::Point, 1), 2);
        assert_eq!(edge1.sub_entity(ReferenceCellType::Interval, 0), 1);
        let edge2 = MixedEntityTopology::new(&t, ReferenceCellType::Interval, 2);
        assert_eq!(edge2.sub_entity(ReferenceCellType::Point, 0), 1);
        assert_eq!(edge2.sub_entity(ReferenceCellType::Point, 1), 3);
        assert_eq!(edge2.sub_entity(ReferenceCellType::Interval, 0), 2);
        let edge3 = MixedEntityTopology::new(&t, ReferenceCellType::Interval, 3);
        assert_eq!(edge3.sub_entity(ReferenceCellType::Point, 0), 2);
        assert_eq!(edge3.sub_entity(ReferenceCellType::Point, 1), 3);
        assert_eq!(edge3.sub_entity(ReferenceCellType::Interval, 0), 3);
        let edge4 = MixedEntityTopology::new(&t, ReferenceCellType::Interval, 4);
        assert_eq!(edge4.sub_entity(ReferenceCellType::Point, 0), 3);
        assert_eq!(edge4.sub_entity(ReferenceCellType::Point, 1), 4);
        assert_eq!(edge4.sub_entity(ReferenceCellType::Interval, 0), 4);
        let edge5 = MixedEntityTopology::new(&t, ReferenceCellType::Interval, 5);
        assert_eq!(edge5.sub_entity(ReferenceCellType::Point, 0), 1);
        assert_eq!(edge5.sub_entity(ReferenceCellType::Point, 1), 4);
        assert_eq!(edge5.sub_entity(ReferenceCellType::Interval, 0), 5);

        let vertex0 = MixedEntityTopology::new(&t, ReferenceCellType::Point, 0);
        assert_eq!(vertex0.sub_entity(ReferenceCellType::Point, 0), 0);
        let vertex1 = MixedEntityTopology::new(&t, ReferenceCellType::Point, 1);
        assert_eq!(vertex1.sub_entity(ReferenceCellType::Point, 0), 1);
        let vertex2 = MixedEntityTopology::new(&t, ReferenceCellType::Point, 2);
        assert_eq!(vertex2.sub_entity(ReferenceCellType::Point, 0), 2);
        let vertex3 = MixedEntityTopology::new(&t, ReferenceCellType::Point, 3);
        assert_eq!(vertex3.sub_entity(ReferenceCellType::Point, 0), 3);
        let vertex4 = MixedEntityTopology::new(&t, ReferenceCellType::Point, 4);
        assert_eq!(vertex4.sub_entity(ReferenceCellType::Point, 0), 4);
    }
}
