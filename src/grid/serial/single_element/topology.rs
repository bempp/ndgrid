//! Implementation of grid topology

use crate::traits::Topology;
use crate::types::Array2D;
use itertools::izip;
use ndelement::reference_cell;
use ndelement::types::ReferenceCellType;
use rlst::{rlst_dynamic_array2, DefaultIteratorMut, RawAccess, Shape};
// use std::collections::HashMap;

/// Topology of a single element grid
pub struct SingleElementTopology {
    dim: usize,
    //    index_map: Vec<usize>,
    //entity_types: Vec<ReferenceCellType>,
    downward_connectivity: Vec<Vec<Array2D<usize>>>,
    //upward_connectivity: Vec<Vec<IntegerArray2>>,
    //    entities_to_vertices: Vec<Vec<Vec<usize>>>,
    //    cells_to_entities: Vec<Vec<Vec<usize>>>,
    //    entities_to_cells: Vec<Vec<Vec<CellLocalIndexPair<usize>>>>,
    //    vertex_indices_to_ids: Vec<usize>,
    //    vertex_ids_to_indices: HashMap<usize, usize>,
    //    edge_indices_to_ids: Vec<usize>,
    //    edge_ids_to_indices: HashMap<usize, usize>,
    //    cell_indices_to_ids: Vec<usize>,
    //    cell_ids_to_indices: HashMap<usize, usize>,
    //    cell_types: [ReferenceCellType; 1],
}

unsafe impl Sync for SingleElementTopology {}

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
        _ => {
            unimplemented!();
        }
    }
}

impl SingleElementTopology {
    /// Create a topology
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        cells: &[usize],
        cell_type: ReferenceCellType,
        _point_ids: Option<&[usize]>,
        _cell_ids: Option<&[usize]>,
    ) -> Self {
        // Cells where faces are mixture of triangles and quads not supported
        if cell_type != ReferenceCellType::Point
            && cell_type != ReferenceCellType::Triangle
            && cell_type != ReferenceCellType::Quadrilateral
            && cell_type != ReferenceCellType::Tetrahedron
            && cell_type != ReferenceCellType::Hexahedron
        {
            panic!("Unsupported cell type for SingleElementTopology: {cell_type:?}");
        }
        let size = reference_cell::entity_counts(cell_type)[0];
        let ncells = cells.len() / size;
        let dim = reference_cell::dim(cell_type);
        let ref_conn = reference_cell::connectivity(cell_type);
        let etypes = reference_cell::entity_types(cell_type);

        // List of entities by dimension
        let mut entities = vec![vec![]; dim - 1];
        for cell_index in 0..ncells {
            let cell = &cells[cell_index * size..(cell_index + 1) * size];
            // Iterate over topological dimension
            for (e_i, rc_i, et_i) in izip!(
                entities.iter_mut(),
                ref_conn.iter().take(dim).skip(1),
                etypes.iter().take(dim).skip(1)
            ) {
                // For each entity of the given dimension
                for (c_ij, et_ij) in izip!(rc_i, et_i) {
                    let mut entity = c_ij[0].iter().map(|i| cell[*i]).collect::<Vec<_>>();
                    orient_entity(*et_ij, &mut entity);
                    if !e_i.contains(&entity) {
                        e_i.push(entity);
                    }
                }
            }
        }
        // Number of entities by dimension
        let mut nentities = vec![0; dim + 1];
        nentities[0] = cells.iter().max().unwrap() + 1;
        nentities[dim] = ncells;
        for d in 1..dim {
            nentities[d] = entities[d - 1].len();
        }

        // Downward connectivity: The entities of dimension dim1 that are subentities of entities of dimension dim0 (with dim0>dim1) (eg edges of a triangle, vertices of a tetrahedron, etc)
        // indices: downward_connectivity[dim0][dim1][[dim1_entity_index, dim0_entity_index]]
        let mut downward_connectivity = nentities
            .iter()
            .enumerate()
            .map(|(i, j)| {
                ref_conn[i][0]
                    .iter()
                    .take(i + 1)
                    .map(|r| rlst_dynamic_array2!(usize, [r.len(), *j]))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        // downward_connectivity[d][d][i] = [i] (ie each entity is a sub-entity of itself)
        for (d, dc) in downward_connectivity.iter_mut().enumerate() {
            for (i, mut j) in dc[d].col_iter_mut().enumerate() {
                j[[0]] = i;
            }
        }

        // downward_connectivity[dim][0] = vertices of each cell
        for (i, mut j) in downward_connectivity[dim][0].col_iter_mut().enumerate() {
            for (k, l) in j.iter_mut().zip(&cells[i * size..(i + 1) * size]) {
                *k = *l;
            }
        }
        // downward_connectivity[i][0] = vertices of entity
        for (es, dc) in entities
            .iter()
            .zip(downward_connectivity.iter_mut().take(dim).skip(1))
        {
            for (e, mut c) in es.iter().zip(dc[0].col_iter_mut()) {
                for (i, j) in e.iter().zip(c.iter_mut()) {
                    *j = *i;
                }
            }
        }

        let mut cell_entities = ref_conn
            .iter()
            .skip(1)
            .map(|i| vec![0; i.len()])
            .collect::<Vec<_>>();
        for cell_index in 0..ncells {
            // Collect indices of each subentity of the cell
            cell_entities[dim - 1][0] = cell_index;
            let cell = &cells[cell_index * size..(cell_index + 1) * size];
            for (e_i, ce_i, rc_i, et_i) in izip!(
                entities.iter(),
                cell_entities.iter_mut(),
                ref_conn.iter().skip(1),
                etypes.iter().skip(1)
            ) {
                for (ce_ij, rc_ij, et_ij) in izip!(ce_i.iter_mut(), rc_i, et_i) {
                    let mut entity = rc_ij[0].iter().map(|i| cell[*i]).collect::<Vec<_>>();
                    orient_entity(*et_ij, &mut entity);
                    *ce_ij = e_i.iter().position(|r| *r == entity).unwrap();
                }
            }
            // Copy these indices into connectivity for each dim
            // Loop over dim0
            for (i, (dc_i, rc_i, ce_i)) in izip!(
                downward_connectivity.iter_mut().skip(2),
                ref_conn.iter().skip(2),
                cell_entities.iter().skip(1)
            )
            .enumerate()
            {
                // Loop over entities of dimension dim0
                for (ce_ij, rc_ij) in izip!(ce_i, rc_i) {
                    // Loop over dim1
                    for (dc_ik, rc_ijk, ce_k) in izip!(
                        dc_i.iter_mut().take(i + 2).skip(1),
                        rc_ij.iter().take(i + 2).skip(1),
                        &cell_entities
                    ) {
                        // Loop over entities of dimension dim1
                        for (l, rc_ijkl) in rc_ijk.iter().enumerate() {
                            dc_ik[[l, *ce_ij]] = ce_k[*rc_ijkl];
                        }
                    }
                }
            }
        }

        Self {
            dim,
            downward_connectivity,
        }
    }
    /// Topological dimension
    pub fn dim(&self) -> usize {
        self.dim
    }
}

pub struct IndexIter {}
impl std::iter::Iterator for IndexIter {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        None
    }
}

use std::iter::Copied;

/// Topology of a cell
pub struct SingleElementCellTopology<'a> {
    topology: &'a SingleElementTopology,
    //entity_type: ReferenceCellType,
    entity_index: usize,
    dim: usize,
}

impl<'t> SingleElementCellTopology<'t> {
    /// Create new
    pub fn new(
        topology: &'t SingleElementTopology,
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
impl<'t> Topology for SingleElementCellTopology<'t> {
    type EntityIndexIter<'a> = Copied<std::slice::Iter<'a, usize>>
    where
        Self: 'a;

    type ConnectedEntityIndexIter<'a> = IndexIter
    where
        Self: 'a;

    fn connected_entity_iter(&self, _dim: usize) -> IndexIter {
        unimplemented!();
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

    fn example_topology_triangle() -> SingleElementTopology {
        //! An example topology
        SingleElementTopology::new(&[0, 1, 2, 2, 1, 3], ReferenceCellType::Triangle, None, None)
    }

    fn example_topology_tetrahedron() -> SingleElementTopology {
        //! An example topology
        SingleElementTopology::new(
            &[0, 1, 2, 3, 4, 0, 2, 3],
            ReferenceCellType::Tetrahedron,
            None,
            None,
        )
    }

    #[test]
    fn test_sub_entities_triangle() {
        //! Test sub-entities of a triangle
        let t = example_topology_triangle();
        let cell0 = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, 0);
        assert_eq!(cell0.sub_entity(0, 0), 0);
        assert_eq!(cell0.sub_entity(0, 1), 1);
        assert_eq!(cell0.sub_entity(0, 2), 2);
        assert_eq!(cell0.sub_entity(1, 0), 0);
        assert_eq!(cell0.sub_entity(1, 1), 1);
        assert_eq!(cell0.sub_entity(1, 2), 2);
        assert_eq!(cell0.sub_entity(2, 0), 0);
        let cell1 = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, 1);
        assert_eq!(cell1.sub_entity(0, 0), 2);
        assert_eq!(cell1.sub_entity(0, 1), 1);
        assert_eq!(cell1.sub_entity(0, 2), 3);
        assert_eq!(cell1.sub_entity(1, 0), 3);
        assert_eq!(cell1.sub_entity(1, 1), 4);
        assert_eq!(cell1.sub_entity(1, 2), 0);
        assert_eq!(cell1.sub_entity(2, 0), 1);

        let edge0 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 0);
        assert_eq!(edge0.sub_entity(0, 0), 1);
        assert_eq!(edge0.sub_entity(0, 1), 2);
        assert_eq!(edge0.sub_entity(1, 0), 0);
        let edge1 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 1);
        assert_eq!(edge1.sub_entity(0, 0), 0);
        assert_eq!(edge1.sub_entity(0, 1), 2);
        assert_eq!(edge1.sub_entity(1, 0), 1);
        let edge2 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 2);
        assert_eq!(edge2.sub_entity(0, 0), 0);
        assert_eq!(edge2.sub_entity(0, 1), 1);
        assert_eq!(edge2.sub_entity(1, 0), 2);
        let edge3 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 3);
        assert_eq!(edge3.sub_entity(0, 0), 1);
        assert_eq!(edge3.sub_entity(0, 1), 3);
        assert_eq!(edge3.sub_entity(1, 0), 3);
        let edge4 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 4);
        assert_eq!(edge4.sub_entity(0, 0), 2);
        assert_eq!(edge4.sub_entity(0, 1), 3);
        assert_eq!(edge4.sub_entity(1, 0), 4);

        let vertex0 = SingleElementCellTopology::new(&t, ReferenceCellType::Point, 0);
        assert_eq!(vertex0.sub_entity(0, 0), 0);
        let vertex1 = SingleElementCellTopology::new(&t, ReferenceCellType::Point, 1);
        assert_eq!(vertex1.sub_entity(0, 0), 1);
        let vertex2 = SingleElementCellTopology::new(&t, ReferenceCellType::Point, 2);
        assert_eq!(vertex2.sub_entity(0, 0), 2);
        let vertex3 = SingleElementCellTopology::new(&t, ReferenceCellType::Point, 3);
        assert_eq!(vertex3.sub_entity(0, 0), 3);
    }

    #[test]
    fn test_sub_entities_tetrahedron() {
        //! Test sub-entities of a tetrahedron
        let t = example_topology_tetrahedron();
        let cell0 = SingleElementCellTopology::new(&t, ReferenceCellType::Tetrahedron, 0);
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
        let cell1 = SingleElementCellTopology::new(&t, ReferenceCellType::Tetrahedron, 1);
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

        let face0 = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, 0);
        assert_eq!(face0.sub_entity(0, 0), 1);
        assert_eq!(face0.sub_entity(0, 1), 2);
        assert_eq!(face0.sub_entity(0, 2), 3);
        assert_eq!(face0.sub_entity(1, 0), 0);
        assert_eq!(face0.sub_entity(1, 1), 1);
        assert_eq!(face0.sub_entity(1, 2), 2);
        assert_eq!(face0.sub_entity(2, 0), 0);
        let face1 = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, 1);
        assert_eq!(face1.sub_entity(0, 0), 0);
        assert_eq!(face1.sub_entity(0, 1), 2);
        assert_eq!(face1.sub_entity(0, 2), 3);
        assert_eq!(face1.sub_entity(1, 0), 0);
        assert_eq!(face1.sub_entity(1, 1), 3);
        assert_eq!(face1.sub_entity(1, 2), 4);
        assert_eq!(face1.sub_entity(2, 0), 1);
        let face2 = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, 2);
        assert_eq!(face2.sub_entity(0, 0), 0);
        assert_eq!(face2.sub_entity(0, 1), 1);
        assert_eq!(face2.sub_entity(0, 2), 3);
        assert_eq!(face2.sub_entity(1, 0), 1);
        assert_eq!(face2.sub_entity(1, 1), 3);
        assert_eq!(face2.sub_entity(1, 2), 5);
        assert_eq!(face2.sub_entity(2, 0), 2);
        let face3 = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, 3);
        assert_eq!(face3.sub_entity(0, 0), 0);
        assert_eq!(face3.sub_entity(0, 1), 1);
        assert_eq!(face3.sub_entity(0, 2), 2);
        assert_eq!(face3.sub_entity(1, 0), 2);
        assert_eq!(face3.sub_entity(1, 1), 4);
        assert_eq!(face3.sub_entity(1, 2), 5);
        assert_eq!(face3.sub_entity(2, 0), 3);
        let face4 = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, 4);
        assert_eq!(face4.sub_entity(0, 0), 2);
        assert_eq!(face4.sub_entity(0, 1), 3);
        assert_eq!(face4.sub_entity(0, 2), 4);
        assert_eq!(face4.sub_entity(1, 0), 0);
        assert_eq!(face4.sub_entity(1, 1), 6);
        assert_eq!(face4.sub_entity(1, 2), 7);
        assert_eq!(face4.sub_entity(2, 0), 4);
        let face5 = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, 5);
        assert_eq!(face5.sub_entity(0, 0), 0);
        assert_eq!(face5.sub_entity(0, 1), 3);
        assert_eq!(face5.sub_entity(0, 2), 4);
        assert_eq!(face5.sub_entity(1, 0), 3);
        assert_eq!(face5.sub_entity(1, 1), 6);
        assert_eq!(face5.sub_entity(1, 2), 8);
        assert_eq!(face5.sub_entity(2, 0), 5);
        let face6 = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, 6);
        assert_eq!(face6.sub_entity(0, 0), 0);
        assert_eq!(face6.sub_entity(0, 1), 2);
        assert_eq!(face6.sub_entity(0, 2), 4);
        assert_eq!(face6.sub_entity(1, 0), 4);
        assert_eq!(face6.sub_entity(1, 1), 7);
        assert_eq!(face6.sub_entity(1, 2), 8);
        assert_eq!(face6.sub_entity(2, 0), 6);

        let edge0 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 0);
        assert_eq!(edge0.sub_entity(0, 0), 2);
        assert_eq!(edge0.sub_entity(0, 1), 3);
        assert_eq!(edge0.sub_entity(1, 0), 0);
        let edge1 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 1);
        assert_eq!(edge1.sub_entity(0, 0), 1);
        assert_eq!(edge1.sub_entity(0, 1), 3);
        assert_eq!(edge1.sub_entity(1, 0), 1);
        let edge2 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 2);
        assert_eq!(edge2.sub_entity(0, 0), 1);
        assert_eq!(edge2.sub_entity(0, 1), 2);
        assert_eq!(edge2.sub_entity(1, 0), 2);
        let edge3 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 3);
        assert_eq!(edge3.sub_entity(0, 0), 0);
        assert_eq!(edge3.sub_entity(0, 1), 3);
        assert_eq!(edge3.sub_entity(1, 0), 3);
        let edge4 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 4);
        assert_eq!(edge4.sub_entity(0, 0), 0);
        assert_eq!(edge4.sub_entity(0, 1), 2);
        assert_eq!(edge4.sub_entity(1, 0), 4);
        let edge5 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 5);
        assert_eq!(edge5.sub_entity(0, 0), 0);
        assert_eq!(edge5.sub_entity(0, 1), 1);
        assert_eq!(edge5.sub_entity(1, 0), 5);
        let edge6 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 6);
        assert_eq!(edge6.sub_entity(0, 0), 3);
        assert_eq!(edge6.sub_entity(0, 1), 4);
        assert_eq!(edge6.sub_entity(1, 0), 6);
        let edge7 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 7);
        assert_eq!(edge7.sub_entity(0, 0), 2);
        assert_eq!(edge7.sub_entity(0, 1), 4);
        assert_eq!(edge7.sub_entity(1, 0), 7);
        let edge8 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 8);
        assert_eq!(edge8.sub_entity(0, 0), 0);
        assert_eq!(edge8.sub_entity(0, 1), 4);
        assert_eq!(edge8.sub_entity(1, 0), 8);

        let vertex0 = SingleElementCellTopology::new(&t, ReferenceCellType::Point, 0);
        assert_eq!(vertex0.sub_entity(0, 0), 0);
        let vertex1 = SingleElementCellTopology::new(&t, ReferenceCellType::Point, 1);
        assert_eq!(vertex1.sub_entity(0, 0), 1);
        let vertex2 = SingleElementCellTopology::new(&t, ReferenceCellType::Point, 2);
        assert_eq!(vertex2.sub_entity(0, 0), 2);
        let vertex3 = SingleElementCellTopology::new(&t, ReferenceCellType::Point, 3);
        assert_eq!(vertex3.sub_entity(0, 0), 3);
        let vertex4 = SingleElementCellTopology::new(&t, ReferenceCellType::Point, 4);
        assert_eq!(vertex4.sub_entity(0, 0), 4);
    }

    #[test]
    fn test_sub_entity_iter_triangle() {
        //! Test sub-entities of a triangle
        let t = example_topology_triangle();
        for index in 0..2 {
            let cell = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, index);
            for dim in 0..3 {
                for (i, j) in cell.sub_entity_iter(dim).enumerate() {
                    assert_eq!(j, cell.sub_entity(dim, i));
                }
            }
        }
        for index in 0..5 {
            let edge = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, index);
            for dim in 0..2 {
                for (i, j) in edge.sub_entity_iter(dim).enumerate() {
                    assert_eq!(j, edge.sub_entity(dim, i));
                }
            }
        }
        for index in 0..4 {
            let vertex = SingleElementCellTopology::new(&t, ReferenceCellType::Point, index);
            for dim in 0..1 {
                for (i, j) in vertex.sub_entity_iter(dim).enumerate() {
                    assert_eq!(j, vertex.sub_entity(dim, i));
                }
            }
        }
    }

    #[test]
    fn test_sub_entity_iter_tetrahedron() {
        //! Test sub-entities of a tetrahedron
        let t = example_topology_tetrahedron();
        for index in 0..2 {
            let cell = SingleElementCellTopology::new(&t, ReferenceCellType::Tetrahedron, index);
            for dim in 0..4 {
                for (i, j) in cell.sub_entity_iter(dim).enumerate() {
                    assert_eq!(j, cell.sub_entity(dim, i));
                }
            }
        }
        for index in 0..7 {
            let triangle = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, index);
            for dim in 0..3 {
                for (i, j) in triangle.sub_entity_iter(dim).enumerate() {
                    assert_eq!(j, triangle.sub_entity(dim, i));
                }
            }
        }
        for index in 0..9 {
            let edge = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, index);
            for dim in 0..2 {
                for (i, j) in edge.sub_entity_iter(dim).enumerate() {
                    assert_eq!(j, edge.sub_entity(dim, i));
                }
            }
        }
        for index in 0..5 {
            let vertex = SingleElementCellTopology::new(&t, ReferenceCellType::Point, index);
            for dim in 0..1 {
                for (i, j) in vertex.sub_entity_iter(dim).enumerate() {
                    assert_eq!(j, vertex.sub_entity(dim, i));
                }
            }
        }
    }
}
