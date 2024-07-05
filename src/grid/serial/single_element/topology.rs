    //! Implementation of grid topology

use crate::traits::Topology;
use crate::types::{CellLocalIndexPair, Ownership, Array2D};
use ndelement::reference_cell;
use ndelement::types::ReferenceCellType;
use std::collections::HashMap;
use rlst::{rlst_dynamic_array2, RawAccess, Shape};

fn all_equal<T: Eq>(a: &[T], b: &[T]) -> bool {
    if a.len() != b.len() {
        false
    } else {
        all_in(a, b)
    }
}

fn all_in<T: Eq>(a: &[T], b: &[T]) -> bool {
    for i in a {
        if !b.contains(i) {
            return false;
        }
    }
    true
}

/// Topology of a single element grid
pub struct SingleElementTopology {
    dim: usize,
//    index_map: Vec<usize>,
    //entity_types: Vec<ReferenceCellType>,
    pub(crate) downward_connectivity: Vec<Vec<Array2D<usize>>>,
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
        ReferenceCellType::Point => {},
        ReferenceCellType::Interval => {
            if vertices[0] > vertices[1] {
                vertices.swap(0, 1);
            }
        },
        _ => { unimplemented!(); }
    }
}

impl SingleElementTopology {
    /// Create a topology
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        cells: &[usize],
        cell_type: ReferenceCellType,
        point_ids: Option<&[usize]>,
        cell_ids: Option<&[usize]>,
    ) -> Self {
        // Cells where faces are mixture of triangles and quads not supported
        if cell_type != ReferenceCellType::Point && cell_type != ReferenceCellType::Triangle && cell_type != ReferenceCellType::Quadrilateral && cell_type != ReferenceCellType::Tetrahedron && cell_type != ReferenceCellType::Hexahedron {
            panic!("Unsupported cell type for SingleElementTopology: {cell_type:?}");
        }
        let size = reference_cell::entity_counts(cell_type)[0];
        let ncells = cells.len() / size;
        let dim = reference_cell::dim(cell_type);
        let ref_conn = reference_cell::connectivity(cell_type);
        let etypes = reference_cell::entity_types(cell_type);

        let mut nentities = vec![0; dim + 1];
        nentities[0] = cells.iter().max().unwrap() + 1;
        nentities[dim] = ncells;
        let mut entities = vec![vec![]; dim - 1];
        for cell_i in 0..ncells {
            let cell = cells[cell_i * size..(cell_i+1)*size].to_vec();
            for d in 1..dim {
                for (c, et) in ref_conn[d].iter().zip(&etypes[d]) {
                    let mut entity = c[0].iter().map(|i| cell[*i]).collect::<Vec<_>>();
                    orient_entity(*et, &mut entity);
                    if !entities[d - 1].contains(&entity) {
                        entities[d - 1].push(entity);
                    }
                }
            }
        }
        for d in 1..dim {
            nentities[d] = entities[d - 1].len();
        }

        let mut downward_connectivity = nentities.iter().enumerate().map(
            |(i, j)| ref_conn[i][0].iter().take(i + 1).map(|r| rlst_dynamic_array2!(usize, [r.len(), *j])).collect::<Vec<_>>()
        ).collect::<Vec<_>>();

        for d in 0..dim + 1 {
            for (i, mut j) in downward_connectivity[d][d].col_iter_mut().enumerate() {
                j[[0]] = i;
            }
        }

        println!("{:?}", entities);
        println!("<<");
        for d in 0..dim + 1 {
            for i in &downward_connectivity[d] {
                println!("{:?}", i.data());
            }
            println!();
        }
        println!(">>");

        Self {
            dim,
            downward_connectivity,
        }
//        panic!();
/*
        

        let mut vertex_indices_to_ids = vec![];
        let mut vertex_ids_to_indices = HashMap::new();
        let mut cell_indices_to_ids = vec![];
        let mut cell_ids_to_indices = HashMap::new();

        let mut index_map = vec![0; ncells];
        let mut vertices = vec![];

        let entity_types = reference_cell::entity_types(cell_type)
            .iter()
            .filter(|t| !t.is_empty())
            .map(|t| t[0])
            .collect::<Vec<_>>();

        let mut cells_to_entities = vec![vec![vec![]; ncells]; dim + 1];
        let mut entities_to_cells = vec![vec![]; dim + 1];
        let mut entities_to_vertices = vec![vec![]; dim];

        entities_to_cells[dim] = vec![vec![]; ncells];

        let mut start = 0;
        for (cell_i, i) in index_map.iter_mut().enumerate() {
            let cell = &cells_input[start..start + size];
            *i = cell_i;
            cell_indices_to_ids.push(grid_cell_indices_to_ids[cell_i]);
            cell_ids_to_indices.insert(grid_cell_indices_to_ids[cell_i], cell_i);
            let mut row = vec![];
            for v in cell {
                if !vertices.contains(v) {
                    entities_to_cells[0].push(vec![]);
                    vertex_indices_to_ids.push(point_indices_to_ids[*v]);
                    vertex_ids_to_indices.insert(point_indices_to_ids[*v], vertices.len());
                    vertices.push(*v);
                }
                row.push(vertices.iter().position(|&r| r == *v).unwrap());
            }

            for (local_index, v) in row.iter().enumerate() {
                entities_to_cells[0][*v].push(CellLocalIndexPair::new(cell_i, local_index));
            }
            entities_to_cells[dim][cell_i] = vec![CellLocalIndexPair::new(cell_i, 0)];

            cells_to_entities[0][cell_i] = row;
            cells_to_entities[dim][cell_i] = vec![cell_i];

            start += size;
        }

        entities_to_vertices[0] = (0..vertices.len()).map(|i| vec![i]).collect::<Vec<_>>();

        let mut edge_indices = HashMap::new();
        let mut edge_indices_to_ids = vec![];
        let mut edge_ids_to_indices = HashMap::new();
        if let Some(e) = &edge_ids {
            for (edge_i, (i, j)) in e.iter().enumerate() {
                let mut v0 = vertex_ids_to_indices[&i[0]];
                let mut v1 = vertex_ids_to_indices[&i[1]];
                if v0 > v1 {
                    std::mem::swap(&mut v0, &mut v1);
                }
                edge_indices.insert((v0, v1), edge_i);
                edge_indices_to_ids.push(*j);
                edge_ids_to_indices.insert(*j, edge_i);
                entities_to_vertices[1].push(vec![v0, v1]);
                entities_to_cells[1].push(vec![]);
            }
        }
        let ref_conn = &reference_cell::connectivity(cell_type)[1];
        for cell_i in 0..ncells {
            for (local_index, rc) in ref_conn.iter().enumerate() {
                let cell = &cells_to_entities[0][cell_i];
                let mut first = cell[rc[0][0]];
                let mut second = cell[rc[0][1]];
                if first > second {
                    std::mem::swap(&mut first, &mut second);
                }
                if let Some(edge_index) = edge_indices.get(&(first, second)) {
                    cells_to_entities[1][cell_i].push(*edge_index);
                    entities_to_cells[1][*edge_index]
                        .push(CellLocalIndexPair::new(cell_i, local_index));
                } else {
                    if edge_ids.is_some() {
                        panic!("Missing id for edge");
                    }
                    let id = entities_to_vertices[1].len();
                    edge_indices.insert((first, second), id);
                    edge_indices_to_ids.push(id);
                    edge_ids_to_indices.insert(id, id);
                    cells_to_entities[1][cell_i].push(entities_to_vertices[1].len());
                    entities_to_cells[1].push(vec![CellLocalIndexPair::new(cell_i, local_index)]);
                    entities_to_vertices[1].push(vec![first, second]);
                }
            }
        }

        for d in 2..dim {
            let mut c_to_e = vec![];
            let ref_conn = &reference_cell::connectivity(cell_type)[d];
            for (cell_i, cell) in cells_to_entities[0].iter().enumerate() {
                let mut entity_ids = vec![];
                for (local_index, rc) in ref_conn.iter().enumerate() {
                    let vertices = rc[0].iter().map(|x| cell[*x]).collect::<Vec<_>>();
                    let mut found = false;
                    for (entity_index, entity) in entities_to_vertices[d].iter().enumerate() {
                        if all_equal(entity, &vertices) {
                            entity_ids.push(entity_index);
                            entities_to_cells[d][entity_index]
                                .push(CellLocalIndexPair::new(cell_i, local_index));
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        entity_ids.push(entities_to_vertices[d].len());
                        entities_to_cells[d]
                            .push(vec![CellLocalIndexPair::new(cell_i, local_index)]);
                        entities_to_vertices[d].push(vertices);
                    }
                }
                c_to_e.push(entity_ids);
            }
            cells_to_entities[d] = c_to_e;
        }

        Self {
            dim,
            index_map,
            entities_to_vertices,
            cells_to_entities,
            entities_to_cells,
            entity_types,
            vertex_indices_to_ids,
            edge_ids_to_indices,
            edge_indices_to_ids,
            vertex_ids_to_indices,
            cell_indices_to_ids,
            cell_ids_to_indices,
            cell_types: [cell_type],
        }
        */
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

struct SingleElementCellTopology<'a> {
    topology: &'a SingleElementTopology,
    entity_type: ReferenceCellType,
    entity_index: usize,
    dim: usize,
}

impl<'t> SingleElementCellTopology<'t> {
    pub fn new(topology: &'t SingleElementTopology,
    entity_type: ReferenceCellType,
    entity_index: usize) -> Self {
        Self {
            topology, entity_type, entity_index, dim: reference_cell::dim(entity_type)
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

    fn connected_entity_iter(&self, dim: usize) -> IndexIter {
        unimplemented!();
    }

    fn sub_entity_iter(&self, dim: usize) -> Copied<std::slice::Iter<'_, usize>> {
        let rows = self.topology.downward_connectivity[self.dim][dim].shape()[0];
        self.topology.downward_connectivity[self.dim][dim].data()[rows * self.entity_index..rows * (self.entity_index + 1)].iter().copied()
    }

    fn sub_entity(&self, dim: usize, index: usize) -> usize {
        self.topology.downward_connectivity[self.dim][dim][[index, self.entity_index]]
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn example_topology() -> SingleElementTopology {
        //! An example topology
        SingleElementTopology::new(
            &[0, 1, 2, 2, 1, 3],
            ReferenceCellType::Triangle,
            None,
            None,
        )
    }

    #[test]
    fn test_sub_entities() {
        //! Test entity counts
        let t = example_topology();
        let cell0 = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, 0);
        //assert_eq!(cell0.sub_entity(0, 0), 0);
        //assert_eq!(cell0.sub_entity(0, 1), 1);
        //assert_eq!(cell0.sub_entity(0, 2), 2);
        //assert_eq!(cell0.sub_entity(1, 0), 0);
        //assert_eq!(cell0.sub_entity(1, 1), 1);
        //assert_eq!(cell0.sub_entity(1, 2), 2);
        assert_eq!(cell0.sub_entity(2, 0), 0);
        let cell1 = SingleElementCellTopology::new(&t, ReferenceCellType::Triangle, 1);
        //assert_eq!(cell1.sub_entity(0, 0), 2);
        //assert_eq!(cell1.sub_entity(0, 1), 1);
        //assert_eq!(cell1.sub_entity(0, 2), 3);
        //assert_eq!(cell1.sub_entity(1, 0), 3);
        //assert_eq!(cell1.sub_entity(1, 1), 4);
        //assert_eq!(cell1.sub_entity(1, 2), 0);
        assert_eq!(cell1.sub_entity(2, 0), 1);

        let edge0 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 0);
        //assert_eq!(edge0.sub_entity(0, 0), 1);
        //assert_eq!(edge0.sub_entity(0, 1), 2);
        assert_eq!(edge0.sub_entity(1, 0), 0);
        let edge1 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 1);
        //assert_eq!(edge1.sub_entity(0, 0), 0);
        //assert_eq!(edge1.sub_entity(0, 1), 2);
        assert_eq!(edge1.sub_entity(1, 0), 1);
        let edge2 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 2);
        //assert_eq!(edge2.sub_entity(0, 0), 0);
        //assert_eq!(edge2.sub_entity(0, 1), 1);
        assert_eq!(edge2.sub_entity(1, 0), 2);
        let edge3 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 3);
        //assert_eq!(edge3.sub_entity(0, 0), 1);
        //assert_eq!(edge3.sub_entity(0, 1), 3);
        assert_eq!(edge3.sub_entity(1, 0), 3);
        let edge4 = SingleElementCellTopology::new(&t, ReferenceCellType::Interval, 4);
        //assert_eq!(edge4.sub_entity(0, 0), 2);
        //assert_eq!(edge4.sub_entity(0, 1), 3);
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
    fn test_sub_entity_iter() {
        //! Test entity counts
        let t = example_topology();
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


/*

    #[test]
    fn test_cell_entities_vertices() {
        //! Test cell vertices
        let t = example_topology();
        for (i, vertices) in [[0, 1, 2], [2, 1, 3]].iter().enumerate() {
            let c = t.cell_to_entities(i, 0).unwrap();
            assert_eq!(c.len(), 3);
            assert_eq!(c, vertices);
        }
    }

    #[test]
    fn test_cell_entities_edges() {
        //! Test cell edges
        let t = example_topology();
        for (i, edges) in [[0, 1, 2], [3, 4, 0]].iter().enumerate() {
            let c = t.cell_to_entities(i, 1).unwrap();
            assert_eq!(c.len(), 3);
            assert_eq!(c, edges);
        }
    }

    #[test]
    fn test_cell_entities_cells() {
        //! Test cells
        let t = example_topology();
        for i in 0..2 {
            let c = t.cell_to_entities(i, 2).unwrap();
            assert_eq!(c.len(), 1);
            assert_eq!(c[0], i);
        }
    }

    #[test]
    fn test_entities_to_cells_vertices() {
        //! Test vertex-to-cell connectivity
        let t = example_topology();
        let c_to_e = (0..t.entity_count(ReferenceCellType::Triangle))
            .map(|i| t.cell_to_entities(i, 0).unwrap())
            .collect::<Vec<_>>();
        let e_to_c = (0..t.entity_count(ReferenceCellType::Point))
            .map(|i| {
                t.entity_to_cells(0, i)
                    .unwrap()
                    .iter()
                    .map(|x| x.cell)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        for (i, cell) in c_to_e.iter().enumerate() {
            for v in *cell {
                assert!(e_to_c[*v].contains(&i));
            }
        }
        for (i, cells) in e_to_c.iter().enumerate() {
            for c in cells {
                assert!(c_to_e[*c].contains(&i));
            }
        }
    }

    #[test]
    fn test_entities_to_cells_edges() {
        //! Test edge-to-cell connectivity
        let t = example_topology();
        let c_to_e = (0..t.entity_count(ReferenceCellType::Triangle))
            .map(|i| t.cell_to_entities(i, 1).unwrap())
            .collect::<Vec<_>>();
        let e_to_c = (0..t.entity_count(ReferenceCellType::Interval))
            .map(|i| {
                t.entity_to_cells(1, i)
                    .unwrap()
                    .iter()
                    .map(|x| x.cell)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        for (i, cell) in c_to_e.iter().enumerate() {
            for v in *cell {
                assert!(e_to_c[*v].contains(&i));
            }
        }
        for (i, cells) in e_to_c.iter().enumerate() {
            for c in cells {
                assert!(c_to_e[*c].contains(&i));
            }
        }
    }
    */
}
