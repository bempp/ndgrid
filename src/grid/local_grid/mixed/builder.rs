//! Grid builder

use super::MixedGrid;
use crate::{
    geometry::MixedGeometry,
    topology::mixed::MixedTopology,
    traits::{Builder, GeometryBuilder, GridBuilder, TopologyBuilder},
    types::Scalar,
};
use itertools::izip;
use ndelement::{
    ciarlet::{CiarletElement, LagrangeElementFamily, lagrange},
    map::IdentityMap,
    reference_cell,
    traits::FiniteElement,
    types::{Continuity, ReferenceCellType},
};
use rlst::rlst_dynamic_array;
use std::collections::{HashMap, HashSet};

/// Grid builder for a grid with a mixture of element types
#[derive(Debug)]
pub struct MixedGridBuilder<T: Scalar> {
    gdim: usize,
    element_indices: HashMap<(ReferenceCellType, usize), usize>,
    elements: Vec<CiarletElement<T, IdentityMap, T>>,
    points_per_cell: Vec<usize>,
    pub(crate) points: Vec<T>,
    cells: Vec<usize>,
    cell_starts: Vec<usize>,
    cell_families: Vec<usize>,
    cell_degrees: Vec<usize>,
    cell_types: Vec<ReferenceCellType>,
    point_indices_to_ids: Vec<usize>,
    point_ids_to_indices: HashMap<usize, usize>,
    cell_indices_to_ids: Vec<usize>,
    cell_ids_to_indices: HashMap<usize, usize>,
    point_indices: HashSet<usize>,
    cell_indices: HashSet<usize>,
}

impl<T: Scalar> MixedGridBuilder<T> {
    /// Create a new grid builder
    pub fn new(gdim: usize) -> Self {
        Self {
            gdim,
            element_indices: HashMap::new(),
            elements: vec![],
            points_per_cell: vec![],
            points: vec![],
            cells: vec![],
            cell_starts: vec![],
            cell_families: vec![],
            cell_degrees: vec![],
            cell_types: vec![],
            point_indices_to_ids: vec![],
            point_ids_to_indices: HashMap::new(),
            cell_indices_to_ids: vec![],
            cell_ids_to_indices: HashMap::new(),
            point_indices: HashSet::new(),
            cell_indices: HashSet::new(),
        }
    }
}

impl<T: Scalar> Builder for MixedGridBuilder<T> {
    type Grid = MixedGrid<T, CiarletElement<T, IdentityMap, T>>;
    type T = T;
    type CellData<'a> = (ReferenceCellType, usize, &'a [usize]);
    type EntityDescriptor = ReferenceCellType;

    fn add_point(&mut self, id: usize, data: &[T]) {
        if data.len() != self.gdim {
            panic!("Point has wrong number of coordinates");
        }
        if self.point_indices.contains(&id) {
            panic!("Cannot add point with duplicate id.");
        }
        self.point_ids_to_indices
            .insert(id, self.point_indices_to_ids.len());
        self.point_indices.insert(id);
        self.point_indices_to_ids.push(id);
        self.points.extend_from_slice(data);
    }

    fn add_cell(&mut self, id: usize, cell_data: (ReferenceCellType, usize, &[usize])) {
        let element_index = *self
            .element_indices
            .entry((cell_data.0, cell_data.1))
            .or_insert_with(|| {
                let n = self.cell_indices_to_ids.len();
                self.elements.push(lagrange::create::<T, T>(
                    cell_data.0,
                    cell_data.1,
                    Continuity::Standard,
                ));
                self.points_per_cell.push(self.elements[n].dim());
                n
            });
        self.cell_families.push(element_index);
        self.cell_starts.push(self.cells.len());
        self.cell_types.push(cell_data.0);
        if reference_cell::dim(cell_data.0) != reference_cell::dim(self.cell_types[0]) {
            panic!("Cannot create grid with cells of different topological dimensions");
        }
        self.cell_degrees.push(cell_data.1);
        if self.cell_indices.contains(&id) {
            panic!("Cannot add cell with duplicate id.");
        }
        assert_eq!(cell_data.2.len(), self.points_per_cell[element_index]);
        self.cell_ids_to_indices
            .insert(id, self.cell_indices_to_ids.len());
        self.cell_indices.insert(id);
        self.cell_indices_to_ids.push(id);
        for id in cell_data.2 {
            self.cells.push(self.point_ids_to_indices[id]);
        }
    }

    fn add_cell_from_nodes_and_type(
        &mut self,
        id: usize,
        nodes: &[usize],
        cell_type: ReferenceCellType,
        cell_degree: usize,
    ) {
        self.add_cell(id, (cell_type, cell_degree, nodes));
    }

    fn create_grid(&self) -> MixedGrid<T, CiarletElement<T, IdentityMap, T>> {
        let cell_vertices =
            self.extract_vertices(&self.cells, &self.cell_types, &self.cell_degrees);

        // Add the vertex ids. But need to make sure that we don't have duplicates.
        // So first generate a set of vertex ids, then convert it to a vector.
        // We finally sort the vector to make sure that grid generation is reproducible across runs
        // as sets do not have a stable ordering.
        let mut vertex_ids = HashSet::<usize>::new();
        for v in &cell_vertices {
            vertex_ids.insert(*v);
        }

        let vertex_ids = {
            let mut tmp = Vec::<usize>::with_capacity(vertex_ids.len());
            tmp.extend(vertex_ids.iter());
            tmp.sort();
            tmp
        };

        let geometry = self.create_geometry(
            &self.point_indices_to_ids,
            &self.points,
            &self.cells,
            &self.cell_types,
            &self.cell_degrees,
        );

        let topology = self.create_topology(
            vertex_ids,
            (0..self.cell_count()).collect::<Vec<_>>(),
            &cell_vertices,
            &self.cell_types,
        );

        self.create_grid_from_topology_geometry(topology, geometry)
    }

    fn point_count(&self) -> usize {
        self.point_indices_to_ids.len()
    }
    fn cell_count(&self) -> usize {
        self.cell_indices_to_ids.len()
    }
    fn point_indices_to_ids(&self) -> &[usize] {
        &self.point_indices_to_ids
    }
    fn cell_indices_to_ids(&self) -> &[usize] {
        unimplemented!();
        // &self.cell_indices_to_ids
    }
    fn cell_points(&self, index: usize) -> &[usize] {
        &self.cells[self.cell_starts[index]
            ..self.cell_starts[index] + self.points_per_cell[self.cell_families[index]]]
    }
    fn cell_vertices(&self, index: usize) -> &[usize] {
        &self.cells[self.cell_starts[index]
            ..self.cell_starts[index] + reference_cell::entity_counts(self.cell_types[index])[0]]
    }
    fn point(&self, index: usize) -> &[T] {
        &self.points[self.gdim * index..self.gdim * (index + 1)]
    }
    fn points(&self) -> &[T] {
        &self.points
    }
    fn cell_type(&self, index: usize) -> ReferenceCellType {
        self.cell_types[index]
    }
    fn cell_degree(&self, index: usize) -> usize {
        self.cell_degrees[index]
    }
    fn gdim(&self) -> usize {
        self.gdim
    }
    fn tdim(&self) -> usize {
        reference_cell::dim(self.cell_types[0])
    }
    fn npts(&self, cell_type: Self::EntityDescriptor, degree: usize) -> usize {
        self.points_per_cell[self.element_indices[&(cell_type, degree)]]
    }
}

impl<T: Scalar> GeometryBuilder for MixedGridBuilder<T> {
    type GridGeometry = MixedGeometry<T, CiarletElement<T, IdentityMap, T>>;
    fn create_geometry(
        &self,
        point_ids: &[usize],
        coordinates: &[Self::T],
        cell_points: &[usize],
        cell_types: &[ReferenceCellType],
        cell_degrees: &[usize],
    ) -> MixedGeometry<T, CiarletElement<T, IdentityMap, T>> {
        let npts = point_ids.len();
        let mut points = rlst_dynamic_array!(T, [self.gdim(), npts]);
        points.data_mut().unwrap().copy_from_slice(coordinates);

        let mut element_families = vec![];
        let mut ef_indices = HashMap::new();
        let mut cell_families = vec![];
        for degree in cell_degrees {
            cell_families.push(*ef_indices.entry(*degree).or_insert_with(|| {
                let n = element_families.len();
                element_families.push(LagrangeElementFamily::<T, T>::new(
                    *degree,
                    Continuity::Standard,
                ));
                n
            }))
        }

        MixedGeometry::<T, CiarletElement<T, IdentityMap, T>>::new(
            cell_types,
            points,
            cell_points,
            &element_families,
            &cell_families,
        )
    }
}

impl<T: Scalar> TopologyBuilder for MixedGridBuilder<T> {
    type GridTopology = MixedTopology;
    fn create_topology(
        &self,
        vertex_ids: Vec<usize>,
        cell_ids: Vec<usize>,
        cells: &[usize],
        _cell_types: &[ReferenceCellType],
    ) -> MixedTopology {
        // Create a map from point ids to the corresponding positions in the points array.
        let vertex_ids_to_pos = {
            let mut tmp = HashMap::<usize, usize>::new();
            for (i, id) in vertex_ids.iter().enumerate() {
                tmp.insert(*id, i);
            }
            tmp
        };

        let cells = {
            let mut new_cells = Vec::<usize>::with_capacity(cells.len());
            for id in cells {
                new_cells.push(vertex_ids_to_pos[id]);
            }
            new_cells
        };

        MixedTopology::new(&cells, &self.cell_types, Some(vertex_ids), Some(cell_ids))
    }

    fn extract_vertices(
        &self,
        cell_points: &[usize],
        cell_types: &[Self::EntityDescriptor],
        cell_degrees: &[usize],
    ) -> Vec<usize> {
        let mut cell_vertices = vec![];
        let mut start = 0;
        for (t, d) in izip!(cell_types, cell_degrees) {
            let e = &self.elements[self.element_indices[&(*t, *d)]];
            for v in 0..reference_cell::entity_counts(*t)[0] {
                for d in e.entity_dofs(0, v).unwrap() {
                    cell_vertices.push(cell_points[start + *d])
                }
            }
            start += e.dim();
        }
        cell_vertices
    }
}

impl<T: Scalar> GridBuilder for MixedGridBuilder<T> {
    fn create_grid_from_topology_geometry(
        &self,
        topology: <Self as TopologyBuilder>::GridTopology,
        geometry: <Self as GeometryBuilder>::GridGeometry,
    ) -> <Self as Builder>::Grid {
        MixedGrid::new(topology, geometry)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    #[should_panic]
    fn test_duplicate_point_id() {
        let mut b = MixedGridBuilder::<f64>::new(3);

        b.add_point(2, &[0.0, 0.0, 0.0]);
        b.add_point(0, &[1.0, 0.0, 0.0]);
        b.add_point(1, &[0.0, 1.0, 0.0]);
        b.add_point(2, &[1.0, 1.0, 0.0]);
    }

    #[test]
    #[should_panic]
    fn test_duplicate_cell_id() {
        let mut b = MixedGridBuilder::<f64>::new(3);

        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0, 0.0]);
        b.add_point(2, &[0.0, 1.0, 0.0]);
        b.add_point(3, &[1.0, 1.0, 0.0]);

        b.add_cell(0, (ReferenceCellType::Triangle, 1, &[0, 1, 2]));
        b.add_cell(0, (ReferenceCellType::Triangle, 1, &[1, 2, 3]));
    }

    #[test]
    fn test_non_contiguous_ids() {
        let mut b = MixedGridBuilder::<f64>::new(3);

        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0, 0.0]);
        b.add_point(2, &[0.0, 1.0, 0.0]);
        b.add_point(4, &[1.0, 1.0, 0.0]);

        b.add_cell(0, (ReferenceCellType::Triangle, 1, &[0, 1, 2]));
        b.add_cell(2, (ReferenceCellType::Triangle, 1, &[1, 2, 4]));

        b.create_grid();
    }
}
