//! Grid builder

use super::SingleElementGrid;
use crate::{
    geometry::SingleElementGeometry,
    topology::single_type::SingleTypeTopology,
    traits::{Builder, GeometryBuilder, GridBuilder, TopologyBuilder},
    types::Scalar,
};
use ndelement::{
    ciarlet::{CiarletElement, LagrangeElementFamily, lagrange},
    map::IdentityMap,
    reference_cell,
    traits::FiniteElement,
    types::{Continuity, ReferenceCellType},
};
use rlst::rlst_dynamic_array;
use std::collections::{HashMap, HashSet};

/// Grid builder for a single element grid
///
/// The following gives an example of creating a new grid consisting
/// of a single triangle.
///
/// ```
/// use ndgrid::traits::Builder;
/// use ndgrid::SingleElementGridBuilder;
/// use ndelement::types::ReferenceCellType;
///
/// // The geometric dimension of our space is 3.
/// let gdim = 3;
///
/// // We are building a two dimensional surface triangle grid within a three dimensional space.
/// // Our grid will have three points and one `Triangle` cell of order 1.
/// let mut builder = SingleElementGridBuilder::new_with_capacity(gdim, 3, 1, (ReferenceCellType::Triangle, 1));
/// builder.add_point(0, &[0.0, 0.0, 0.0]);
/// builder.add_point(1, &[1.0, 0.0, 0.0]);
/// builder.add_point(2, &[0.0, 1.0, 0.0]);
/// builder.add_cell(0, &[0, 1, 2]);
///
/// let grid = builder.create_grid();
/// ```
#[derive(Debug)]
pub struct SingleElementGridBuilder<T: Scalar> {
    gdim: usize,
    element_data: (ReferenceCellType, usize),
    element: CiarletElement<T, IdentityMap, T>,
    points_per_cell: usize,
    pub(crate) points: Vec<T>,
    cells: Vec<usize>,
    point_indices_to_ids: Vec<usize>,
    point_ids_to_indices: HashMap<usize, usize>,
    cell_indices_to_ids: Vec<usize>,
    cell_ids_to_indices: HashMap<usize, usize>,
    point_indices: HashSet<usize>,
    cell_indices: HashSet<usize>,
}

impl<T: Scalar> SingleElementGridBuilder<T> {
    /// Create a new grid builder
    pub fn new(gdim: usize, data: (ReferenceCellType, usize)) -> Self {
        Self::new_with_capacity(gdim, 0, 0, data)
    }

    /// Create a new grid builder with capacaty for a given number of points and cells
    pub fn new_with_capacity(
        gdim: usize,
        npoints: usize,
        ncells: usize,
        data: (ReferenceCellType, usize),
    ) -> Self {
        let element = lagrange::create::<T, T>(data.0, data.1, Continuity::Standard);
        let points_per_cell = element.dim();
        Self {
            gdim,
            element_data: data,
            element,
            points_per_cell,
            points: Vec::with_capacity(npoints * gdim),
            cells: Vec::with_capacity(ncells * points_per_cell),
            point_indices_to_ids: Vec::with_capacity(npoints),
            point_ids_to_indices: HashMap::new(),
            cell_indices_to_ids: Vec::with_capacity(ncells),
            cell_ids_to_indices: HashMap::new(),
            point_indices: HashSet::new(),
            cell_indices: HashSet::new(),
        }
    }
}

impl<T: Scalar> Builder for SingleElementGridBuilder<T> {
    type Grid = SingleElementGrid<T, CiarletElement<T, IdentityMap, T>>;
    type T = T;
    type CellData<'a> = &'a [usize];
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

    fn add_cell(&mut self, id: usize, cell_data: &[usize]) {
        if self.cell_indices.contains(&id) {
            panic!("Cannot add cell with duplicate id.");
        }
        assert_eq!(cell_data.len(), self.points_per_cell);
        self.cell_ids_to_indices
            .insert(id, self.cell_indices_to_ids.len());
        self.cell_indices.insert(id);
        self.cell_indices_to_ids.push(id);
        for id in cell_data {
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
        if (cell_type, cell_degree) != self.element_data {
            panic!("Invalid cell type.");
        }
        self.add_cell(id, nodes);
    }

    fn create_grid(&self) -> SingleElementGrid<T, CiarletElement<T, IdentityMap, T>> {
        let cell_vertices =
            self.extract_vertices(&self.cells, &[self.element_data.0], &[self.element_data.1]);

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
            &[self.element_data.0],
            &[self.element_data.1],
        );

        let topology = self.create_topology(
            vertex_ids,
            (0..self.cell_count()).collect::<Vec<_>>(),
            &cell_vertices,
            &[self.element_data.0],
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
        &self.cell_indices_to_ids
    }
    fn cell_points(&self, index: usize) -> &[usize] {
        &self.cells[self.points_per_cell * index..self.points_per_cell * (index + 1)]
    }
    fn cell_vertices(&self, index: usize) -> &[usize] {
        &self.cells[self.points_per_cell * index
            ..self.points_per_cell * index + reference_cell::entity_counts(self.element_data.0)[0]]
    }
    fn point(&self, index: usize) -> &[T] {
        &self.points[self.gdim * index..self.gdim * (index + 1)]
    }
    fn points(&self) -> &[T] {
        &self.points
    }
    fn cell_type(&self, _index: usize) -> ReferenceCellType {
        self.element_data.0
    }
    fn cell_degree(&self, _index: usize) -> usize {
        self.element_data.1
    }
    fn gdim(&self) -> usize {
        self.gdim
    }
    fn tdim(&self) -> usize {
        reference_cell::dim(self.element_data.0)
    }

    fn npts(&self, _cell_type: Self::EntityDescriptor, _degree: usize) -> usize {
        self.points_per_cell
    }
}

impl<T: Scalar> GeometryBuilder for SingleElementGridBuilder<T> {
    type GridGeometry = SingleElementGeometry<T, CiarletElement<T, IdentityMap, T>>;
    fn create_geometry(
        &self,
        point_ids: &[usize],
        coordinates: &[Self::T],
        cell_points: &[usize],
        _cell_types: &[ReferenceCellType],
        _cell_degrees: &[usize],
    ) -> SingleElementGeometry<T, CiarletElement<T, IdentityMap, T>> {
        let npts = point_ids.len();
        let mut points = rlst_dynamic_array!(T, [self.gdim(), npts]);
        points.data_mut().unwrap().copy_from_slice(coordinates);

        let family = LagrangeElementFamily::<T, T>::new(self.element_data.1, Continuity::Standard);

        SingleElementGeometry::<T, CiarletElement<T, IdentityMap, T>>::new(
            self.element_data.0,
            points,
            cell_points,
            &family,
        )
    }
}

impl<T: Scalar> TopologyBuilder for SingleElementGridBuilder<T> {
    type GridTopology = SingleTypeTopology;
    fn create_topology(
        &self,
        vertex_ids: Vec<usize>,
        cell_ids: Vec<usize>,
        cells: &[usize],
        _cell_types: &[ReferenceCellType],
    ) -> SingleTypeTopology {
        // Create a map from point ids to the corresponding positions in the points array.
        let vertex_ids_to_pos = {
            let mut tmp = HashMap::<usize, usize>::new();
            for (i, id) in vertex_ids.iter().enumerate() {
                tmp.insert(*id, i);
            }
            tmp
        };

        // Remap the cell definitions to use internal indices instead of ids.
        let cells = {
            let mut new_cells = Vec::<usize>::with_capacity(cells.len());
            for id in cells {
                new_cells.push(vertex_ids_to_pos[id]);
            }
            new_cells
        };

        SingleTypeTopology::new(
            &cells,
            self.element_data.0,
            Some(vertex_ids),
            Some(cell_ids),
        )
    }

    fn extract_vertices(
        &self,
        cell_points: &[usize],
        _cell_types: &[Self::EntityDescriptor],
        _cell_degrees: &[usize],
    ) -> Vec<usize> {
        let mut vertices = vec![];
        for v in 0..reference_cell::entity_counts(self.element_data.0)[0] {
            for d in self.element.entity_dofs(0, v).unwrap() {
                vertices.push(*d);
            }
        }
        let mut cell_vertices = vec![];
        let mut start = 0;
        let npoints = self.element.dim();
        while start < cell_points.len() {
            for v in &vertices {
                cell_vertices.push(cell_points[start + v])
            }
            start += npoints;
        }
        cell_vertices
    }
}

impl<T: Scalar> GridBuilder for SingleElementGridBuilder<T> {
    fn create_grid_from_topology_geometry(
        &self,
        topology: <Self as TopologyBuilder>::GridTopology,
        geometry: <Self as GeometryBuilder>::GridGeometry,
    ) -> <Self as Builder>::Grid {
        SingleElementGrid::new(topology, geometry)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    #[should_panic]
    fn test_duplicate_point_id() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));

        b.add_point(2, &[0.0, 0.0, 0.0]);
        b.add_point(0, &[1.0, 0.0, 0.0]);
        b.add_point(1, &[0.0, 1.0, 0.0]);
        b.add_point(2, &[1.0, 1.0, 0.0]);
    }

    #[test]
    #[should_panic]
    fn test_duplicate_cell_id() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));

        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0, 0.0]);
        b.add_point(2, &[0.0, 1.0, 0.0]);
        b.add_point(3, &[1.0, 1.0, 0.0]);

        b.add_cell(0, &[0, 1, 2]);
        b.add_cell(0, &[1, 2, 3]);
    }

    #[test]
    fn test_non_contiguous_ids() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));

        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0, 0.0]);
        b.add_point(2, &[0.0, 1.0, 0.0]);
        b.add_point(4, &[1.0, 1.0, 0.0]);

        b.add_cell(0, &[0, 1, 2]);
        b.add_cell(2, &[1, 2, 4]);

        b.create_grid();
    }
}
