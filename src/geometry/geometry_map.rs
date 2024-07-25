//! Geometry map
use crate::{traits::GeometryMap as GeometryMapTrait, types::{Array2D, ArrayND, RealScalar}};
use ndelement::{
    ciarlet::CiarletElement,
    reference_cell,
    traits::{ElementFamily, FiniteElement},
    types::ReferenceCellType,
};
use rlst::{rlst_dynamic_array2, RawAccessMut, Shape, rlst_dynamic_array4, RandomAccessByRef};

/// Single element geometry
pub struct GeometryMap<'a, T: RealScalar> {
    geometry_points: &'a Array2D<T>,
    entities: &'a Array2D<usize>,
    tdim: usize,
    gdim: usize,
    table: ArrayND<4, T>,
}

impl<'a, T: RealScalar> GeometryMap<'a, T> {
    /// Create new
    pub fn new<A2D: RandomAccessByRef<2, Item=T> + Shape<2>>(
        element: &impl FiniteElement<CellType=ReferenceCellType, T=T>,
        points: &A2D,
        geometry_points: &'a Array2D<T>,
        entities: &'a Array2D<usize>,
    ) -> Self {
        let tdim = reference_cell::dim(element.cell_type());
        let gdim = points.shape()[0];
        let npoints = points.shape()[1];

        let mut table = rlst_dynamic_array4!(T, element.tabulate_array_shape(1, npoints));
        element.tabulate(points, 1, &mut table);

        Self {
            geometry_points,
            entities,
            tdim,
            gdim,
            table,
        }
    }
}

impl<'a, T: RealScalar> GeometryMapTrait for GeometryMap<'a, T> {
    type T = T;

    fn entity_topology_dimension(&self) -> usize{
        self.tdim
    }
    fn geometry_dimension(&self) -> usize{
        self.gdim
    }
    fn point_count(&self) -> usize{
        self.table.shape()[1]
    }
    fn points(&self, cell_index: usize, value: &mut [T]){
        todo!();/*
        let gdim = geometry.dim();
        let npts = table.shape()[1];
        assert_eq!(points.len(), geometry.dim() * npts);

        let cell = geometry.index_map()[cell_index];

        for component in points.iter_mut() {
            *component = T::from(0.0).unwrap();
        }
        for (i, v) in geometry.cell_points(cell).unwrap().iter().enumerate() {
            for point_index in 0..npts {
                let t = unsafe { *table.get_unchecked([0, point_index, i, 0]) };
                for j in 0..gdim {
                    points[j * npts + point_index] += *geometry.coordinate(*v, j).unwrap() * t;
                }
            }
        }
        */
    }
    fn jacobian(&self, cell_index: usize, value: &mut [T]){
        unimplemented!();
    }
    fn normal(&self, cell_index: usize, value: &mut [T]){
        if self.tdim + 1 != self.gdim {
            panic!("Can only compute normal for entities where tdim + 1 == gdim");
        }
        println!("{} {}", self.tdim, self.gdim);
        unimplemented!();
    }
}
