//! Binding for C
#![allow(missing_docs)]
#![allow(clippy::missing_safety_doc)]

#[derive(Debug, PartialEq, Clone, Copy)]
#[repr(u8)]
pub enum DType {
    F32 = 0,
    F64 = 1,
}

mod grid {
    use super::DType;
    use crate::grid::SingleElementGrid;
    use crate::traits::Grid;
    use crate::types::RealScalar;
    use ndelement::{
        ciarlet::CiarletElement,
        ciarlet::LagrangeElementFamily,
        traits::{ElementFamily, FiniteElement},
        types::{Continuity, ReferenceCellType},
    };
    use std::ffi::c_void;
    use std::slice::from_raw_parts;

    #[derive(Debug, PartialEq, Clone, Copy)]
    #[repr(u8)]
    pub enum GridType {
        SerialSingleElementGrid = 0,
    }

    #[repr(C)]
    pub struct GridWrapper {
        pub grid: *const c_void,
        pub dtype: DType,
        pub gtype: GridType,
    }

    impl Drop for GridWrapper {
        fn drop(&mut self) {
            let Self { grid, dtype, gtype } = self;
            match gtype {
                GridType::SerialSingleElementGrid => match dtype {
                    DType::F32 => drop(unsafe {
                        Box::from_raw(*grid as *mut SingleElementGrid<f32, CiarletElement<f32>>)
                    }),
                    DType::F64 => drop(unsafe {
                        Box::from_raw(*grid as *mut SingleElementGrid<f64, CiarletElement<f64>>)
                    }),
                },
            }
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_free_grid(g: *mut GridWrapper) {
        assert!(!g.is_null());
        unsafe { drop(Box::from_raw(g)) }
    }

    #[allow(clippy::too_many_arguments)]
    unsafe fn single_element_grid_new_from_raw_data<T: RealScalar>(
        coordinates: *const T,
        npts: usize,
        gdim: usize,
        cells: *const usize,
        ncells: usize,
        cell_type: u8,
        geometry_degree: usize,
        dtype: DType,
    ) -> *const GridWrapper {
        let coordinates = from_raw_parts(coordinates, npts * gdim);
        let points_per_cell =
            LagrangeElementFamily::<T>::new(geometry_degree, Continuity::Standard)
                .element(ReferenceCellType::from(cell_type).unwrap())
                .dim();
        let cells = from_raw_parts(cells, ncells * points_per_cell);

        let grid = Box::into_raw(Box::new(
            SingleElementGrid::<T, CiarletElement<T>>::new_from_raw_data(
                coordinates,
                gdim,
                cells,
                ReferenceCellType::from(cell_type).unwrap(),
                geometry_degree,
            ),
        )) as *const c_void;
        Box::into_raw(Box::new(GridWrapper {
            grid,
            dtype,
            gtype: GridType::SerialSingleElementGrid,
        }))
    }

    #[no_mangle]
    pub unsafe extern "C" fn single_element_grid_new_from_raw_data_f32(
        coordinates: *const f32,
        npts: usize,
        gdim: usize,
        cells: *const usize,
        ncells: usize,
        cell_type: u8,
        geometry_degree: usize,
    ) -> *const GridWrapper {
        single_element_grid_new_from_raw_data::<f32>(
            coordinates,
            npts,
            gdim,
            cells,
            ncells,
            cell_type,
            geometry_degree,
            DType::F32,
        )
    }

    #[no_mangle]
    pub unsafe extern "C" fn single_element_grid_new_from_raw_data_f64(
        coordinates: *const f64,
        npts: usize,
        gdim: usize,
        cells: *const usize,
        ncells: usize,
        cell_type: u8,
        geometry_degree: usize,
    ) -> *const GridWrapper {
        single_element_grid_new_from_raw_data::<f64>(
            coordinates,
            npts,
            gdim,
            cells,
            ncells,
            cell_type,
            geometry_degree,
            DType::F64,
        )
    }

    unsafe fn extract_grid<G: Grid>(grid: *const GridWrapper) -> *const G {
        (*grid).grid as *const G
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_tdim(grid: *const GridWrapper) -> usize {
        match (*grid).gtype {
            GridType::SerialSingleElementGrid => match (*grid).dtype {
                DType::F32 => (*extract_grid::<SingleElementGrid<f32, CiarletElement<f32>>>(grid))
                    .topology_dim(),
                DType::F64 => (*extract_grid::<SingleElementGrid<f64, CiarletElement<f64>>>(grid))
                    .topology_dim(),
            },
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_gdim(grid: *const GridWrapper) -> usize {
        match (*grid).gtype {
            GridType::SerialSingleElementGrid => match (*grid).dtype {
                DType::F32 => (*extract_grid::<SingleElementGrid<f32, CiarletElement<f32>>>(grid))
                    .geometry_dim(),
                DType::F64 => (*extract_grid::<SingleElementGrid<f64, CiarletElement<f64>>>(grid))
                    .geometry_dim(),
            },
        }
    }

    #[no_mangle]
    pub unsafe extern "C" fn grid_dtype(grid: *const GridWrapper) -> u8 {
        (*grid).dtype as u8
    }
}

mod shape {
    use super::grid::{GridType, GridWrapper};
    use super::DType;
    use crate::shapes::regular_sphere;
    use crate::types::RealScalar;
    use std::ffi::c_void;

    unsafe fn shapes_regular_sphere<T: RealScalar>(level: u32, dtype: DType) -> *const GridWrapper {
        let grid = Box::into_raw(Box::new(regular_sphere::<T>(level))) as *const c_void;
        Box::into_raw(Box::new(GridWrapper {
            grid,
            dtype,
            gtype: GridType::SerialSingleElementGrid,
        }))
    }
    #[no_mangle]
    pub unsafe extern "C" fn shapes_regular_sphere_f64(level: u32) -> *const GridWrapper {
        shapes_regular_sphere::<f64>(level, DType::F64)
    }
    #[no_mangle]
    pub unsafe extern "C" fn shapes_regular_sphere_f32(level: u32) -> *const GridWrapper {
        shapes_regular_sphere::<f32>(level, DType::F32)
    }
}
