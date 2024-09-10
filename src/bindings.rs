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
    use ndelement::ciarlet::CiarletElement;
    use std::ffi::c_void;

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
