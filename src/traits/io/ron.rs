//! RON I/O
use crate::traits::Grid;
#[cfg(feature = "mpi")]
use crate::traits::ParallelGrid;
#[cfg(feature = "mpi")]
use mpi::traits::Communicator;
use std::fs;

pub trait ConvertToSerializable {
    //! Convert to/from a RON string
    type SerializableType: serde::Serialize;
    /// Convert to ron
    fn to_serializable(&self) -> Self::SerializableType;
    /// Convert from ron
    fn from_serializable(ron: Self::SerializableType) -> Self;
}

pub trait RONExport: Grid {
    //! Grid export for RON

    /// Generate the RON string for a grid
    fn to_ron_string(&self) -> String;

    /// Export as RON
    fn export_as_ron(&self, filename: &str) {
        let ron_s = self.to_ron_string();
        fs::write(filename, ron_s).expect("Unable to write file");
    }
}

pub trait RONImport: Sized + Grid {
    //! Grid import for RON

    /// Generate the RON string for a grid
    fn from_ron_string(s: String) -> Self;

    /// Export as RON
    fn import_from_ron(filename: &str) -> Self {
        let content = fs::read_to_string(filename).expect("Unable to read file");
        Self::from_ron_string(content)
    }
}

#[derive(Debug, serde::Serialize, serde::Deserialize)]
/// Summary I/O data for a parallel grid
pub struct ParallelGridSummaryData {
    mpi_ranks: i32,
}

#[cfg(feature = "mpi")]
pub trait RONExportParallel<'a, C: Communicator + 'a>: ParallelGrid<C>
where
    Self::LocalGrid<'a>: RONExport,
    Self: 'a,
{
    //! Parallel grid export for RON

    /// Export as RON
    fn export_as_ron(&'a self, filename: &str) {
        let parts = filename.split('.').collect::<Vec<_>>();
        assert!(parts.len() > 1);
        let sub_filename = format!(
            "{}.{}.{}",
            parts[0..parts.len() - 1].join("."),
            self.comm().rank(),
            parts[parts.len() - 1]
        );

        self.local_grid().export_as_ron(&sub_filename);
        if self.comm().rank() == 0 {
            let grid_data = ParallelGridSummaryData {
                mpi_ranks: self.comm().size(),
            };
            fs::write(filename, ron::to_string(&grid_data).unwrap()).expect("Unable to write file");
        }
    }
}

#[cfg(feature = "mpi")]
pub trait RONImportParallel<'a, C: Communicator + 'a>: Sized + ParallelGrid<C>
where
    Self::LocalGrid<'a>: RONImport,
    Self: 'a,
{
    //! Parallel grid import for RON

    /// Create parallel grid from comm and local_grid
    fn create_from_ron_info(comm: &'a C, local_grid: Self::LocalGrid<'a>) -> Self;

    /// Export as RON
    fn import_from_ron(comm: &'a C, filename: &str) -> Self {
        let parts = filename.split('.').collect::<Vec<_>>();
        assert!(parts.len() > 1);
        let sub_filename = format!(
            "{}.{}.{}",
            parts[0..parts.len() - 1].join("."),
            comm.rank(),
            parts[parts.len() - 1]
        );

        let content = fs::read_to_string(filename).expect("Unable to read file");
        let summary: ParallelGridSummaryData = ron::from_str(&content).unwrap();

        if summary.mpi_ranks != comm.size() {
            panic!("Incorrect number of MPI ranks");
        }

        let local_grid = Self::LocalGrid::import_from_ron(&sub_filename);
        Self::create_from_ron_info(comm, local_grid)
    }
}
