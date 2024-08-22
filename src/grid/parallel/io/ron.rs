//! I/O to/from RON
use super::super::ParallelGrid;
use crate::{
    traits::{RONExport, RONImport, ConvertToSerializable},
    types::RealScalar,
};
use ndelement::ciarlet::CiarletElement;

impl<'a, C: Communicator, G: Grid + Sync + ConvertToSerializable> RONExport for ParallelGrid<'a, C, G> {
    fn to_ron_string(&self) -> String {
        ron::to_string(&self.local_grid().to_serializable()).unwrap()
    }
}
