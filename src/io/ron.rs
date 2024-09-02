//! RON I/O
use crate::traits::{ConvertToSerializable, Grid, RONExport, RONImport};
#[cfg(feature = "mpi")]
use crate::traits::{ParallelGrid, RONExportParallel};
#[cfg(feature = "mpi")]
use mpi::traits::Communicator;

impl<G: Grid + ConvertToSerializable> RONExport for G {
    fn to_ron_string(&self) -> String {
        ron::to_string(&self.to_serializable()).unwrap()
    }
}

impl<G: Grid + ConvertToSerializable> RONImport for G
where
    for<'a> G::SerializableType: serde::Deserialize<'a>,
{
    fn from_ron_string(s: String) -> Self {
        Self::from_serializable(ron::from_str(&s).unwrap())
    }
}

#[cfg(feature = "mpi")]
impl<'a, C: Communicator + 'a, G: ParallelGrid<C>> RONExportParallel<'a, C> for G
where
    Self::LocalGrid<'a>: RONExport,
    Self: 'a,
{
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::shapes::regular_sphere;
    use crate::traits::Grid;
    use crate::SingleElementGrid;
    use ndelement::ciarlet::CiarletElement;
    use ndelement::types::ReferenceCellType;

    #[test]
    fn test_ron_export_and_import() {
        let g = regular_sphere::<f64>(1);
        let n = g.entity_count(ReferenceCellType::Interval);
        g.export_as_ron("_test_export.ron");

        let g2 = SingleElementGrid::<f64, CiarletElement<f64>>::import_from_ron("_test_export.ron");

        assert_eq!(g2.entity_count(ReferenceCellType::Interval), n);
    }
}
