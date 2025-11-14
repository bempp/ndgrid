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
impl<'a, C: Communicator + 'a, G: ParallelGrid<C = C>> RONExportParallel<'a, C> for G
where
    Self::LocalGrid: RONExport,
    Self: 'a,
{
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        MixedGrid, MixedGridBuilder, SingleElementGrid,
        shapes::regular_sphere,
        traits::{Builder, Grid},
    };
    use ndelement::{ciarlet::CiarletElement, map::IdentityMap, types::ReferenceCellType};

    #[test]
    fn test_ron_export_and_import() {
        let g = regular_sphere::<f64>(1);
        let n = g.entity_count(ReferenceCellType::Interval);
        g.export_as_ron("_test_export.ron");

        let g2 = SingleElementGrid::<f64, CiarletElement<f64, IdentityMap>>::import_from_ron(
            "_test_export.ron",
        );

        assert_eq!(g2.entity_count(ReferenceCellType::Interval), n);
    }

    #[test]
    fn test_ron_export_and_import_mixed() {
        let mut b = MixedGridBuilder::<f64>::new(2);
        b.add_point(0, &[0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0]);
        b.add_point(2, &[0.0, 1.0]);
        b.add_point(3, &[1.0, 1.0]);
        b.add_point(4, &[2.0, 0.5]);
        b.add_point(5, &[1.6, 0.9]);
        b.add_point(6, &[1.0, 0.5]);
        b.add_point(7, &[1.6, 0.1]);

        b.add_cell(0, (ReferenceCellType::Quadrilateral, 1, &[0, 1, 2, 3]));
        b.add_cell(1, (ReferenceCellType::Triangle, 2, &[1, 4, 3, 5, 6, 7]));
        let g = b.create_grid();

        let n = g.entity_count(ReferenceCellType::Interval);
        g.export_as_ron("_test_export_mixed.ron");

        let g2 = MixedGrid::<f64, CiarletElement<f64, IdentityMap>>::import_from_ron(
            "_test_export_mixed.ron",
        );

        assert_eq!(g2.entity_count(ReferenceCellType::Interval), n);
    }
}
