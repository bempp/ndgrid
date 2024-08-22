//! I/O as RON
use super::super::SingleElementGrid;
use crate::{
    traits::{ConvertToSerializable, RONExport, RONImport},
    types::RealScalar,
};
use ndelement::ciarlet::CiarletElement;

impl<T: RealScalar> RONExport for SingleElementGrid<T, CiarletElement<T>> {
    fn to_ron_string(&self) -> String {
        ron::to_string(&self.to_serializable()).unwrap()
    }
}

impl<T: RealScalar + for<'de> serde::Deserialize<'de>> RONImport
    for SingleElementGrid<T, CiarletElement<T>>
{
    fn from_ron_string(s: String) -> Self {
        Self::from_serializable(ron::from_str(&s).unwrap())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::shapes::regular_sphere;
    use crate::traits::Grid;
    use ndelement::types::ReferenceCellType;

    #[test]
    fn test_ron_export_and_import() {
        let g = regular_sphere::<f64>(1);
        let n = g.entity_count(ReferenceCellType::Interval);
        g.export_as_ron("_test_export.ron".to_string());

        let g2 = SingleElementGrid::<f64, CiarletElement<f64>>::import_from_ron(
            "_test_export.ron".to_string(),
        );

        assert_eq!(g2.entity_count(ReferenceCellType::Interval), n);
    }
}
