//! RON I/O
use std::fs;

pub trait RONExport {
    //! Grid export for RON

    /// Generate the RON string for a grid
    fn to_ron_string(&self) -> String;

    /// Export as RON
    fn export_as_ron(&self, filename: String) {
        let ron_s = self.to_ron_string();
        fs::write(filename, ron_s).expect("Unable to write file");
    }
}

pub trait RONImport: Sized {
    //! Grid import for RON

    /// Generate the RON string for a grid
    fn from_ron_string(s: String) -> Self;

    /// Export as RON
    fn import_from_ron(filename: String) -> Self {
        let content = fs::read_to_string(filename).expect("Unable to read file");
        Self::from_ron_string(content)
    }
}

#[cfg(feature = "serde")]
pub trait ConvertToSerializable {
    //! Convert to/from a RON string
    type SerializableType: serde::Serialize;
    /// Convert to ron
    fn to_serializable(&self) -> Self::SerializableType;
    /// Convert from ron
    fn from_serializable(ron: Self::SerializableType) -> Self;
}
