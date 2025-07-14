//! Gmsh I/O
use std::fs;
use crate::traits::{Grid, Builder};

pub trait GmshExport: Grid {
    //! Grid export for Gmsh

    /// Generate the Gmsh string for a grid
    fn to_gmsh_string(&self) -> String;

    /// Export as Gmsh
    fn export_as_gmsh(&self, filename: &str) {
        let gmsh_s = self.to_gmsh_string();
        fs::write(filename, gmsh_s).expect("Unable to write file");
    }
}

pub trait GmshImport: Builder {
    //! Grid import for Gmsh

    /// Generate grid from a Gmsh string
    fn import_from_gmsh_string(&mut self, s: String);

    /// Export as Gmsh
    fn import_from_gmsh(&mut self, filename: &str) {
        let content = fs::read_to_string(filename).expect("Unable to read file");
        self.import_from_gmsh_string(content);
    }
}
