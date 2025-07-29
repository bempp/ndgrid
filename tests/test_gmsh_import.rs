use ndelement::types::ReferenceCellType;
use ndgrid::{
    traits::{Builder, GmshImport},
    SingleElementGridBuilder,
};
use std::path::PathBuf;

fn relative_file(filename: &str) -> String {
    let file = PathBuf::from(file!());
    let dir = file.parent().unwrap();
    format!("{}/{filename}", dir.display())
}

#[test]
fn test_gmsh_import() {
    let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));
    b.import_from_gmsh(&relative_file("test_mesh.msh"));
    let _g = b.create_grid();
}
