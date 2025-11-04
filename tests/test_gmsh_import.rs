use ndelement::types::ReferenceCellType;
use ndgrid::{
    SingleElementGridBuilder,
    traits::{Builder, GmshImport},
};
use std::path::PathBuf;

fn relative_file(filename: &str) -> String {
    let file = PathBuf::from(file!());
    let dir = file.parent().unwrap();
    format!("{}/{filename}", dir.display())
}

#[test]
fn test_gmsh_import_v1() {
    let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));
    b.import_from_gmsh(&relative_file("test_mesh_v1.msh"));
    let _g = b.create_grid();
}

#[test]
fn test_gmsh_import_ascii_v2() {
    let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));
    b.import_from_gmsh(&relative_file("test_mesh_ascii_v2.msh"));
    let _g = b.create_grid();
}

#[test]
fn test_gmsh_import_binary_v2() {
    let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));
    b.import_from_gmsh(&relative_file("test_mesh_binary_v2.msh"));
    let _g = b.create_grid();
}

#[test]
fn test_gmsh_import_ascii_v4() {
    let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));
    b.import_from_gmsh(&relative_file("test_mesh_ascii_v4.msh"));
    let _g = b.create_grid();
}

#[test]
fn test_gmsh_import_binary_v4() {
    let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));
    b.import_from_gmsh(&relative_file("test_mesh_binary_v4.msh"));
    let _g = b.create_grid();
}
