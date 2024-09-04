//! Gmsh I/O

use crate::traits::{Entity, Geometry, GmshExport, Grid, Point, Topology};
use itertools::izip;
use ndelement::types::ReferenceCellType;
use num::Zero;

fn get_permutation_to_gmsh(cell_type: ReferenceCellType, degree: usize) -> Vec<usize> {
    match cell_type {
        ReferenceCellType::Triangle => match degree {
            1 => vec![0, 1, 2],
            2 => vec![0, 1, 2, 5, 3, 4],
            3 => vec![0, 1, 2, 7, 8, 3, 4, 6, 5, 9],
            4 => vec![0, 1, 2, 9, 10, 11, 3, 4, 5, 8, 7, 6, 12, 13, 14],
            5 => vec![
                0, 1, 2, 11, 12, 13, 14, 3, 4, 5, 6, 10, 9, 8, 7, 15, 16, 17, 18, 19, 20,
            ],
            _ => {
                panic!("Unsupported degree");
            }
        },
        ReferenceCellType::Quadrilateral => match degree {
            1 => vec![0, 1, 3, 2],
            2 => vec![0, 1, 3, 2, 4, 6, 7, 5, 8],
            _ => {
                panic!("Unsupported degree");
            }
        },
        ReferenceCellType::Tetrahedron => match degree {
            1 => vec![0, 1, 2, 3],
            _ => {
                panic!("Unsupported degree");
            }
        },
        ReferenceCellType::Hexahedron => match degree {
            1 => vec![0, 1, 3, 2, 4, 5, 7, 6],
            _ => {
                panic!("Unsupported degree");
            }
        },
        _ => {
            panic!("Unsupported cell type.");
        }
    }
}

fn get_gmsh_cell(cell_type: ReferenceCellType, degree: usize) -> usize {
    match cell_type {
        ReferenceCellType::Triangle => match degree {
            1 => 2,
            2 => 9,
            3 => 21,
            4 => 23,
            5 => 25,
            _ => {
                panic!("Unsupported degree");
            }
        },
        ReferenceCellType::Quadrilateral => match degree {
            1 => 3,
            2 => 10,
            _ => {
                panic!("Unsupported degree");
            }
        },
        ReferenceCellType::Tetrahedron => match degree {
            1 => 4,
            _ => {
                panic!("Unsupported degree");
            }
        },
        ReferenceCellType::Hexahedron => match degree {
            1 => 5,
            _ => {
                panic!("Unsupported degree");
            }
        },
        _ => {
            panic!("Unsupported cell type.");
        }
    }
}

impl<G: Grid<EntityDescriptor = ReferenceCellType>> GmshExport for G {
    fn to_gmsh_string(&self) -> String {
        let tdim = self.topology_dim();
        let gdim = self.geometry_dim();
        let node_count = self.entity_count(ReferenceCellType::Point);

        let mut gmsh_s = String::from("");
        gmsh_s.push_str("$MeshFormat\n");
        gmsh_s.push_str("4.1 0 8\n");
        gmsh_s.push_str("$EndMeshFormat\n");
        gmsh_s.push_str("$Nodes\n");
        gmsh_s.push_str(&format!("1 {node_count} 1 {node_count}\n"));
        gmsh_s.push_str(&format!("{tdim} 1 0 {node_count}\n"));
        for i in 0..node_count {
            gmsh_s.push_str(&format!("{}\n", i + 1));
        }
        let mut coords = vec![G::T::zero(); gdim];
        for node in self.entity_iter(0) {
            node.geometry().points().next().unwrap().coords(&mut coords);
            println!("{coords:?}");
            for (n, c) in coords.iter().enumerate() {
                if n != 0 {
                    gmsh_s.push(' ');
                }
                gmsh_s.push_str(&format!("{c}"));
            }
            gmsh_s.push('\n');
        }
        gmsh_s.push_str("$EndNodes\n");
        gmsh_s.push_str("$Elements\n");

        let cell_count = self
            .entity_types(tdim)
            .iter()
            .map(|t| self.entity_count(*t))
            .sum::<usize>();

        let mut elements = vec![];
        let mut cells_by_element = vec![];
        for (i, cell) in self.entity_iter(tdim).enumerate() {
            let element = (cell.entity_type(), cell.geometry().degree());
            if !elements.contains(&element) {
                elements.push(element);
                cells_by_element.push(vec![]);
            }
            cells_by_element[elements.iter().position(|i| *i == element).unwrap()].push(i);
        }

        gmsh_s.push_str(&format!("{} {cell_count} 1 {cell_count}\n", elements.len()));

        for ((cell_type, degree), cells) in izip!(elements, cells_by_element) {
            let gmsh_perm = get_permutation_to_gmsh(cell_type, degree);

            gmsh_s.push_str(&format!(
                "{tdim} 1 {} {}\n",
                get_gmsh_cell(cell_type, degree),
                cells.len()
            ));

            for (i, index) in cells.iter().enumerate() {
                gmsh_s.push_str(&format!("{}", i + 1));
                let entity = self.entity(tdim, *index).unwrap();
                let topology = entity.topology();
                for j in &gmsh_perm {
                    gmsh_s.push_str(&format!(" {}", topology.sub_entity(0, *j) + 1));
                }
                gmsh_s.push('\n');
            }
        }

        gmsh_s.push_str("$EndElements\n");

        gmsh_s
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{shapes::regular_sphere, traits::Builder, SingleElementGridBuilder};

    #[test]
    fn test_regular_sphere_gmsh_io() {
        let g = regular_sphere::<f64>(2);
        g.export_as_gmsh("_test_io_sphere.msh");
    }

    #[test]
    fn test_gmsh_output_quads() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));
        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0, 0.0]);
        b.add_point(2, &[0.0, 1.0, 0.0]);
        b.add_point(3, &[1.0, 1.0, 0.0]);
        b.add_point(4, &[0.0, 0.0, 1.0]);
        b.add_point(5, &[1.0, 0.0, 1.0]);
        b.add_point(6, &[0.0, 1.0, 1.0]);
        b.add_point(7, &[1.0, 1.0, 1.0]);
        b.add_cell(0, &[0, 2, 1, 3]);
        b.add_cell(1, &[0, 1, 4, 5]);
        b.add_cell(2, &[0, 4, 2, 6]);
        b.add_cell(3, &[1, 3, 5, 7]);
        b.add_cell(4, &[2, 6, 3, 7]);
        b.add_cell(5, &[4, 5, 6, 7]);
        let g = b.create_grid();
        g.export_as_gmsh("_test_io_cube.msh");
    }

    #[test]
    fn test_gmsh_output_tetrahedra() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Tetrahedron, 1));
        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0, 0.0]);
        b.add_point(2, &[0.0, 1.0, 0.0]);
        b.add_point(3, &[1.0, 1.0, 0.0]);
        b.add_point(4, &[0.0, 0.0, 1.0]);
        b.add_point(5, &[1.0, 0.0, 1.0]);
        b.add_point(6, &[0.0, 1.0, 1.0]);
        b.add_point(7, &[1.0, 1.0, 1.0]);
        b.add_cell(0, &[0, 1, 5, 7]);
        b.add_cell(1, &[0, 2, 6, 7]);
        b.add_cell(2, &[0, 4, 5, 7]);
        b.add_cell(3, &[0, 1, 3, 7]);
        b.add_cell(4, &[0, 2, 3, 7]);
        b.add_cell(5, &[0, 4, 6, 7]);
        let g = b.create_grid();
        g.export_as_gmsh("_test_io_tetrahedra.msh");
    }
    #[test]
    fn test_gmsh_output_hexahedra() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Hexahedron, 1));
        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0, 0.0]);
        b.add_point(2, &[0.0, 1.0, 0.0]);
        b.add_point(3, &[1.0, 1.0, 0.0]);
        b.add_point(4, &[0.0, 0.0, 1.0]);
        b.add_point(5, &[1.0, 0.0, 1.0]);
        b.add_point(6, &[0.0, 1.0, 1.0]);
        b.add_point(7, &[1.0, 1.0, 1.0]);
        b.add_point(8, &[0.0, 0.0, 2.0]);
        b.add_point(9, &[1.0, 0.0, 2.0]);
        b.add_point(10, &[0.0, 1.0, 2.0]);
        b.add_point(11, &[1.0, 1.0, 1.5]);
        b.add_cell(1, &[0, 1, 2, 3, 4, 5, 6, 7]);
        b.add_cell(2, &[4, 5, 6, 7, 8, 9, 10, 11]);
        let g = b.create_grid();
        g.export_as_gmsh("_test_io_hexahedra.msh");
    }
}
