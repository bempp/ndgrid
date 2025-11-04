//! Gmsh I/O

use crate::traits::{Builder, Entity, Geometry, GmshExport, GmshImport, Grid, Point};
use itertools::izip;
use ndelement::types::ReferenceCellType;
use num::Zero;
use num::traits::FromBytes;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::str::FromStr;

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

fn interpret_gmsh_cell(gmsh_cell: usize) -> (ReferenceCellType, usize) {
    match gmsh_cell {
        2 => (ReferenceCellType::Triangle, 1),
        9 => (ReferenceCellType::Triangle, 2),
        21 => (ReferenceCellType::Triangle, 3),
        23 => (ReferenceCellType::Triangle, 4),
        25 => (ReferenceCellType::Triangle, 5),
        3 => (ReferenceCellType::Quadrilateral, 1),
        10 => (ReferenceCellType::Quadrilateral, 2),
        4 => (ReferenceCellType::Tetrahedron, 1),
        5 => (ReferenceCellType::Hexahedron, 1),
        _ => {
            panic!("Unsupported cell type.");
        }
    }
}

/// Parse gmsh binary data
fn parse<T: FromBytes>(data: &[u8], gmsh_data_size: usize, is_le: bool) -> T
where
    <T as FromBytes>::Bytes: Sized,
{
    let mut buf: T::Bytes = unsafe { std::mem::zeroed() };
    let target_data_size = std::mem::size_of::<T::Bytes>();
    let buf_bytes: &mut [u8] = unsafe {
        std::slice::from_raw_parts_mut((&mut buf as *mut _) as *mut u8, target_data_size)
    };

    if is_le {
        buf_bytes[..gmsh_data_size].copy_from_slice(data);
        T::from_le_bytes(&buf)
    } else {
        buf_bytes[target_data_size - gmsh_data_size..].copy_from_slice(data);
        T::from_be_bytes(&buf)
    }
}

impl<G: Grid<EntityDescriptor = ReferenceCellType>> GmshExport for G {
    fn to_gmsh_string(&self) -> String {
        let tdim = self.topology_dim();
        let gdim = self.geometry_dim();

        let mut points = HashMap::new();
        for cell_type in self.cell_types() {
            for cell in self.entity_iter(*cell_type) {
                for point in cell.geometry().points() {
                    let mut p = vec![G::T::zero(); gdim];
                    point.coords(&mut p);
                    points.insert(point.index(), p);
                }
            }
        }
        let mut points = points.iter().collect::<Vec<_>>();
        points.sort_by(|i, j| i.0.cmp(j.0));
        let node_count = points.len();

        let mut gmsh_s = String::from("");
        gmsh_s.push_str("$MeshFormat\n");
        gmsh_s.push_str("4.1 0 8\n");
        gmsh_s.push_str("$EndMeshFormat\n");
        gmsh_s.push_str("$Nodes\n");
        gmsh_s.push_str(&format!("1 {node_count} 1 {node_count}\n"));
        gmsh_s.push_str(&format!("{tdim} 1 0 {node_count}\n"));
        for (i, _) in &points {
            gmsh_s.push_str(&format!("{}\n", *i + 1));
        }
        for (_, coords) in &points {
            for n in 0..3 {
                if n != 0 {
                    gmsh_s.push(' ');
                }
                gmsh_s.push_str(&format!(
                    "{}",
                    if let Some(c) = coords.get(n) {
                        *c
                    } else {
                        G::T::zero()
                    }
                ));
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
        for t in self.cell_types() {
            for (index, cell) in self.entity_iter(*t).enumerate() {
                let element = (cell.entity_type(), cell.geometry().degree());
                if !elements.contains(&element) {
                    elements.push(element);
                    cells_by_element.push(vec![]);
                }
                cells_by_element[elements.iter().position(|i| *i == element).unwrap()].push(index);
            }
        }

        gmsh_s.push_str(&format!("{} {cell_count} 1 {cell_count}\n", elements.len()));

        let mut next_index = 1;
        for ((cell_type, degree), cells) in izip!(elements, cells_by_element) {
            let gmsh_perm = get_permutation_to_gmsh(cell_type, degree);

            gmsh_s.push_str(&format!(
                "{tdim} 1 {} {}\n",
                get_gmsh_cell(cell_type, degree),
                cells.len()
            ));

            for index in cells.iter() {
                gmsh_s.push_str(&format!("{}", {
                    let idx = next_index;
                    next_index += 1;
                    idx
                },));
                let entity = self.entity(cell_type, *index).unwrap();
                let point_indices = entity
                    .geometry()
                    .points()
                    .map(|i| i.index())
                    .collect::<Vec<_>>();
                for j in &gmsh_perm {
                    gmsh_s.push_str(&format!(" {}", point_indices[*j] + 1));
                }
                gmsh_s.push('\n');
            }
        }

        gmsh_s.push_str("$EndElements\n");

        gmsh_s
    }
}

/// Get a section from a gmsh string
fn gmsh_section(s: &str, section: &str) -> String {
    let a = s.split(&format!("${section}\n")).collect::<Vec<_>>();
    if a.len() <= 1 {
        panic!("Section not found: {section}");
    }
    String::from(a[1].split(&format!("\n$End{section}")).collect::<Vec<_>>()[0])
}

impl<T, B> GmshImport for B
where
    T: FromStr + FromBytes + std::fmt::Debug,
    <T as FromBytes>::Bytes: Sized,
    B: Builder<T = T, EntityDescriptor = ReferenceCellType>,
{
    fn import_from_gmsh_v1(&mut self, s: String) {
        // Load nodes
        let nodes = gmsh_section(&s, "Nodes");
        let mut nodes = nodes.lines();

        let Some(num_nodes) = nodes.next() else {
            panic!("Could not read num nodes");
        };
        let num_nodes = num_nodes
            .parse::<usize>()
            .expect("Could not parse num nodes");

        for _ in 0..num_nodes {
            let Some(line) = nodes.next() else {
                panic!("Unable to read node line");
            };
            let line = line.split(" ").collect::<Vec<&str>>();

            let (tag, coords) = line.split_at(1);
            let tag = tag[0].parse::<usize>().expect("Could not parse node tag");
            let coords = coords
                .iter()
                .map(|c| {
                    if let Ok(d) = T::from_str(c) {
                        d
                    } else {
                        panic!("Could not parse coords")
                    }
                })
                .collect::<Vec<_>>();

            println!("{coords:?}");
            self.add_point(tag, &coords[..self.gdim()]);
            println!("ADDED");
        }

        // Load elements
        let elements = gmsh_section(&s, "Elements");
        let mut elements = elements.lines();

        let Some(num_elements) = elements.next() else {
            panic!("Could not read num nodes");
        };
        let num_elements = num_elements
            .parse::<usize>()
            .expect("Could not parse num nodes");

        for _ in 0..num_elements {
            let Some(line) = elements.next() else {
                panic!("Unable to read element line");
            };
            let line = line.split(" ").collect::<Vec<&str>>();

            let (a, node_tags) = line.split_at(5);
            let [tag, element_type, ..] = a
                .iter()
                .map(|i| {
                    i.parse::<usize>()
                        .expect("Could not parse element tag and type")
                })
                .collect::<Vec<_>>()[..]
            else {
                panic!("Unrecognised format");
            };

            let node_tags = node_tags
                .iter()
                .map(|i| {
                    i.parse::<usize>()
                        .expect("Could not parse nodes for element")
                })
                .collect::<Vec<_>>();
            let (cell_type, degree) = interpret_gmsh_cell(element_type);
            let gmsh_perm = get_permutation_to_gmsh(cell_type, degree);

            let mut cell = vec![0; node_tags.len()];
            for (i, j) in gmsh_perm.iter().enumerate() {
                cell[*j] = node_tags[i];
            }

            self.add_cell_from_nodes_and_type(tag, &cell, cell_type, degree);
        }
    }

    fn import_from_gmsh_string_v2(&mut self, s: String) {
        // Load nodes
        let nodes = gmsh_section(&s, "Nodes");
        let mut nodes = nodes.lines();

        let Some(num_nodes) = nodes.next() else {
            panic!("Could not read num nodes");
        };
        let num_nodes = num_nodes
            .parse::<usize>()
            .expect("Could not parse num nodes");

        for _ in 0..num_nodes {
            let Some(line) = nodes.next() else {
                panic!("Unable to read node line");
            };
            let line = line.split(" ").collect::<Vec<&str>>();

            let (tag, coords) = line.split_at(1);
            let tag = tag[0].parse::<usize>().expect("Could not parse node tag");
            let coords = coords
                .iter()
                .map(|c| {
                    if let Ok(d) = T::from_str(c) {
                        d
                    } else {
                        panic!("Could not parse coords")
                    }
                })
                .collect::<Vec<_>>();

            self.add_point(tag, &coords[..self.gdim()]);
        }

        // Load elements
        let elements = gmsh_section(&s, "Elements");
        let mut elements = elements.lines();

        let Some(num_elements) = elements.next() else {
            panic!("Could not read num nodes");
        };
        let num_elements = num_elements
            .parse::<usize>()
            .expect("Could not parse num nodes");

        for _ in 0..num_elements {
            let Some(line) = elements.next() else {
                panic!("Unable to read element line");
            };
            let line = line.split(" ").collect::<Vec<&str>>();

            let (a, rem_line) = line.split_at(3);
            let [tag, element_type, num_tags] = a
                .iter()
                .map(|i| {
                    i.parse::<usize>()
                        .expect("Could not parse element tag and type")
                })
                .collect::<Vec<_>>()[..]
            else {
                panic!("Unrecognised format");
            };

            let (_, node_tags) = rem_line.split_at(num_tags);
            let node_tags = node_tags
                .iter()
                .map(|i| {
                    i.parse::<usize>()
                        .expect("Could not parse nodes for element")
                })
                .collect::<Vec<_>>();
            let (cell_type, degree) = interpret_gmsh_cell(element_type);
            let gmsh_perm = get_permutation_to_gmsh(cell_type, degree);

            let mut cell = vec![0; node_tags.len()];
            for (i, j) in gmsh_perm.iter().enumerate() {
                cell[*j] = node_tags[i];
            }

            self.add_cell_from_nodes_and_type(tag, &cell, cell_type, degree);
        }
    }

    fn import_from_gmsh_binary_v2(
        &mut self,
        mut reader: BufReader<File>,
        data_size: usize,
        is_le: bool,
    ) {
        let mut line = String::new();
        let mut buf = Vec::new();

        macro_rules! read_exact {
            ($size: expr, $msg: expr) => {{
                buf.resize($size, 0);
                reader.read_exact(&mut buf).expect($msg);
            }};
        }

        const GMSH_INT_SIZE: usize = 4;

        loop {
            let Ok(num_bytes) = reader.read_line(&mut line) else {
                continue;
            };

            // EOF reached.
            if num_bytes == 0 {
                break;
            }

            match line.trim() {
                // Load all nodes.
                "$Nodes" => {
                    line.clear();
                    let Ok(_) = reader.read_line(&mut line) else {
                        panic!("Unable to read num nodes");
                    };
                    let num_nodes = line
                        .trim()
                        .parse::<usize>()
                        .expect("Could not parse num nodes");

                    for _ in 0..num_nodes {
                        read_exact!(GMSH_INT_SIZE, "Unable to read node tag");
                        let tag = parse::<usize>(&buf, GMSH_INT_SIZE, is_le);

                        read_exact!(3 * data_size, "Unable to read node coords");
                        let coords = buf
                            .chunks(data_size)
                            .map(|i| parse::<T>(i, data_size, is_le))
                            .collect::<Vec<_>>();

                        self.add_point(tag, &coords[..self.gdim()]);
                    }

                    line.clear();
                }
                // Load all elements.
                "$Elements" => {
                    line.clear();
                    let Ok(_) = reader.read_line(&mut line) else {
                        panic!("Unable to read num elements")
                    };
                    let num_elements = line
                        .trim()
                        .parse::<usize>()
                        .expect("Could not parse num elements");

                    for _ in 0..num_elements {
                        read_exact!(3 * GMSH_INT_SIZE, "Unable to element tag and type");
                        let [elm_type, _num_elm_follow, num_tags] = buf
                            .chunks(GMSH_INT_SIZE)
                            .map(|i| parse::<usize>(i, GMSH_INT_SIZE, is_le))
                            .collect::<Vec<_>>()[..]
                        else {
                            panic!("Could not parse element tag and type")
                        };

                        read_exact!(GMSH_INT_SIZE, "Unable to read num tags");
                        let tag = parse::<usize>(&buf, 4, is_le);

                        // Skip tags
                        read_exact!(num_tags * GMSH_INT_SIZE, "Unable to read element tags");

                        let (cell_type, degree) = interpret_gmsh_cell(elm_type);
                        let gmsh_perm = get_permutation_to_gmsh(cell_type, degree);

                        read_exact!(
                            gmsh_perm.len() * GMSH_INT_SIZE,
                            "Unable to read element node tags"
                        );
                        let node_tags = buf
                            .chunks(GMSH_INT_SIZE)
                            .map(|i| parse::<usize>(i, GMSH_INT_SIZE, is_le))
                            .collect::<Vec<_>>();

                        let mut cell = vec![0; node_tags.len()];
                        for (i, j) in gmsh_perm.iter().enumerate() {
                            cell[*j] = node_tags[i];
                        }

                        self.add_cell_from_nodes_and_type(tag, &cell, cell_type, degree);
                    }
                    line.clear();
                }
                _ => {
                    line.clear();
                }
            }
        }
    }

    fn import_from_gmsh_string_v4(&mut self, s: String) {
        // Load nodes
        let nodes = gmsh_section(&s, "Nodes");
        let nodes = nodes.lines().collect::<Vec<_>>();

        let [num_entity_blocks, _num_nodes, _min_node_tag, _max_node_tag] = nodes[0]
            .trim()
            .split(" ")
            .map(|i| i.parse::<usize>().unwrap())
            .collect::<Vec<_>>()[..]
        else {
            panic!("Unrecognised gmsh format for node blocks");
        };

        let mut line_n = 1;
        for _ in 0..num_entity_blocks {
            let [_entity_dim, _entity_tag, parametric, num_nodes_in_block] = nodes[line_n]
                .trim()
                .split(" ")
                .map(|i| i.parse::<usize>().unwrap())
                .collect::<Vec<_>>()[..]
            else {
                panic!("Unrecognised gmsh format for nodes");
            };
            if parametric == 1 {
                unimplemented!("Parametric nodes currently not supported")
            }
            line_n += 1;
            let tags = &nodes[line_n..line_n + num_nodes_in_block];
            let coords = &nodes[line_n + num_nodes_in_block..line_n + 2 * num_nodes_in_block];
            for (t, c) in izip!(tags, coords) {
                let pt = c
                    .trim()
                    .split(" ")
                    .map(|i| {
                        if let Ok(j) = T::from_str(i) {
                            j
                        } else {
                            panic!("Could not parse coordinate");
                        }
                    })
                    .collect::<Vec<_>>();
                self.add_point(t.parse::<usize>().unwrap(), &pt[..self.gdim()]);
            }
            line_n += num_nodes_in_block * 2;
        }

        // Load elements
        let elements = gmsh_section(&s, "Elements");
        let elements = elements.lines().collect::<Vec<_>>();

        let [
            num_entity_blocks,
            _num_elements,
            _min_element_tag,
            _max_element_tag,
        ] = elements[0]
            .trim()
            .split(" ")
            .map(|i| i.parse::<usize>().unwrap())
            .collect::<Vec<_>>()[..]
        else {
            panic!("Unrecognised gmsh format");
        };

        let mut line_n = 1;
        for _ in 0..num_entity_blocks {
            let [
                _entity_dim,
                _entity_tag,
                element_type,
                num_elements_in_block,
            ] = elements[line_n]
                .trim()
                .split(" ")
                .map(|i| i.parse::<usize>().unwrap())
                .collect::<Vec<_>>()[..]
            else {
                panic!("Unrecognised gmsh format");
            };
            let (cell_type, degree) = interpret_gmsh_cell(element_type);
            let gmsh_perm = get_permutation_to_gmsh(cell_type, degree);

            line_n += 1;
            for line in &elements[line_n..line_n + num_elements_in_block] {
                let line = line
                    .trim()
                    .split(" ")
                    .map(|i| i.parse::<usize>().unwrap())
                    .collect::<Vec<_>>();
                let mut cell = vec![0; line.len() - 1];
                for (i, j) in gmsh_perm.iter().enumerate() {
                    cell[*j] = line[i + 1];
                }
                self.add_cell_from_nodes_and_type(line[0], &cell, cell_type, degree);
            }

            line_n += num_elements_in_block;
        }
    }

    fn import_from_gmsh_binary_v4(
        &mut self,
        mut reader: BufReader<File>,
        data_size: usize,
        is_le: bool,
    ) {
        let mut line = String::new();
        let mut buf = Vec::new();

        macro_rules! read_exact {
            ($size: expr, $msg: expr) => {{
                buf.resize($size, 0);
                reader.read_exact(&mut buf).expect($msg);
            }};
        }

        const GMSH_INT_SIZE: usize = 4;

        loop {
            let Ok(num_bytes) = reader.read_line(&mut line) else {
                continue;
            };

            // EOF reached.
            if num_bytes == 0 {
                break;
            }

            match line.trim() {
                // Load all nodes.
                "$Nodes" => {
                    read_exact!(4 * data_size, "Unable to read node section info");
                    let [num_entity_blocks, _num_nodes, _min_node_tag, _max_node_tag] = buf
                        .chunks(data_size)
                        .map(|i| parse::<usize>(i, data_size, is_le))
                        .collect::<Vec<_>>()[..]
                    else {
                        panic!("Could not parse node section info")
                    };

                    for _ in 0..num_entity_blocks {
                        read_exact!(3 * GMSH_INT_SIZE, "Unable to read node entity block info");
                        let [_entity_dim, _entity_tag, parametric] = buf
                            .chunks(GMSH_INT_SIZE)
                            .map(|i| parse::<usize>(i, GMSH_INT_SIZE, is_le))
                            .collect::<Vec<_>>()[..]
                        else {
                            panic!("Could not parse node entity block info")
                        };

                        if parametric == 1 {
                            unimplemented!("Parametric nodes currently not supported")
                        }

                        read_exact!(data_size, "Unable to read num nodes in block");
                        let num_nodes_in_block = parse::<usize>(&buf, data_size, is_le);

                        read_exact!(num_nodes_in_block * data_size, "Unable to read node tags");
                        let tags = buf
                            .chunks(data_size)
                            .map(|i| parse::<usize>(i, data_size, is_le))
                            .collect::<Vec<_>>();

                        read_exact!(
                            3 * num_nodes_in_block * data_size,
                            "Unable to read node coords"
                        );
                        let coords = buf
                            .chunks(3 * data_size)
                            .map(|i| {
                                i.chunks(data_size)
                                    .map(|j| parse::<T>(j, data_size, is_le))
                                    .collect::<Vec<_>>()
                            })
                            .collect::<Vec<_>>();

                        for (t, c) in izip!(tags, coords) {
                            self.add_point(t, &c[..self.gdim()]);
                        }
                    }

                    line.clear();
                }
                // Load all elements.
                "$Elements" => {
                    read_exact!(4 * data_size, "Unable to read element section info");
                    let [
                        num_entity_blocks,
                        _num_elements,
                        _min_element_tag,
                        _max_element_tag,
                    ] = buf
                        .chunks(data_size)
                        .map(|i| parse::<usize>(i, data_size, is_le))
                        .collect::<Vec<_>>()[..]
                    else {
                        panic!("Could not parse element section info")
                    };

                    for _ in 0..num_entity_blocks {
                        read_exact!(
                            3 * GMSH_INT_SIZE,
                            "Unable to read element entity block info"
                        );
                        let [_entity_dim, _entity_tag, element_type] = buf
                            .chunks(GMSH_INT_SIZE)
                            .map(|i| parse::<usize>(i, GMSH_INT_SIZE, is_le))
                            .collect::<Vec<_>>()[..]
                        else {
                            panic!("Could not parse element entity block info")
                        };

                        read_exact!(data_size, "Unable to read num elements in block");
                        let num_elements_in_block = parse::<usize>(&buf, data_size, is_le);

                        let (cell_type, degree) = interpret_gmsh_cell(element_type);
                        let gmsh_perm = get_permutation_to_gmsh(cell_type, degree);

                        for _ in 0..num_elements_in_block {
                            read_exact!(data_size, "Unable to read element tag");
                            let tag = parse::<usize>(&buf, data_size, is_le);

                            read_exact!(data_size * gmsh_perm.len(), "Unable to read node tags");
                            let node_tags = buf
                                .chunks(data_size)
                                .map(|i| parse::<usize>(i, data_size, is_le))
                                .collect::<Vec<_>>();

                            let mut cell = vec![0; node_tags.len()];
                            for (i, j) in gmsh_perm.iter().enumerate() {
                                cell[*j] = node_tags[i];
                            }

                            self.add_cell_from_nodes_and_type(tag, &cell, cell_type, degree);
                        }
                    }
                    line.clear();
                }
                _ => {
                    line.clear();
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        MixedGridBuilder, SingleElementGridBuilder,
        shapes::regular_sphere,
        traits::{Builder, Topology},
    };
    use approx::*;

    #[test]
    fn test_regular_sphere_gmsh_io() {
        let g = regular_sphere::<f64>(2);
        g.export_as_gmsh("_test_io_sphere.msh");
    }

    #[test]
    fn test_export_quads() {
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
    fn test_export_tetrahedra() {
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
    fn test_hexahedra() {
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

        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Hexahedron, 1));
        b.import_from_gmsh("_test_io_hexahedra.msh");
        let g2 = b.create_grid();

        let mut p1 = [0.0; 3];
        let mut p2 = [0.0; 3];
        for (v1, v2) in izip!(
            g.entity_iter(ReferenceCellType::Point),
            g2.entity_iter(ReferenceCellType::Point)
        ) {
            for (pt1, pt2) in izip!(v1.geometry().points(), v2.geometry().points()) {
                pt1.coords(&mut p1);
                pt2.coords(&mut p2);
                for (c1, c2) in izip!(&p1, &p2) {
                    assert_relative_eq!(c1, c2, epsilon = 1e-10);
                }
            }
        }
        for (h1, h2) in izip!(
            g.entity_iter(ReferenceCellType::Hexahedron),
            g2.entity_iter(ReferenceCellType::Hexahedron)
        ) {
            for (v1, v2) in izip!(
                h1.topology().sub_entity_iter(ReferenceCellType::Point),
                h2.topology().sub_entity_iter(ReferenceCellType::Point)
            ) {
                assert_eq!(v1, v2);
            }
        }
    }

    #[test]
    fn test_high_order_triangles() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 5));
        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[5.0, 0.0, 0.2]);
        b.add_point(2, &[0.0, 5.0, 0.4]);
        b.add_point(3, &[4.0, 1.0, -0.1]);
        b.add_point(4, &[3.0, 2.0, 0.2]);
        b.add_point(5, &[2.0, 3.0, -0.3]);
        b.add_point(6, &[1.0, 4.0, -0.5]);
        b.add_point(7, &[0.0, 1.0, 0.6]);
        b.add_point(8, &[0.0, 2.0, 0.2]);
        b.add_point(9, &[0.0, 3.0, 0.1]);
        b.add_point(10, &[0.0, 4.0, -0.2]);
        b.add_point(11, &[1.0, 0.0, -0.3]);
        b.add_point(12, &[2.0, 0.0, -0.4]);
        b.add_point(13, &[3.0, 0.0, -0.5]);
        b.add_point(14, &[4.0, 0.0, -0.2]);
        b.add_point(15, &[1.0, 1.0, 0.1]);
        b.add_point(16, &[2.0, 1.0, 0.1]);
        b.add_point(17, &[3.0, 1.0, 0.1]);
        b.add_point(18, &[2.0, 1.0, 0.2]);
        b.add_point(19, &[2.0, 2.0, 0.1]);
        b.add_point(20, &[3.0, 1.0, 0.1]);
        b.add_cell(
            1,
            &[
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
            ],
        );
        let g = b.create_grid();
        g.export_as_gmsh("_test_io_high_order_triangle.msh");

        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 5));
        b.import_from_gmsh("_test_io_high_order_triangle.msh");
        let g2 = b.create_grid();

        let mut p1 = [0.0; 3];
        let mut p2 = [0.0; 3];
        for (v1, v2) in izip!(
            g.entity_iter(ReferenceCellType::Point),
            g2.entity_iter(ReferenceCellType::Point)
        ) {
            v1.geometry().points().next().unwrap().coords(&mut p1);
            v2.geometry().points().next().unwrap().coords(&mut p2);
        }
        for (v1, v2) in izip!(
            g.entity_iter(ReferenceCellType::Point),
            g2.entity_iter(ReferenceCellType::Point)
        ) {
            v1.geometry().points().next().unwrap().coords(&mut p1);
            v2.geometry().points().next().unwrap().coords(&mut p2);
            for (c1, c2) in izip!(&p1, &p2) {
                assert_relative_eq!(c1, c2, epsilon = 1e-10);
            }
        }
        for (h1, h2) in izip!(
            g.entity_iter(ReferenceCellType::Triangle),
            g2.entity_iter(ReferenceCellType::Triangle)
        ) {
            for (v1, v2) in izip!(
                h1.topology().sub_entity_iter(ReferenceCellType::Point),
                h2.topology().sub_entity_iter(ReferenceCellType::Point)
            ) {
                assert_eq!(v1, v2);
            }
        }
    }

    #[test]
    fn test_export_mixed() {
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
        g.export_as_gmsh("_test_io_mixed.msh");

        let mut b = MixedGridBuilder::<f64>::new(2);
        b.import_from_gmsh("_test_io_mixed.msh");
        let _g = b.create_grid();
    }

    #[test]
    fn test_import_mixed() {
        let mut b = MixedGridBuilder::<f64>::new(2);
        b.import_from_gmsh("meshes/mixed.msh");
        let _g = b.create_grid();
    }
    #[test]
    fn test_import_mixed_triangle() {
        let mut b = MixedGridBuilder::<f64>::new(3);
        b.import_from_gmsh("meshes/sphere_triangle.msh");
        let _g = b.create_grid();
    }

    #[test]
    fn test_import_triangle() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));
        b.import_from_gmsh("meshes/sphere_triangle.msh");
        let _g = b.create_grid();
    }

    #[test]
    fn test_import_quadrilateral() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));
        b.import_from_gmsh("meshes/cube_quadrilateral.msh");
        let _g = b.create_grid();
    }

    #[test]
    fn test_import_tetrahedron() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Tetrahedron, 1));
        b.import_from_gmsh("meshes/cube_tetrahedron.msh");
        let _g = b.create_grid();
    }

    #[test]
    #[should_panic]
    fn test_import_wrong_cell() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));
        b.import_from_gmsh("meshes/cube_tetrahedron.msh");
        let _g = b.create_grid();
    }

    #[test]
    #[should_panic]
    fn test_import_wrong_degree() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 2));
        b.import_from_gmsh("meshes/sphere_triangle.msh");
        let _g = b.create_grid();
    }

    #[test]
    fn test_import_triangle_bin() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));
        b.import_from_gmsh("meshes/sphere_triangle_bin.msh");
        let _g = b.create_grid();
    }

    #[test]
    fn test_import_quadrilateral_bin() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));
        b.import_from_gmsh("meshes/cube_quadrilateral_bin.msh");
        let _g = b.create_grid();
    }

    #[test]
    fn test_import_tetrahedron_bin() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Tetrahedron, 1));
        b.import_from_gmsh("meshes/cube_tetrahedron_bin.msh");
        let _g = b.create_grid();
    }

    #[test]
    #[should_panic]
    fn test_import_wrong_cell_bin() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));
        b.import_from_gmsh("meshes/cube_tetrahedron_bin.msh");
        let _g = b.create_grid();
    }

    #[test]
    #[should_panic]
    fn test_import_wrong_degree_bin() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 2));
        b.import_from_gmsh("meshes/sphere_triangle_bin.msh");
        let _g = b.create_grid();
    }
}
