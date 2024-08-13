//! Parallel grid builder

use crate::traits::Builder;
use ndelement::types::ReferenceCellType;

use coupe::{KMeans, Partition, Point3D};
use itertools::izip;

pub trait ParallelBuilderFunctions: Builder {
    //! Parallel builder functions

    /// Partition the cells
    fn partition_cells(&self, nprocesses: usize) -> Vec<usize>;

    /// Assign vertex owners
    fn assign_vertex_owners(&self, cell_owners: &[usize]) -> Vec<usize>;

    /// Compute the vertices and cells that each process needs access to
    fn get_vertices_and_cells(&self, cell_owners: &[usize]) -> (Vec<Vec<usize>>, Vec<Vec<usize>>);
}

pub trait ParallelBuilder: Builder + ParallelBuilderFunctions {
    //! MPI parallel grid builder
    // TODO: move to traits/builder.rs or traits/parallel.rs

    // TODO: remove this test
    fn test(&self) {
        let cell_owners = self.partition_cells(4);
        println!("Cell owners:   {cell_owners:?}");

        let vertex_owners = self.assign_vertex_owners(&cell_owners);
        println!("Vertex owners: {vertex_owners:?}");

        let (vertices_per_proc, cells_per_proc) = self.get_vertices_and_cells(&cell_owners);
        for i in 0..4 {
            println!("[{i}] Vertices:  {:?}", vertices_per_proc[i]);
            println!("    Cells:     {:?}", cells_per_proc[i]);
        }

        // Distribute point coordinates for each cell
        // Distribute list of points for each cell
        // Create local geometry

        // Distribute vertex indices
        // Distribute vertex owners
        // Distribute cell indices
        // Distribute cell vertex indices
        // Distribute cell owners

        // Create local topology
        // Assign global DOFs to local DOFs
        // Collect vertex and cell global DOFs for ghosts

        // Edges (and other intermediates
    }
}

impl<B: Builder<EntityDescriptor = ReferenceCellType>> ParallelBuilderFunctions for B {
    fn partition_cells(&self, nprocesses: usize) -> Vec<usize> {
        let mut partition = vec![];
        for i in 0..nprocesses {
            while partition.len() < self.cell_count() * (i + 1) / nprocesses {
                partition.push(i);
            }
        }

        let midpoints = (0..self.cell_count())
            .map(|i| {
                let v = self.cell_points(i);
                Point3D::new(
                    if self.gdim() > 0 {
                        num::cast::<Self::T, f64>(
                            v.iter().map(|j| self.point(*j)[0]).sum::<Self::T>(),
                        )
                        .unwrap()
                            / v.len() as f64
                    } else {
                        0.0
                    },
                    if self.gdim() > 1 {
                        num::cast::<Self::T, f64>(
                            v.iter().map(|j| self.point(*j)[1]).sum::<Self::T>(),
                        )
                        .unwrap()
                            / v.len() as f64
                    } else {
                        0.0
                    },
                    if self.gdim() > 2 {
                        num::cast::<Self::T, f64>(
                            v.iter().map(|j| self.point(*j)[2]).sum::<Self::T>(),
                        )
                        .unwrap()
                            / v.len() as f64
                    } else {
                        0.0
                    },
                )
            })
            .collect::<Vec<_>>();

        let weights = vec![1.0; self.point_count()];

        KMeans {
            part_count: nprocesses,
            delta_threshold: 0.0,
            ..Default::default()
        }
        .partition(&mut partition, (&midpoints, &weights))
        .unwrap();

        partition
    }

    fn assign_vertex_owners(&self, cell_owners: &[usize]) -> Vec<usize> {
        let mut vertex_owners = vec![cell_owners.iter().max().unwrap() + 1; self.point_count()];

        for (i, owner) in cell_owners.iter().enumerate() {
            for v in self.cell_vertices(i) {
                vertex_owners[*v] = usize::min(vertex_owners[*v], *owner);
            }
        }

        vertex_owners
    }

    fn get_vertices_and_cells(&self, cell_owners: &[usize]) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
        let mut vertex_to_cells = (0..self.point_count()).map(|_| vec![]).collect::<Vec<_>>();

        for i in 0..self.cell_count() {
            for v in self.cell_vertices(i) {
                if !vertex_to_cells[*v].contains(&i) {
                    vertex_to_cells[*v].push(i);
                }
            }
        }

        let nprocesses = *cell_owners.iter().max().unwrap() + 1;
        let mut vertices = (0..nprocesses).map(|_| vec![]).collect::<Vec<_>>();
        let mut cells = (0..nprocesses).map(|_| vec![]).collect::<Vec<_>>();

        for (i, owner) in cell_owners.iter().enumerate() {
            for v in self.cell_vertices(i) {
                for cell in &vertex_to_cells[*v] {
                    if !cells[*owner].contains(cell) {
                        cells[*owner].push(*cell);
                    }
                }
            }
        }

        for (vs, cs) in izip!(vertices.iter_mut(), &cells) {
            for c in cs {
                for v in self.cell_vertices(*c) {
                    if !vs.contains(v) {
                        vs.push(*v);
                    }
                }
            }
        }

        vertices[0].push(0);
        cells[0].push(0);

        (vertices, cells)
    }
}
impl<B: Builder<EntityDescriptor = ReferenceCellType>> ParallelBuilder for B {}
