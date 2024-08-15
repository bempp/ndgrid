//! Parallel grid builder

use crate::{traits::{Builder, GeometryBuilder, TopologyBuilder}, types::Ownership};
use mpi::{
    point_to_point::{Destination, Source},
    request::WaitGuard,
    traits::{Buffer, Communicator, Equivalence},
};
use std::collections::HashMap;
use coupe::{KMeans, Partition, Point3D};
use itertools::izip;

pub trait ParallelBuilderFunctions: Builder + GeometryBuilder + TopologyBuilder {
    //! Parallel builder functions
    // TODO: replace f64 with T

    /// Partition the cells
    fn partition_cells(&self, nprocesses: usize) -> Vec<usize>;

    /// Assign vertex owners
    fn assign_vertex_owners(&self, cell_owners: &[usize]) -> Vec<usize>;

    /// Compute the vertices and cells that each process needs access to
    fn get_vertices_points_and_cells(
        &self,
        cell_owners: &[usize],
    ) -> (Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<Vec<usize>>);

    fn distribute_cells<C: Communicator>(
        &self,
        comm: &C,
        cells_per_proc: &[Vec<usize>],
        cell_ownwes: &[usize],
    ) -> (
        Vec<usize>,
        Vec<usize>,
        Vec<Self::EntityDescriptor>,
        Vec<usize>,
        Vec<usize>,
    );
    fn receive_cells<C: Communicator>(
        &self,
        comm: &C,
        root_rank: i32,
    ) -> (
        Vec<usize>,
        Vec<usize>,
        Vec<Self::EntityDescriptor>,
        Vec<usize>,
        Vec<usize>,
    );

    fn distribute_points<C: Communicator>(
        &self,
        comm: &C,
        points_per_proc: &[Vec<usize>],
    ) -> (Vec<usize>, Vec<Self::T>);
    fn receive_points<C: Communicator>(
        &self,
        comm: &C,
        root_rank: i32,
    ) -> (Vec<usize>, Vec<Self::T>);

    fn distribute_vertices<C: Communicator>(
        &self,
        comm: &C,
        vertices_per_proc: &[Vec<usize>],
        vertex_owners: &[usize],
    ) -> (Vec<usize>, Vec<usize>);
    fn receive_vertices<C: Communicator>(
        &self,
        comm: &C,
        root_rank: i32,
    ) -> (Vec<usize>, Vec<usize>);

    /// Reorder cells so that owned cells are first
    fn reorder_cells(
        &self,
        rank: usize,
        cell_indices: &[usize],
        cell_points: &[usize],
        cell_types: &[Self::EntityDescriptor],
        cell_degrees: &[usize],
        cell_owners: &[usize],
    ) -> (
        Vec<usize>,
        Vec<usize>,
        Vec<Self::EntityDescriptor>,
        Vec<usize>,
        Vec<usize>,
    );

    /// Reorder vertices so that owned vertices are first
    fn reorder_vertices(
        &self,
        rank: usize,
        vertex_indices: &[usize],
        vertex_owners: &[usize],
    ) -> (Vec<usize>, Vec<usize>);

    fn communicate_owners<C: Communicator>(&self, comm: &C, owners: &[usize], indices: &[usize]) -> Vec<Ownership>;
}

pub trait ParallelBuilder:
    Builder + ParallelBuilderFunctions + GeometryBuilder + TopologyBuilder
{
    //! MPI parallel grid builder
    // TODO: move to traits/builder.rs or traits/parallel.rs

    // TODO: remove this test
    fn test(&self) {
        use mpi::{environment::Universe, topology::Communicator};

        let universe: Universe = mpi::initialize().unwrap();
        let comm = universe.world();
        let rank = comm.rank();

        let (geometry, topology, vertex_owners, cell_owners) = if rank == 0 {
            let cell_owners = self.partition_cells(comm.size() as usize);
            let vertex_owners = self.assign_vertex_owners(&cell_owners);
            let (vertices_per_proc, points_per_proc, cells_per_proc) =
                self.get_vertices_points_and_cells(&cell_owners);

            let (point_indices, coordinates) = self.distribute_points(&comm, &points_per_proc);
            let (vertex_indices, vertex_owners) =
                self.distribute_vertices(&comm, &vertices_per_proc, &vertex_owners);
            let (cell_indices, cell_points, cell_types, cell_degrees, cell_owners) =
                self.distribute_cells(&comm, &cells_per_proc, &cell_owners);

            let geometry = self.create_geometry(
                &point_indices,
                &coordinates,
                &cell_points,
                &cell_types,
                &cell_degrees,
            );
            let topology =
                self.create_topology(&vertex_indices, &cell_indices, &cell_points, &cell_types);

            let vertex_owners = self.communicate_owners(&comm, &vertex_owners, &vertex_indices);
            let cell_owners = self.communicate_owners(&comm, &cell_owners, &cell_indices);

            (geometry, topology, vertex_owners, cell_owners)
        } else {
            let (point_indices, coordinates) = self.receive_points(&comm, 0);
            let (vertex_indices, vertex_owners) = self.receive_vertices(&comm, 0);
            let (cell_indices, cell_points, cell_types, cell_degrees, cell_owners) =
                self.receive_cells(&comm, 0);

            let geometry = self.create_geometry(
                &point_indices,
                &coordinates,
                &cell_points,
                &cell_types,
                &cell_degrees,
            );
            let topology =
                self.create_topology(&vertex_indices, &cell_indices, &cell_points, &cell_types);

            let vertex_owners = self.communicate_owners(&comm, &vertex_owners, &vertex_indices);
            let cell_owners = self.communicate_owners(&comm, &cell_owners, &cell_indices);

            (geometry, topology, vertex_owners, cell_owners)
        };

        println!("[{rank}] {vertex_owners:?}");
        println!("[{rank}] {cell_owners:?}");

        // Assign global DOFs to local DOFs
        // Collect vertex and cell global DOFs for ghosts

        // Edges (and other intermediates)
    }
}

impl<B: Builder + GeometryBuilder + TopologyBuilder> ParallelBuilderFunctions for B
where
    Vec<B::T>: Buffer,
    B::T: Equivalence,
    Vec<B::EntityDescriptor>: Buffer,
    B::EntityDescriptor: Equivalence,
{
    fn communicate_owners<C: Communicator>(&self, comm: &C, owners: &[usize], indices: &[usize]) -> Vec<Ownership> {
        let rank = comm.rank() as usize;
        let mut ownership = vec![Ownership::Undefined; owners.len()];
        let mut to_send = vec![vec![]; comm.size() as usize];
        let mut to_receive = vec![vec![]; comm.size() as usize];
        let mut ids_to_indices = HashMap::new();
        for (j, (o, i)) in izip!(owners, indices).enumerate() {
            if *o == rank {
                ids_to_indices.insert(i, j);
                ownership[j] = Ownership::Owned;
            } else {
                to_send[*o].push(*i);
                to_receive[*o].push(j);
            }
        }
        mpi::request::scope(|scope| {
            for p in 0..comm.size() {
                if p != comm.rank() {
                    let process = comm.process_at_rank(p);
                    let _ = WaitGuard::from(process.immediate_send(scope, &to_send[p as usize]));
                }
            }
        });
        let mut received = vec![];
        for p in 0..comm.size() {
            if p != comm.rank() {
                let process = comm.process_at_rank(p);
                let (r, _status) = process.receive_vec::<usize>();
                received.push(r);
            } else {
                received.push(vec![]);
            }
        }
        
        let mut to_send_back = received.iter().map(|r| r.iter().map(
            |i| ids_to_indices[i]
        ).collect::<Vec<_>>()).collect::<Vec<_>>();

        mpi::request::scope(|scope| {
            for p in 0..comm.size() {
                if p != comm.rank() {
                    let process = comm.process_at_rank(p);
                    let _ = WaitGuard::from(process.immediate_send(scope, &to_send_back[p as usize]));
                }
            }
        });

        for p in 0..comm.size() {
            if p != comm.rank() {
                let process = comm.process_at_rank(p);
                let (r, _status) = process.receive_vec::<usize>();
                for (i, j) in izip!(&to_receive[p as usize], &r) {
                    ownership[*i] = Ownership::Ghost(p as usize, *j);
                }
            }
        }

        ownership
    }

    fn reorder_vertices(
        &self,
        rank: usize,
        vertex_indices: &[usize],
        vertex_owners: &[usize],
    ) -> (Vec<usize>, Vec<usize>) {
        let mut indices = vec![];
        let mut owners = vec![];
        let mut indices2 = vec![];
        let mut owners2 = vec![];
        for (i, o) in izip!(vertex_indices, vertex_owners) {
            if *o == rank {
                indices.push(*i);
                owners.push(*o);
            } else {
                indices2.push(*i);
                owners2.push(*o);
            }
        }
        indices.extend_from_slice(&indices2);
        owners.extend_from_slice(&owners2);
        (indices, owners)
    }
    fn reorder_cells(
        &self,
        rank: usize,
        cell_indices: &[usize],
        cell_points: &[usize],
        cell_types: &[Self::EntityDescriptor],
        cell_degrees: &[usize],
        cell_owners: &[usize],
    ) -> (
        Vec<usize>,
        Vec<usize>,
        Vec<Self::EntityDescriptor>,
        Vec<usize>,
        Vec<usize>,
    ) {
        let mut indices = vec![];
        let mut points = vec![];
        let mut types = vec![];
        let mut degrees = vec![];
        let mut owners = vec![];
        let mut indices2 = vec![];
        let mut points2 = vec![];
        let mut types2 = vec![];
        let mut degrees2 = vec![];
        let mut owners2 = vec![];
        let mut start = 0;
        for (i, t, d, o) in izip!(
            cell_indices,
            cell_types,
            cell_degrees,
            cell_owners
        ) {
            let npts = self.npts(*t, *d);
            if *o == rank {
                indices.push(*i);
                points.extend_from_slice(&cell_points[start..start+npts]);
                types.push(*t);
                degrees.push(*d);
                owners.push(*o);
            } else {
                indices2.push(*i);
                points2.extend_from_slice(&cell_points[start..start+npts]);
                types2.push(*t);
                degrees2.push(*d);
                owners2.push(*o);
            }
            start += npts;
        }
        indices.extend_from_slice(&indices2);
        points.extend_from_slice(&points2);
        types.extend_from_slice(&types2);
        degrees.extend_from_slice(&degrees2);
        owners.extend_from_slice(&owners2);
        (indices, points, types, degrees, owners)
    }
    fn distribute_vertices<C: Communicator>(
        &self,
        comm: &C,
        vertices_per_proc: &[Vec<usize>],
        vertex_owners: &[usize],
    ) -> (Vec<usize>, Vec<usize>) {
        let rank = comm.rank() as usize;

        let vertex_owners_per_proc = vertices_per_proc
            .iter()
            .map(|vs| vs.iter().map(|v| vertex_owners[*v]).collect::<Vec<_>>())
            .collect::<Vec<_>>();

        mpi::request::scope(|scope| {
            for i in 0..comm.size() {
                if i != comm.rank() {
                    let process = comm.process_at_rank(i);
                    let _ = WaitGuard::from(
                        process.immediate_send(scope, &vertices_per_proc[i as usize]),
                    );
                    let _ = WaitGuard::from(
                        process.immediate_send(scope, &vertex_owners_per_proc[i as usize]),
                    );
                }
            }
        });

        self.reorder_vertices(rank,
            &vertices_per_proc[rank],
            &vertex_owners_per_proc[rank],
        )
    }
    fn receive_vertices<C: Communicator>(
        &self,
        comm: &C,
        root_rank: i32,
    ) -> (Vec<usize>, Vec<usize>) {
        let root_process = comm.process_at_rank(root_rank);
        let (vertices, _status) = root_process.receive_vec::<usize>();
        let (vertex_owners, _status) = root_process.receive_vec::<usize>();

        self.reorder_vertices(comm.rank() as usize, &vertices, &vertex_owners)
    }
    fn distribute_points<C: Communicator>(
        &self,
        comm: &C,
        points_per_proc: &[Vec<usize>],
    ) -> (Vec<usize>, Vec<Self::T>) {
        let rank = comm.rank() as usize;
        let mut coords_rank = vec![];
        let mut coords = vec![];
        for (p, points) in points_per_proc.iter().enumerate() {
            let mut coords_i = vec![];
            for i in points {
                coords_i.extend_from_slice(self.point(*i));
            }
            if p == rank {
                coords_rank = coords_i;
                coords.push(vec![]);
            } else {
                coords.push(coords_i);
            }
        }
        mpi::request::scope(|scope| {
            for i in 0..comm.size() {
                if i != comm.rank() {
                    let process = comm.process_at_rank(i);
                    let _ = WaitGuard::from(
                        process.immediate_send(scope, &points_per_proc[i as usize]),
                    );
                    let _ = WaitGuard::from(process.immediate_send(scope, &coords[i as usize]));
                }
            }
        });

        (points_per_proc[rank].clone(), coords_rank)
    }
    fn receive_points<C: Communicator>(
        &self,
        comm: &C,
        root_rank: i32,
    ) -> (Vec<usize>, Vec<Self::T>) {
        let root_process = comm.process_at_rank(root_rank);
        let (indices, _status) = root_process.receive_vec::<usize>();
        let (coords, _status) = root_process.receive_vec::<Self::T>();

        (indices, coords)
    }
    fn distribute_cells<C: Communicator>(
        &self,
        comm: &C,
        cells_per_proc: &[Vec<usize>],
        cell_owners: &[usize],
    ) -> (
        Vec<usize>,
        Vec<usize>,
        Vec<Self::EntityDescriptor>,
        Vec<usize>,
        Vec<usize>,
    ) {
        let rank = comm.rank() as usize;
        let mut cell_points_rank = vec![];
        let mut cell_types_rank = vec![];
        let mut cell_degrees_rank = vec![];
        let mut cell_owners_rank = vec![];

        let mut cell_points = vec![];
        let mut cell_types = vec![];
        let mut cell_degrees = vec![];
        let mut cell_owners_per_proc = vec![];
        for (p, cells) in cells_per_proc.iter().enumerate() {
            let mut cell_points_i = vec![];
            let mut cell_types_i = vec![];
            let mut cell_degrees_i = vec![];
            let mut cell_owners_i = vec![];
            for cell in cells {
                let pts = self.cell_points(*cell);
                cell_points_i.extend_from_slice(pts);
                cell_types_i.push(self.cell_type(*cell));
                cell_degrees_i.push(self.cell_degree(*cell));
                cell_owners_i.push(cell_owners[*cell]);
            }
            if p == rank {
                cell_points_rank = cell_points_i;
                cell_types_rank = cell_types_i;
                cell_degrees_rank = cell_degrees_i;
                cell_owners_rank = cell_owners_i;
                cell_points.push(vec![]);
                cell_types.push(vec![]);
                cell_degrees.push(vec![]);
                cell_owners_per_proc.push(vec![]);
            } else {
                cell_points.push(cell_points_i);
                cell_types.push(cell_types_i);
                cell_degrees.push(cell_degrees_i);
                cell_owners_per_proc.push(cell_owners_i);
            }
        }
        mpi::request::scope(|scope| {
            for i in 0..comm.size() {
                if i != comm.rank() {
                    let process = comm.process_at_rank(i);
                    let _ =
                        WaitGuard::from(process.immediate_send(scope, &cells_per_proc[i as usize]));
                    let _ =
                        WaitGuard::from(process.immediate_send(scope, &cell_points[i as usize]));
                    let _ = WaitGuard::from(process.immediate_send(scope, &cell_types[i as usize]));
                    let _ =
                        WaitGuard::from(process.immediate_send(scope, &cell_degrees[i as usize]));
                    let _ = WaitGuard::from(
                        process.immediate_send(scope, &cell_owners_per_proc[i as usize]),
                    );
                }
            }
        });

        self.reorder_cells(rank,
            &cells_per_proc[rank],
            &cell_points_rank,
            &cell_types_rank,
            &cell_degrees_rank,
            &cell_owners_rank,
        )
    }
    fn receive_cells<C: Communicator>(
        &self,
        comm: &C,
        root_rank: i32,
    ) -> (
        Vec<usize>,
        Vec<usize>,
        Vec<Self::EntityDescriptor>,
        Vec<usize>,
        Vec<usize>,
    ) {
        let root_process = comm.process_at_rank(root_rank);
        let (cell_indices, _status) = root_process.receive_vec::<usize>();
        let (cell_points, _status) = root_process.receive_vec::<usize>();
        let (cell_types, _status) = root_process.receive_vec::<Self::EntityDescriptor>();
        let (cell_degrees, _status) = root_process.receive_vec::<usize>();
        let (cell_owners, _status) = root_process.receive_vec::<usize>();

        self.reorder_cells(comm.rank() as usize,
            &cell_indices,
            &cell_points,
            &cell_types,
            &cell_degrees,
            &cell_owners,
        )
    }

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

    fn get_vertices_points_and_cells(
        &self,
        cell_owners: &[usize],
    ) -> (Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<Vec<usize>>) {
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
        let mut points = (0..nprocesses).map(|_| vec![]).collect::<Vec<_>>();
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

        for (vs, ps, cs) in izip!(vertices.iter_mut(), points.iter_mut(), &cells) {
            for c in cs {
                for v in self.cell_vertices(*c) {
                    if !vs.contains(v) {
                        vs.push(*v);
                    }
                }
                for v in self.cell_points(*c) {
                    if !ps.contains(v) {
                        ps.push(*v);
                    }
                }
            }
        }

        (vertices, points, cells)
    }
}
impl<B: Builder + GeometryBuilder + TopologyBuilder> ParallelBuilder for B
where
    Vec<B::T>: Buffer,
    B::T: Equivalence,
    Vec<B::EntityDescriptor>: Buffer,
    B::EntityDescriptor: Equivalence,
{
}
