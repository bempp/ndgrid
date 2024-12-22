//! Parallel grid builder

use super::ParallelGrid;
use crate::{
    traits::{
        Builder, Entity, GeometryBuilder, Grid, GridBuilder, ParallelBuilder, Topology,
        TopologyBuilder,
    },
    types::Ownership,
};
use coupe::{KMeans, Partition, Point3D};
use itertools::izip;
use mpi::{
    point_to_point::{Destination, Source},
    request::WaitGuard,
    traits::{Buffer, Communicator, Equivalence},
};
use std::collections::HashMap;

impl<B: Builder + GeometryBuilder + TopologyBuilder + GridBuilder> ParallelBuilder for B
where
    Vec<B::T>: Buffer,
    B::T: Equivalence,
    Vec<B::EntityDescriptor>: Buffer,
    B::EntityDescriptor: Equivalence,
    B::Grid: Sync,
{
    type ParallelGrid<'a, C: Communicator + 'a>
        = ParallelGrid<'a, C, B::Grid>
    where
        Self: 'a;
    fn create_parallel_grid<'a, C: Communicator>(
        &self,
        comm: &'a C,
    ) -> ParallelGrid<'a, C, B::Grid> {
        // Partition the cells via a KMeans algorithm. The midpoint of each cell is used and distributed
        // across processes via coupe. The returned array assigns each cell a process.
        let cell_owners = self.partition_cells(comm.size() as usize);
        // Each vertex is assigned the minimum process that has a cell that contains it.
        // Note that the array is of the size of the points in the grid and contains garbage information
        // for points that are not vertices.
        let vertex_owners = self.assign_vertex_owners(&cell_owners);

        // This distributes cells, vertices and points to processes.
        // Each process gets now only the cell that it owns via `cell_owners` but also its neighbours.
        // Then all the corresponding vertices and points are also added to the corresponding cell.
        // The layout of the Vec<Vec<usize>> structures is that the outer dimension is the process and
        // the inner dimension is the index, i.e. points_per_proc[0][j] is the jth point on process 0.
        let (vertices_per_proc, points_per_proc, cells_per_proc) =
            self.get_vertices_points_and_cells(&cell_owners);

        // Distribute the points and corresponding indices to all processes.
        let (point_indices, coordinates) = self.distribute_points(comm, &points_per_proc);

        // Returns the vertex indices on the current process and the corresponding owning processes.
        // Vertices are reordered so that owned vertices are first.
        let (vertex_indices, vertex_owners) =
            self.distribute_vertices(comm, &vertices_per_proc, &vertex_owners);

        // Distribute the cell information. Cells are reordered so that owned cells are first.
        let (cell_indices, cell_points, cell_types, cell_degrees, cell_owners) =
            self.distribute_cells(comm, &cells_per_proc, &cell_owners);

        // This is executed on all ranks and creates the local grid.
        self.create_parallel_grid_internal(
            comm,
            point_indices,
            coordinates,
            vertex_indices,
            vertex_owners,
            cell_indices,
            cell_points,
            cell_types,
            cell_degrees,
            cell_owners,
        )
    }
    fn receive_parallel_grid<'a, C: Communicator>(
        &self,
        comm: &'a C,
        root_rank: i32,
    ) -> ParallelGrid<'a, C, B::Grid> {
        let (point_indices, coordinates) = self.receive_points(comm, root_rank);
        let (vertex_indices, vertex_owners) = self.receive_vertices(comm, root_rank);
        let (cell_indices, cell_points, cell_types, cell_degrees, cell_owners) =
            self.receive_cells(comm, 0);

        self.create_parallel_grid_internal(
            comm,
            point_indices,
            coordinates,
            vertex_indices,
            vertex_owners,
            cell_indices,
            cell_points,
            cell_types,
            cell_degrees,
            cell_owners,
        )
    }
}

trait ParallelBuilderFunctions: Builder + GeometryBuilder + TopologyBuilder + GridBuilder {
    //! Parallel builder functions
    //!
    //! These functions are included in a trait so that they can be implemented for an arbitrary builder B below

    /// Intrernal function to create parallel grid
    #[allow(clippy::too_many_arguments)]
    fn create_parallel_grid_internal<'a, C: Communicator>(
        &self,
        comm: &'a C,
        point_indices: Vec<usize>,
        coordinates: Vec<Self::T>,
        vertex_indices: Vec<usize>,
        vertex_owners: Vec<usize>,
        cell_indices: Vec<usize>,
        cell_points: Vec<usize>,
        cell_types: Vec<Self::EntityDescriptor>,
        cell_degrees: Vec<usize>,
        cell_owners: Vec<usize>,
    ) -> ParallelGrid<'a, C, Self::Grid>
    where
        Self::Grid: Sync,
    {
        let geometry = self.create_geometry(
            &point_indices,
            &coordinates,
            &cell_points,
            &cell_types,
            &cell_degrees,
        );
        let topology = self.create_topology(
            vertex_indices.clone(),
            cell_indices.clone(),
            &cell_points,
            &cell_types,
        );

        let serial_grid = self.create_grid_from_topology_geometry(topology, geometry);

        // Create global indices and proper ownership information with ghost process and local index
        // on ghost process.
        // After execution `vertex_global_indices` and `cell_global_indices` have the global indices
        // of the `vertex_indices` and `cell_indices` arrays.
        // The global indices are not necessarily identical to the original grid indices.
        let (vertex_global_indices, vertex_owners) = self
            .assign_global_indices_and_communicate_owners(comm, &vertex_owners, &vertex_indices);
        let (cell_global_indices, cell_owners) =
            self.assign_global_indices_and_communicate_owners(comm, &cell_owners, &cell_indices);

        let mut owners = vec![vertex_owners];
        let mut global_indices = vec![vertex_global_indices];
        for dim in 1..self.tdim() {
            let (entity_global_indices, entity_owners) = self
                .assign_sub_entity_global_indices_and_owners(
                    comm,
                    &serial_grid,
                    &owners[0],
                    &global_indices[0],
                    &cell_owners,
                    dim,
                );
            owners.push(entity_owners);
            global_indices.push(entity_global_indices);
        }
        owners.push(cell_owners);
        global_indices.push(cell_global_indices);

        ParallelGrid::new(comm, serial_grid, owners, global_indices)
    }

    /// Assign global indices to entities that are neither vertices nor cells and communicator ownership information
    fn assign_sub_entity_global_indices_and_owners<C: Communicator>(
        &self,
        comm: &C,
        grid: &Self::Grid,
        vertex_ownership: &[Ownership],
        vertex_global_indices: &[usize],
        cell_ownership: &[Ownership],
        dim: usize,
    ) -> (Vec<usize>, Vec<Ownership>) {
        let rank = comm.rank() as usize;
        let mut entity_owners = vec![
            comm.size() as usize;
            grid.entity_types(dim)
                .iter()
                .map(|t| grid.entity_count(*t))
                .sum()
        ];
        let mut entity_indices = vec![0; entity_owners.len()];
        let mut vertices_to_local_index = HashMap::new();
        let mut to_ask_for = vec![vec![]; comm.size() as usize];
        let mut to_ask_for_sizes = vec![vec![]; comm.size() as usize];
        let mut to_ask_for_indices = vec![vec![]; comm.size() as usize];
        let mut global_index = if rank == 0 {
            0
        } else {
            let process = comm.process_at_rank(comm.rank() - 1);
            let (global_index, _status) = process.receive::<usize>();
            global_index
        };
        for (index, entity) in grid.entity_iter(dim).enumerate() {
            let local_v = entity.topology().sub_entity_iter(0).collect::<Vec<_>>();
            let mut global_v = local_v
                .iter()
                .map(|i| vertex_global_indices[*i])
                .collect::<Vec<_>>();
            global_v.sort();
            vertices_to_local_index.insert(global_v.clone(), index);
            let mut in_charge = comm.size() as usize;
            for v in &local_v {
                in_charge = usize::min(
                    in_charge,
                    match vertex_ownership[*v] {
                        Ownership::Owned => rank,
                        Ownership::Ghost(p, _) => p,
                        _ => {
                            panic!("Unsupported ownership: {:?}", vertex_ownership[*v]);
                        }
                    },
                );
            }
            if in_charge == rank {
                entity_owners[index] = rank;
                entity_indices[index] = global_index;
                global_index += 1;
            } else {
                to_ask_for[in_charge].extend_from_slice(&global_v);
                to_ask_for_sizes[in_charge].push(global_v.len());
                to_ask_for_indices[in_charge].push(index);
            }
        }
        for (cell, ownership) in izip!(grid.entity_iter(grid.topology_dim()), cell_ownership) {
            for entity in cell.topology().sub_entity_iter(dim) {
                let mut global_v = grid
                    .entity(dim, entity)
                    .unwrap()
                    .topology()
                    .sub_entity_iter(0)
                    .map(|i| vertex_global_indices[i])
                    .collect::<Vec<_>>();
                global_v.sort();
                if let Ownership::Ghost(p, _) = ownership {
                    if let Some(i) = vertices_to_local_index.get(&global_v) {
                        entity_owners[*i] = usize::min(entity_owners[*i], *p);
                    }
                }
            }
        }

        mpi::request::scope(|scope| {
            if comm.rank() + 1 < comm.size() {
                let process = comm.process_at_rank(comm.rank() + 1);
                let _ = WaitGuard::from(process.immediate_send(scope, &global_index));
            }
            for p in 0..comm.size() {
                if p != comm.rank() {
                    let process = comm.process_at_rank(p);
                    let _ = WaitGuard::from(process.immediate_send(scope, &to_ask_for[p as usize]));
                    let _ = WaitGuard::from(
                        process.immediate_send(scope, &to_ask_for_sizes[p as usize]),
                    );
                }
            }
        });
        let mut send_back_owners = vec![];
        let mut send_back_indices = vec![];
        for p in 0..comm.size() {
            if p != comm.rank() {
                let process = comm.process_at_rank(p);
                let (asked, _status) = process.receive_vec::<usize>();
                let (sizes, _status) = process.receive_vec::<usize>();

                let mut owners = vec![];
                let mut indices = vec![];
                let mut start = 0;
                for s in sizes {
                    owners.push(entity_owners[vertices_to_local_index[&asked[start..start + s]]]);
                    indices.push(entity_indices[vertices_to_local_index[&asked[start..start + s]]]);
                    start += s;
                }
                send_back_owners.push(owners);
                send_back_indices.push(indices);
            } else {
                send_back_owners.push(vec![]);
                send_back_indices.push(vec![]);
            }
        }
        mpi::request::scope(|scope| {
            for p in 0..comm.size() {
                if p != comm.rank() {
                    let process = comm.process_at_rank(p);
                    let _ = WaitGuard::from(
                        process.immediate_send(scope, &send_back_owners[p as usize]),
                    );
                    let _ = WaitGuard::from(
                        process.immediate_send(scope, &send_back_indices[p as usize]),
                    );
                }
            }
        });

        for p in 0..comm.size() {
            if p != comm.rank() {
                let process = comm.process_at_rank(p);
                let (owners, _status) = process.receive_vec::<usize>();
                let (indices, _status) = process.receive_vec::<usize>();
                for (local_i, o, i) in izip!(&to_ask_for_indices[p as usize], &owners, &indices) {
                    entity_owners[*local_i] = *o;
                    entity_indices[*local_i] = *i;
                }
            }
        }
        self.assign_global_indices_and_communicate_owners(comm, &entity_owners, &entity_indices)
    }
    /// Assign global indices to vertices or cells and communicator ownership information
    fn assign_global_indices_and_communicate_owners<C: Communicator>(
        &self,
        comm: &C,
        owners: &[usize],
        indices: &[usize],
    ) -> (Vec<usize>, Vec<Ownership>) {
        let rank = comm.rank() as usize;
        let mut ownership = vec![Ownership::Undefined; owners.len()];
        let mut global_indices = vec![0; owners.len()];
        let mut to_send = vec![vec![]; comm.size() as usize];
        let mut to_receive = vec![vec![]; comm.size() as usize];
        let mut ids_to_indices = HashMap::new();
        // Global indices start at 0 for the first process and then are contiguously increased across processes.
        // This is very inefficient since it effectively serializes the setup.
        let mut global_index = if rank == 0 {
            0
        } else {
            let process = comm.process_at_rank(comm.rank() - 1);
            let (global_index, _status) = process.receive::<usize>();
            global_index
        };
        for (j, (o, i)) in izip!(owners, indices).enumerate() {
            if *o == rank {
                ids_to_indices.insert(i, j);
                ownership[j] = Ownership::Owned;
                global_indices[j] = global_index;
                global_index += 1
            } else {
                // These are the indices that are sent to the other process.
                to_send[*o].push(*i);
                // These are the corresponding local indices of the indices
                // that are sent out.
                to_receive[*o].push(j);
            }
        }
        if comm.rank() + 1 < comm.size() {
            mpi::request::scope(|scope| {
                let process = comm.process_at_rank(comm.rank() + 1);
                let _ = WaitGuard::from(process.immediate_send(scope, &global_index));
            });
        }
        mpi::request::scope(|scope| {
            for p in 0..comm.size() {
                if p != comm.rank() {
                    let process = comm.process_at_rank(p);
                    // Send my indices to the actual owning process.
                    let _ = WaitGuard::from(process.immediate_send(scope, &to_send[p as usize]));
                }
            }
        });
        let mut received = vec![];
        for p in 0..comm.size() {
            if p != comm.rank() {
                let process = comm.process_at_rank(p);
                // Receive the indices that are owned by me but are sent to me from other processes.
                let (r, _status) = process.receive_vec::<usize>();
                received.push(r);
            } else {
                received.push(vec![]);
            }
        }

        // I am the owning process of the indices that are sent to me. I need send back to the other
        // process which local indices the indices are associated with that they sent to me.
        let to_send_back = received
            .iter()
            .map(|r| r.iter().map(|i| ids_to_indices[i]).collect::<Vec<_>>())
            .collect::<Vec<_>>();

        // I am also sending back the global indices of my owned indices.
        let global_to_send_back = to_send_back
            .iter()
            .map(|r| r.iter().map(|i| global_indices[*i]).collect::<Vec<_>>())
            .collect::<Vec<_>>();

        // Now do the actual send operations.
        mpi::request::scope(|scope| {
            for p in 0..comm.size() {
                if p != comm.rank() {
                    let process = comm.process_at_rank(p);
                    let _ =
                        WaitGuard::from(process.immediate_send(scope, &to_send_back[p as usize]));
                    let _ = WaitGuard::from(
                        process.immediate_send(scope, &global_to_send_back[p as usize]),
                    );
                }
            }
        });

        // And do the actual receive operations.
        for p in 0..comm.size() {
            if p != comm.rank() {
                let process = comm.process_at_rank(p);
                let (r, _status) = process.receive_vec::<usize>();
                let (g, _status) = process.receive_vec::<usize>();
                for (i, j, k) in izip!(&to_receive[p as usize], &r, &g) {
                    ownership[*i] = Ownership::Ghost(p as usize, *j);
                    global_indices[*i] = *k;
                }
            }
        }

        // Return the global indices and the associated ownership.
        (global_indices, ownership)
    }

    /// Reorder vertices so that owned vertices are first
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
    /// Reorder cells so that owned cells are first
    #[allow(clippy::type_complexity)]
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
        for (i, t, d, o) in izip!(cell_indices, cell_types, cell_degrees, cell_owners) {
            let npts = self.npts(*t, *d);
            if *o == rank {
                indices.push(*i);
                points.extend_from_slice(&cell_points[start..start + npts]);
                types.push(*t);
                degrees.push(*d);
                owners.push(*o);
            } else {
                indices2.push(*i);
                points2.extend_from_slice(&cell_points[start..start + npts]);
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
    /// Distribute vertices to all processes
    fn distribute_vertices<C: Communicator>(
        &self,
        comm: &C,
        vertices_per_proc: &[Vec<usize>],
        vertex_owners: &[usize],
    ) -> (Vec<usize>, Vec<usize>) {
        let rank = comm.rank() as usize;

        // Assign for each vertex on a proc the corresponding processess that owns the vertex.
        let vertex_owners_per_proc = vertices_per_proc
            .iter()
            .map(|vs| vs.iter().map(|v| vertex_owners[*v]).collect::<Vec<_>>())
            .collect::<Vec<_>>();

        // Send everything around.
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

        self.reorder_vertices(
            rank,
            &vertices_per_proc[rank],
            &vertex_owners_per_proc[rank],
        )
    }
    /// Receive vertices from root process
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
    /// Distribute points to all processes
    fn distribute_points<C: Communicator>(
        &self,
        comm: &C,
        points_per_proc: &[Vec<usize>],
    ) -> (Vec<usize>, Vec<Self::T>)
    where
        Vec<Self::T>: Buffer,
    {
        let rank = comm.rank() as usize;
        let mut coords_rank = vec![];
        let mut coords = vec![];
        // Get the coordinates of all points for all processes.
        // coords is a Vec<Vec<T>> where outer vec is the process
        // and the inner vec is the coordinates of the points.
        // For points on `rank` it just pushes an empty vec and stores
        // the coordinates in `coords_rank`.
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
        // Send the points to the other processes, one after another.
        // This could be much better handled via a broadcast operation.
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

        // Returns the indices of the local coordinates and the corresponding points.
        (points_per_proc[rank].clone(), coords_rank)
    }
    /// Receive points from root process
    fn receive_points<C: Communicator>(
        &self,
        comm: &C,
        root_rank: i32,
    ) -> (Vec<usize>, Vec<Self::T>)
    where
        Self::T: Equivalence,
    {
        let root_process = comm.process_at_rank(root_rank);
        let (indices, _status) = root_process.receive_vec::<usize>();
        let (coords, _status) = root_process.receive_vec::<Self::T>();

        (indices, coords)
    }
    /// Distribute cells to all processes
    #[allow(clippy::type_complexity)]
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
    )
    where
        Vec<Self::EntityDescriptor>: Buffer,
    {
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

        self.reorder_cells(
            rank,
            &cells_per_proc[rank],
            &cell_points_rank,
            &cell_types_rank,
            &cell_degrees_rank,
            &cell_owners_rank,
        )
    }
    /// Receive cells from root process
    #[allow(clippy::type_complexity)]
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
    )
    where
        Self::EntityDescriptor: Equivalence,
    {
        let root_process = comm.process_at_rank(root_rank);
        let (cell_indices, _status) = root_process.receive_vec::<usize>();
        let (cell_points, _status) = root_process.receive_vec::<usize>();
        let (cell_types, _status) = root_process.receive_vec::<Self::EntityDescriptor>();
        let (cell_degrees, _status) = root_process.receive_vec::<usize>();
        let (cell_owners, _status) = root_process.receive_vec::<usize>();

        self.reorder_cells(
            comm.rank() as usize,
            &cell_indices,
            &cell_points,
            &cell_types,
            &cell_degrees,
            &cell_owners,
        )
    }

    /// Partition the cells
    fn partition_cells(&self, nprocesses: usize) -> Vec<usize> {
        // Create an initial partitioning that roughly assigns the same
        // number of cells to each process.
        let mut partition = vec![];
        for i in 0..nprocesses {
            while partition.len() < self.cell_count() * (i + 1) / nprocesses {
                partition.push(i);
            }
        }

        // Compute the midpoints of each cell. If the geometric dimension is smaller than 3
        // then the remaining dimensions are just set to 0.
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

        // Assign equal weights to each cell
        let weights = vec![1.0; self.point_count()];

        // Run the Coupe KMeans algorithm to create the partitioning. Maybe replace
        // by Metis at some point.
        KMeans {
            delta_threshold: 0.0,
            ..Default::default()
        }
        .partition(&mut partition, (&midpoints, &weights))
        .unwrap();

        // Return the new partitioning.
        partition
    }

    /// Assign vertex owners
    fn assign_vertex_owners(&self, cell_owners: &[usize]) -> Vec<usize> {
        // Initialise empty array with number of processes as value and length same as number of points.
        let mut vertex_owners = vec![cell_owners.iter().max().unwrap() + 1; self.point_count()];

        // Each vertex is owned by the minimum process that owns a cell that contains the vertex.
        for (i, owner) in cell_owners.iter().enumerate() {
            for v in self.cell_vertices(i) {
                vertex_owners[*v] = usize::min(vertex_owners[*v], *owner);
            }
        }

        // Return the vertex owners. Note that the array contains garbage information
        // for the points of the grids that are not vertices.
        vertex_owners
    }

    /// Compute the vertices and cells that each process needs access to
    #[allow(clippy::type_complexity)]
    fn get_vertices_points_and_cells(
        &self,
        cell_owners: &[usize],
    ) -> (Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<Vec<usize>>) {
        // Each point in the grid is associated with a list of cells that contain it.
        let mut vertex_to_cells = (0..self.point_count()).map(|_| vec![]).collect::<Vec<_>>();

        // For each cell, add the cell to the list of cells that contain each of its vertices.
        for i in 0..self.cell_count() {
            for v in self.cell_vertices(i) {
                if !vertex_to_cells[*v].contains(&i) {
                    vertex_to_cells[*v].push(i);
                }
            }
        }

        // Initialise empty arrays for vertices, points and cells for each process.
        let nprocesses = *cell_owners.iter().max().unwrap() + 1;
        let mut vertices = (0..nprocesses).map(|_| vec![]).collect::<Vec<_>>();
        let mut points = (0..nprocesses).map(|_| vec![]).collect::<Vec<_>>();
        let mut cells = (0..nprocesses).map(|_| vec![]).collect::<Vec<_>>();

        // This assigns cell to processes. It iterates through a cell and associates
        // not only the cell itself to the process that owns it but also adds all the
        // neighboring cells to the same process. This means that processes own not
        // only the original cells from `cell_owners` but also all the neighbouring cells.
        // Also this means that a cell can live on multiple processes.
        // TODO! Use sets instead of vectors for the cells.
        for (i, owner) in cell_owners.iter().enumerate() {
            for v in self.cell_vertices(i) {
                for cell in &vertex_to_cells[*v] {
                    if !cells[*owner].contains(cell) {
                        cells[*owner].push(*cell);
                    }
                }
            }
        }

        // Now go through and add the points and vertices to the same processes that also the cells are
        // assigned to.
        // TODO! This is again superexpensive. Replace the vectors with sets.
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

impl<B: Builder + GeometryBuilder + TopologyBuilder + GridBuilder> ParallelBuilderFunctions for B
where
    Vec<B::T>: Buffer,
    B::T: Equivalence,
    Vec<B::EntityDescriptor>: Buffer,
    B::EntityDescriptor: Equivalence,
{
}
