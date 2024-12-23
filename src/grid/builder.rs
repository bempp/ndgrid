//! Parallel grid builder

use super::ParallelGrid;
use crate::{
    traits::{
        Builder, Entity, GeometryBuilder, Grid, GridBuilder, ParallelBuilder, Topology,
        TopologyBuilder,
    },
    types::Ownership,
};
use itertools::{izip, Chunk, Itertools};
use mpi::{
    collective::SystemOperation,
    point_to_point::{Destination, Source},
    request::WaitGuard,
    traits::{Buffer, Communicator, CommunicatorCollectives, Equivalence, Root},
};
use std::{
    cell,
    collections::{HashMap, HashSet},
};

// A simple struct to hold chunked data. The data is a long
// array and `idx_bounds` is a vector sunch that `chunks[idx_bounds[i]..idx_bounds[i+1]]`
// is one chunk. The last element of `idx_bounds` is the length of the data array.
struct ChunkedData<T: Equivalence + Copy> {
    data: Vec<T>,
    idx_bounds: Vec<usize>,
}

impl<T: Equivalence + Copy> ChunkedData<T> {
    // Return a vector with the number of elements per chunk.
    fn counts(&self) -> Vec<usize> {
        self.idx_bounds
            .iter()
            .tuple_windows()
            .map(|(a, b)| b - a)
            .collect()
    }
}

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
    fn create_parallel_grid_root<'a, C: Communicator>(
        &self,
        comm: &'a C,
    ) -> ParallelGrid<'a, C, B::Grid> {
        // Partition the cells via a KMeans algorithm. The midpoint of each cell is used and distributed
        // across processes via coupe. The returned array assigns each cell a process.
        let cell_owners = self.partition_cells(comm.size() as usize);
        // Each vertex is assigned the minimum process that has a cell that contains it.
        // Note that the array is of the size of the points in the grid and contains garbage information
        // for points that are not vertices.
        let vertex_owners = self.assign_vertex_owners(comm.size() as usize, &cell_owners);

        // This distributes cells, vertices and points to processes.
        // Each process gets now only the cell that it owns via `cell_owners` but also its neighbours.
        // Then all the corresponding vertices and points are also added to the corresponding cell.
        // Each return value is a tuple consisting of a counts array that specifies how many indices are
        // assigned to each process and the actual indices.
        let (vertices_per_proc, points_per_proc, cells_per_proc) =
            self.get_vertices_points_and_cells(comm.size() as usize, &cell_owners);

        // Compute the chunks array for the coordinates associated with the points.
        // The idx array for the coords is simply `self.gdim()` times the idx array for the points.
        let coords_per_proc = ChunkedData {
            data: {
                let mut tmp =
                    Vec::with_capacity(self.gdim() * points_per_proc.idx_bounds.last().unwrap());
                points_per_proc
                    .data
                    .iter()
                    .for_each(|point| tmp.extend(self.point(*point).iter()));
                tmp
            },
            idx_bounds: points_per_proc
                .idx_bounds
                .iter()
                .map(|x| self.gdim() * x)
                .collect(),
        };

        // This compputes for each process the vertex owners.
        let vertex_owners_per_proc = ChunkedData::<usize> {
            data: {
                let mut tmp = Vec::with_capacity(*vertices_per_proc.idx_bounds.last().unwrap());
                vertices_per_proc
                    .data
                    .iter()
                    .for_each(|v| tmp.push(vertex_owners[*v]));
                tmp
            },
            idx_bounds: vertices_per_proc.idx_bounds.clone(),
        };

        // We now compute the cell information for each process.

        // First we need the cell type.
        let cell_types_per_proc = ChunkedData {
            data: {
                let mut tmp = Vec::with_capacity(*cells_per_proc.idx_bounds.last().unwrap());
                cells_per_proc
                    .data
                    .iter()
                    .for_each(|c| tmp.push(self.cell_type(*c)));
                tmp
            },
            idx_bounds: cells_per_proc.idx_bounds.clone(),
        };
        // Now we need the cell degrees.
        let cell_degrees_per_proc = ChunkedData {
            data: {
                let mut tmp = Vec::with_capacity(*cells_per_proc.idx_bounds.last().unwrap());
                cells_per_proc
                    .data
                    .iter()
                    .for_each(|c| tmp.push(self.cell_degree(*c)));
                tmp
            },
            idx_bounds: cells_per_proc.idx_bounds.clone(),
        };
        // Now need the cell owners.
        let cell_owners_per_proc = ChunkedData {
            data: {
                let mut tmp = Vec::with_capacity(*cells_per_proc.idx_bounds.last().unwrap());
                cells_per_proc
                    .data
                    .iter()
                    .for_each(|c| tmp.push(cell_owners[*c]));
                tmp
            },
            idx_bounds: cells_per_proc.idx_bounds.clone(),
        };
        // Finally the cell points. These are more messy since each cell might have different numbers of points. Hence,
        // we need to compute a new counts array.

        let cell_points_per_proc = {
            // First compute the total number of points needed so that we can pre-allocated the array.
            let total_points = cells_per_proc
                .data
                .iter()
                .map(|c| self.npts(self.cell_type(*c), self.cell_degree(*c)))
                .sum();
            let mut data = Vec::with_capacity(total_points);
            let mut counts = vec![0; comm.size() as usize];
            for (proc, (&chunk_ind_start, &chunk_ind_end)) in
                cells_per_proc.idx_bounds.iter().tuple_windows().enumerate()
            {
                for c in cells_per_proc.data[chunk_ind_start..chunk_ind_end].iter() {
                    let npts = self.npts(self.cell_type(*c), self.cell_degree(*c));
                    counts[proc] += npts;
                    data.extend(self.cell_points(*c).iter());
                }
            }
            // Now turn the counts into an `idx_bounds` array.
            let mut idx_bounds = counts
                .iter()
                .scan(0, |state, &x| {
                    let old = *state;
                    *state += x;
                    Some(old)
                })
                .collect_vec();
            idx_bounds.push(idx_bounds.last().unwrap() + counts.last().unwrap());

            // Finally return the chunks
            ChunkedData { data, idx_bounds }
        };

        // Now we send everything across the processes. Every _per_proc variable is scattered across
        // the processes. After the scatter operation we have on each process the following data.
        // - `point_indices` contains the indices of the points that are owned by the current process.
        // - `coordinates` contains the coordinates of the points that are owned by the current process.
        // - `vertex_indices` contains the indices of the vertices that are owned by the current process.
        // - `vertex_owners` contains the owners of the vertices that are owned by the current process.
        // - `cell_indices` contains the indices of the cells that are owned by the current process.
        // - `cell_points` contains the indices of the points of the cells that are owned by the current process.
        // - `cell_types` contains the types of the cells that are owned by the current process.
        // - `cell_degrees` contains the degrees of the cells that are owned by the current process.
        // - `cell_owners` contains the owners of the cells that are owned by the current process.

        let point_indices = scatterv_root(comm, &points_per_proc);
        let coordinates = scatterv_root(comm, &coords_per_proc);
        let vertex_indices = scatterv_root(comm, &vertices_per_proc);
        let vertex_owners = scatterv_root(comm, &vertex_owners_per_proc);
        let cell_indices = scatterv_root(comm, &cells_per_proc);
        let cell_points = scatterv_root(comm, &cell_points_per_proc);
        let cell_types = scatterv_root(comm, &cell_types_per_proc);
        let cell_degrees = scatterv_root(comm, &cell_degrees_per_proc);
        let cell_owners = scatterv_root(comm, &cell_owners_per_proc);

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
    fn create_parallel_grid<'a, C: Communicator>(
        &self,
        comm: &'a C,
        root_rank: i32,
    ) -> ParallelGrid<'a, C, B::Grid> {
        // First we receive all the data.
        let point_indices = scatterv(comm, 0);
        let coordinates = scatterv(comm, 0);
        let vertex_indices = scatterv(comm, 0);
        let vertex_owners = scatterv(comm, 0);
        let cell_indices = scatterv(comm, 0);
        let cell_points = scatterv(comm, 0);
        let cell_types = scatterv(comm, 0);
        let cell_degrees = scatterv(comm, 0);
        let cell_owners = scatterv(comm, 0);

        // Now we reassemble the grid.

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
        let rank = comm.rank() as usize;
        // First we want to reorder everything so that owned data comes first in the arrays.

        // Reorder the vertices.
        let (vertex_indices, vertex_owners) = {
            let mut new_indices = Vec::<usize>::with_capacity(vertex_indices.len());
            let mut new_owners = Vec::<usize>::with_capacity(vertex_owners.len());

            for (&v, &o) in izip!(&vertex_indices, &vertex_owners).filter(|(_, &o)| o == rank) {
                new_indices.push(v);
                new_owners.push(o);
            }

            for (&v, &o) in izip!(&vertex_indices, &vertex_owners).filter(|(_, &o)| o != rank) {
                new_indices.push(v);
                new_owners.push(o);
            }

            (new_indices, new_owners)
        };

        // Reorder the cell information. However, things are a bit messy with the cell points as each cell
        // may have a different number of points attached. Hence, we need to compute a new counts array for the cell points.

        let cell_point_owners = {
            let mut new_owners = Vec::<usize>::with_capacity(cell_points.len());
            // We now iterate through the cells and their owners and assign owners accordingly for the points.
            for (c, o) in izip!(&cell_indices, &cell_owners) {
                let npts = self.npts(cell_types[*c], cell_degrees[*c]);
                for _ in 0..npts {
                    new_owners.push(*o);
                }
            }
            new_owners
        };

        // First we reorder the cells, degrees, types and owners.

        let (cell_indices, cell_types, cell_degrees, cell_owners) = {
            let mut new_indices = Vec::<usize>::with_capacity(cell_indices.len());
            let mut new_types = Vec::<Self::EntityDescriptor>::with_capacity(cell_types.len());
            let mut new_degrees = Vec::<usize>::with_capacity(cell_degrees.len());
            let mut new_owners = Vec::<usize>::with_capacity(cell_owners.len());

            for (cell_index, cell_type, cell_degree, cell_owner) in
                izip!(&cell_indices, &cell_types, &cell_degrees, &cell_owners,)
                    .filter(|(_, _, _, &o)| o == rank)
            {
                new_indices.push(*cell_index);
                new_types.push(*cell_type);
                new_degrees.push(*cell_degree);
                new_owners.push(*cell_owner);
            }

            for (cell_index, cell_type, cell_degree, cell_owner) in
                izip!(&cell_indices, &cell_types, &cell_degrees, &cell_owners,)
                    .filter(|(_, _, _, &o)| o != rank)
            {
                new_indices.push(*cell_index);
                new_types.push(*cell_type);
                new_degrees.push(*cell_degree);
                new_owners.push(*cell_owner);
            }

            (new_indices, new_types, new_degrees, new_owners)
        };

        // Finally, we reorder the cell points.
        let cell_points = {
            let mut new_points = Vec::<usize>::with_capacity(cell_points.len());
            for &p in izip!(&cell_points, &cell_point_owners)
                .filter(|(_, &o)| o == rank)
                .map(|(p, _)| p)
            {
                new_points.push(p);
            }

            for &p in izip!(&cell_points, &cell_point_owners)
                .filter(|(_, &o)| o != rank)
                .map(|(p, _)| p)
            {
                new_points.push(p);
            }
            new_points
        };

        // Everything is properly sorted now. Can generate the local grids.
        // Need to check those routines. Implementation is not efficient right now.
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
        let (vertex_global_indices, vertex_owners) =
            synchronize_indices(comm, &vertex_indices, &vertex_owners);
        let (cell_global_indices, cell_owners) =
            synchronize_indices(comm, &cell_indices, &cell_owners);

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
            // Get vertices of entity.
            let local_v = entity.topology().sub_entity_iter(0).collect::<Vec<_>>();
            // Assign the global dofs to the entity vertices.
            let mut global_v = local_v
                .iter()
                .map(|i| vertex_global_indices[*i])
                .collect::<Vec<_>>();
            global_v.sort();
            // Use the sorted global vertices as key that maps to the local index of the entity.
            vertices_to_local_index.insert(global_v.clone(), index);
            let mut in_charge = comm.size() as usize;
            // The ownership of the entity is assigned to minimum process that owns one of the
            // vertices of the entity.
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
            // If the entity is local simply assign the global index.
            if in_charge == rank {
                entity_owners[index] = rank;
                entity_indices[index] = global_index;
                global_index += 1;
            } else {
                // Otherwise, we need to ask the process that owns the entity for its global index.
                to_ask_for[in_charge].extend_from_slice(&global_v);
                to_ask_for_sizes[in_charge].push(global_v.len());
                to_ask_for_indices[in_charge].push(index);
            }
        }
        // Iterate through all the cells. This updates entity ownership to being owned by a cell
        // it is adjacent to if the cell has a lower process number than the current owner.
        // I don't think this is necessary. A vertex is already owned by the same process as the smallest
        // process cell it is adjacent to, and the smallest ownership transfers to the entity because of the above.
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

    /// Reorder data so that data on the current rank comes first.
    fn reorder<T: Equivalence + Copy>(
        &self,
        rank: usize,
        data: &[T],
        owners: &[usize],
    ) -> (Vec<T>, Vec<usize>) {
        assert_eq!(data.len(), owners.len());
        let mut new_data = Vec::<T>::with_capacity(data.len());
        let mut new_owners = Vec::<usize>::with_capacity(owners.len());

        // First push all local elements.
        for (d, o) in izip!(data, owners).filter(|(_, &o)| o == rank) {
            new_data.push(*d);
            new_owners.push(*o);
        }
        // Now push all the other elements.
        for (d, o) in izip!(data, owners).filter(|(_, &o)| o != rank) {
            new_data.push(*d);
            new_owners.push(*o);
        }

        (new_data, new_owners)
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
        points_per_proc: (Vec<usize>, Vec<usize>),
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
        cells_per_proc: &(Vec<usize>, Vec<usize>),
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

        let total_cell_count = cells_per_proc.0.iter().sum();
        // Need also extended counts that have as last element the total cell count.
        let mut extended_counts = cells_per_proc.0;
        extended_counts.push(total_cell_count);

        let mut cell_types = Vec::with_capacity(total_cell_count);
        let mut cell_degrees = Vec::with_capacity(total_cell_count);
        let mut new_cell_owners = Vec::with_capacity(total_cell_count);

        cells_per_proc.1.iter().for_each(|c| {
            cell_types.push(self.cell_type(*c));
            cell_degrees.push(self.cell_degree(*c));
            new_cell_owners.push(cell_owners[*c]);
        });

        let cell_owners = new_cell_owners;

        // Cell points are a bit more complicated as there may be different numbers of points per cell.
        // Hence, have to adapt the corresponding counts.

        let mut cell_points = Vec::<usize>::default();
        let mut cell_point_count = Vec::<usize>::default();

        for (first, last) in extended_counts.iter().tuple_windows() {
            let mut current_count = 0;
            for cell in cells_per_proc.1[*first..*last].iter() {
                let pts = self.cell_points(*cell);
                cell_points.extend_from_slice(pts);
                current_count += pts.len();
            }
            cell_point_count.push(current_count);
        }

        // Have to send the cell data to the other processes now.
        let cell_types = scatterv_root(comm, &cells_per_proc.0, &cell_types);
        let cell_degrees = scatterv_root(comm, &cells_per_proc.0, &cell_degrees);
        let cell_owners = scatterv_root(comm, &cells_per_proc.0, &cell_owners);
        let cell_points = scatterv_root(comm, &cell_point_count, &cell_points);

        // TODO! Fix up the order below.
        self.reorder_cells(
            rank,
            &cell_owners,
            &cell_points,
            &cell_types,
            &cell_degrees,
            &cell_owners,
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
        use coupe::{KMeans, Partition, Point3D};

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
    fn assign_vertex_owners(&self, nprocesses: usize, cell_owners: &[usize]) -> Vec<usize> {
        // Initialise empty array with number of processes as value and length same as number of points.
        let mut vertex_owners = vec![nprocesses; self.point_count()];

        // Each vertex is owned by the minimum process that owns a cell that contains the vertex.
        for (i, owner) in cell_owners.iter().enumerate() {
            for v in self.cell_vertices(i) {
                vertex_owners[*v] = std::cmp::min(vertex_owners[*v], *owner);
            }
        }

        // Return the vertex owners. Note that the array contains garbage information
        // for the points of the grids that are not vertices.
        vertex_owners
    }

    /// Compute the vertices and cells that each process needs access to.
    /// The function returns three tuples, vertices, points, and cells.
    /// The first vector of each tuple contains the number of entities on each process
    /// and is of length `nprocesses`. The second tuple contains a flattened array
    /// of all the entity indices aligned according to the counts in the first array.
    #[allow(clippy::type_complexity)]
    fn get_vertices_points_and_cells(
        &self,
        nprocesses: usize,
        cell_owners: &[usize],
    ) -> (ChunkedData<usize>, ChunkedData<usize>, ChunkedData<usize>) {
        // The following maps each vertex index to a set of cells that contain it.
        let mut vertex_to_cells = HashMap::<usize, HashSet<usize>>::new();

        // For each cell, go through its vertices and add the cell to the cell-list of the vertex.
        for i in 0..self.cell_count() {
            for v in self.cell_vertices(i) {
                vertex_to_cells.entry(*v).or_default().insert(i);
            }
        }

        // This instantiates a vector with `nprocesses` empty sets.
        let map_creator =
            || -> Vec<HashSet<usize>> { (0..nprocesses).map(|_| HashSet::new()).collect_vec() };

        // Each of these is a map from process to a set of indices.
        let mut vertices = map_creator();
        let mut points = map_creator();
        let mut cells = map_creator();

        // This assigns cell to processes. It iterates through a cell and associates
        // not only the cell itself to the process that owns it but also adds all the
        // neighboring cells to the same process. This means that processes own not
        // only the original cells from `cell_owners` but also all the neighbouring cells.
        // Also this means that a cell can live on multiple processes.
        for (i, owner) in cell_owners.iter().enumerate() {
            for v in self.cell_vertices(i) {
                // Add all cells that are connected to the vertex to the same process.
                cells[*owner].extend(vertex_to_cells.get(v).expect("Vertex not found."));
            }
        }

        // Now go through and add the points and vertices to the same processes that also the cells are
        // assigned to.
        for (proc, cells) in cells.iter().enumerate() {
            for cell in cells {
                for v in self.cell_vertices(*cell) {
                    vertices[proc].insert(*v);
                }
                for p in self.cell_points(*cell) {
                    points[proc].insert(*p);
                }
            }
        }

        // The following function flattens an array of sets and returns the counts and the flattened array
        // as ChunkedData struct.
        let flatten = |data: Vec<HashSet<usize>>| -> ChunkedData<usize> {
            let mut idx_bounds = data
                .iter()
                .scan(0, |acc, x| {
                    let old = *acc;
                    *acc += x.len();
                    Some(old)
                })
                .collect_vec();
            idx_bounds.push(idx_bounds.last().unwrap() + data.last().unwrap().len());
            let mut tmp = Vec::<usize>::with_capacity(*idx_bounds.last().unwrap());
            for s in data {
                tmp.extend(s.iter());
            }
            ChunkedData {
                data: tmp,
                idx_bounds,
            }
        };

        (flatten(vertices), flatten(points), flatten(cells))
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

fn scatterv_root<T: Equivalence + Copy>(
    comm: &impl Communicator,
    chunks_per_proc: &ChunkedData<T>,
) -> Vec<T> {
    let rank = comm.rank() as usize;
    let size = comm.size() as usize;

    let send_counts = chunks_per_proc
        .counts()
        .iter()
        .map(|&x| x as i32)
        .collect_vec();

    let mut recv_count: i32 = 0;
    let mut recvbuf: Vec<T> = Vec::<T>::with_capacity(send_counts[rank] as usize);
    // This avoids having the pre-initialise the array. We simply transmute the spare capacity
    // into a valid reference and later manually set the length of the array to the full capacity.
    let recvbuf_ref: &mut [T] = unsafe { std::mem::transmute(recvbuf.spare_capacity_mut()) };

    // The idx-bounds are one too long as their last element is the total number of
    // elements. We don't want this for the displacements.
    let displacements = chunks_per_proc.idx_bounds[0..size]
        .iter()
        .map(|&x| x as i32)
        .collect_vec();

    // Now scatter the counts to each process.
    comm.this_process()
        .scatter_into_root(&send_counts, &mut recv_count);

    // We now prepare the send partition of the variable length data.
    let send_partition =
        mpi::datatype::Partition::new(&chunks_per_proc.data, &send_counts[..], &displacements[..]);

    // And now we send the partition.
    comm.this_process()
        .scatter_varcount_into_root(&send_partition, recvbuf_ref);

    unsafe { recvbuf.set_len(send_counts[rank] as usize) };

    recvbuf
}

// Receiev the scattered data from `root`.
fn scatterv<T: Equivalence + Copy>(comm: &impl Communicator, root: usize) -> Vec<T> {
    let mut recv_count: i32 = 0;

    // First we need to receive the number of elements that we are about to get.
    comm.process_at_rank(root as i32)
        .scatter_into(&mut recv_count);

    // We prepare an unitialized buffer to receive the data.
    let mut recvbuf: Vec<T> = Vec::<T>::with_capacity(recv_count as usize);
    // This avoids having the pre-initialise the array. We simply transmute the spare capacity
    // into a valid reference and later manually set the length of the array to the full capacity.
    let recvbuf_ref: &mut [T] = unsafe { std::mem::transmute(recvbuf.spare_capacity_mut()) };

    // And finally we receive the data.
    comm.process_at_rank(root as i32)
        .scatter_varcount_into(recvbuf_ref);

    // Don't forget to manually set the length of the vector to the correct value.
    unsafe { recvbuf.set_len(recv_count as usize) };
    recvbuf
}

// This routine synchronizes indices across processes.
// - indices: The set of all indices.
// - owners: The owning rank of each index.
// The output is a tuple consisting of the associated global indices and the ownership information.
fn synchronize_indices(
    comm: &impl Communicator,
    indices: &[usize],
    owners: &[usize],
) -> (Vec<usize>, Vec<Ownership>) {
    let rank = comm.rank() as usize;
    let size = comm.size() as usize;

    // First count the local indices.
    let number_of_local_indices = izip!(indices, owners).filter(|(_, &o)| o == rank).count();
    let mut current_global_index: usize = 0;

    // We now do an exclusive scan to get the right start point for the number of global indices.

    comm.exclusive_scan_into(
        &number_of_local_indices,
        &mut current_global_index,
        SystemOperation::sum(),
    );

    // Each process now has the start offset for its global indices. We now iterate through
    // the indices. We assign a global index if the index is owned by the current process.
    // If not we save the index, its owning process and the local index into the `indices` array.
    let mut global_indices = vec![0; indices.len()];

    let mut ghosts = Vec::<usize>::default();
    let mut local_indices_of_ghosts = Vec::<usize>::default();
    let mut owners_of_ghosts = Vec::<usize>::default();

    for (pos, (index, owner)) in izip!(indices, owners).enumerate() {
        if *owner == rank {
            global_indices[pos] = current_global_index;
            current_global_index += 1;
        } else {
            ghosts.push(*index);
            local_indices_of_ghosts.push(pos);
            owners_of_ghosts.push(*owner);
        }
    }

    // We want to resort the ghost by process. This makes the communication a lot easier.
    let sort_indices = (0..ghosts.len())
        .sorted_by_key(|&i| owners_of_ghosts[i])
        .collect_vec();
    let ghosts = sort_indices.iter().map(|&i| ghosts[i]).collect::<Vec<_>>();
    let local_indices_of_ghosts = sort_indices
        .iter()
        .map(|&i| local_indices_of_ghosts[i])
        .collect::<Vec<_>>();
    let owners_of_ghosts = sort_indices
        .iter()
        .map(|&i| owners_of_ghosts[i])
        .collect::<Vec<_>>();

    // Now we have to send the ghost indices to the owning process. For this
    // each process first needs to know how many global indices it gets.

    let mut counts = vec![0 as usize; size];
    for o in owners_of_ghosts.iter() {
        counts[*o] += 1;
    }

    // Now send the counts around via an all-to-all communication.

    let (recv_counts, recv_data) = all_to_all_varcount(comm, &counts, &ghosts);

    // We now have the ghosts on the owning processes. We iterate through and assign the corresponding local indices
    // and send those back to the original processes.
    // The difficulty is to assign the local index to each received ghost index. We need a map from indices to local
    // indices. That is best done with a HashMap.
    let send_back_local_indices = {
        let mut map = HashMap::<usize, usize>::new();
        for (pos, &index) in indices.iter().enumerate() {
            map.insert(index, pos);
        }

        recv_data.iter().map(|i| *map.get(i).unwrap()).collect_vec()
    };

    // Now we can actually send the indices back. The nice thing is that the recv_counts turn into send_counts.

    let (_, remote_local_indices_of_ghosts) =
        all_to_all_varcount(comm, &recv_counts, &send_back_local_indices);

    // We now have to do exactly the same thing with the global indices so that the original process knows the
    // global indices of its ghosts.

    let send_back_global_indices = {
        let mut map = HashMap::<usize, usize>::new();
        for (pos, &index) in indices.iter().enumerate() {
            map.insert(index, global_indices[pos]);
        }

        recv_data.iter().map(|i| *map.get(i).unwrap()).collect_vec()
    };

    let (_, remote_global_indices_of_ghosts) =
        all_to_all_varcount(comm, &recv_counts, &send_back_global_indices);

    // We can now iterate through the ghosts and assign the correct global indices.
    for (pos, ghost_local_index) in local_indices_of_ghosts.iter().enumerate() {
        global_indices[*ghost_local_index] = remote_global_indices_of_ghosts[pos];
    }

    // Finally, we setup the ownership information

    let mut ownership = vec![Ownership::Owned; indices.len()];

    for (pos, ghost_local_index) in local_indices_of_ghosts.iter().enumerate() {
        ownership[*ghost_local_index] =
            Ownership::Ghost(owners_of_ghosts[pos], remote_local_indices_of_ghosts[pos]);
    }

    (global_indices, ownership)
}

// Synchronize the entity ids across processes.
// - entity_length: The number of vertices in each entity.
// - entities: An array of size `number_of_entities * entity_length` containing all the global vertex ids for all entities.
// - owners: An array of size `number_of_entities` containing the owning process of each entity.
fn synchronize_entities(
    comm: &impl Communicator,
    entity_length: usize,
    entities: Vec<usize>,
    owners: Vec<usize>,
) -> (Vec<usize>, Vec<Ownership>) {
    todo!()
}

// Performs an all-to-all communication.
// Returns the receive counts from each processor
// and the received data.
fn all_to_all_varcount<T: Equivalence>(
    comm: &impl Communicator,
    counts: &[usize],
    data: &[T],
) -> (Vec<usize>, Vec<T>) {
    // We need the counts as i32 types.

    let counts = counts.iter().map(|&x| x as i32).collect_vec();

    // First send around the counts via an all-to-all
    let mut recv_counts = vec![0 as i32; comm.size() as usize];
    comm.all_to_all_into(&counts, &mut recv_counts);

    // Now we can prepare the actual data. We have to allocate the data and compute the send partition and the receive partition.

    let mut receive_data = Vec::<T>::with_capacity(recv_counts.iter().sum::<i32>() as usize);
    let receive_buf: &mut [T] = unsafe { std::mem::transmute(receive_data.spare_capacity_mut()) };

    let send_displacements = counts
        .iter()
        .scan(0, |acc, &x| {
            let old = *acc;
            *acc += x;
            Some(old)
        })
        .collect_vec();

    let receive_displacements = recv_counts
        .iter()
        .scan(0, |acc, &x| {
            let old = *acc;
            *acc += x;
            Some(old)
        })
        .collect_vec();

    let send_partition = mpi::datatype::Partition::new(data, counts, send_displacements);
    let mut receive_partition =
        mpi::datatype::PartitionMut::new(receive_buf, &recv_counts[..], receive_displacements);

    comm.all_to_all_varcount_into(&send_partition, &mut receive_partition);

    unsafe { receive_data.set_len(recv_counts.iter().sum::<i32>() as usize) };

    (
        recv_counts.iter().map(|i| *i as usize).collect_vec(),
        receive_data,
    )
}
