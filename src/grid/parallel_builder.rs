//! Parallel grid builder

use super::ParallelGridImpl;
use crate::{
    traits::{
        Builder, Entity, GeometryBuilder, Grid, GridBuilder, ParallelBuilder, Topology,
        TopologyBuilder,
    },
    types::{GraphPartitioner, Ownership},
};
use itertools::{Itertools, izip};
use mpi::{
    collective::SystemOperation,
    traits::{Buffer, Communicator, CommunicatorCollectives, Equivalence, Root},
};
use std::collections::{HashMap, HashSet};

#[cfg(feature = "coupe")]
use coupe::{KMeans, Partition, Point3D};

#[cfg(feature = "scotch")]
use scotch::{Architecture, Graph, Strategy, graph};

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
        = ParallelGridImpl<'a, C, B::Grid>
    where
        Self: 'a;
    fn create_parallel_grid_root<'a, C: Communicator>(
        &self,
        comm: &'a C,
        partitioner: GraphPartitioner,
    ) -> ParallelGridImpl<'a, C, B::Grid> {
        // If we have only a single process we can just create a serial grid.
        if comm.size() == 1 {
            let serial_grid = self.create_grid();

            // Need to fill ownership and global indices now.

            let mut owners = HashMap::new();
            let mut global_indices = HashMap::new();

            for dim in 0..self.tdim() {
                for t in serial_grid.entity_types(dim) {
                    let entity_count = serial_grid.entity_iter(*t).count();
                    owners.insert(*t, vec![Ownership::Owned; entity_count]);
                    global_indices.insert(*t, (0..entity_count).collect_vec());
                }
            }

            return ParallelGridImpl::new(comm, serial_grid, owners, global_indices);
        }
        // Partition the cells
        let cell_owners = self.partition_cells(comm.size() as usize, partitioner);

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
    ) -> ParallelGridImpl<'a, C, B::Grid> {
        // First we receive all the data.
        let root_rank = root_rank as usize;
        let point_indices = scatterv(comm, root_rank);
        let coordinates = scatterv(comm, root_rank);
        let vertex_indices = scatterv(comm, root_rank);
        let vertex_owners = scatterv(comm, root_rank);
        let cell_indices = scatterv(comm, root_rank);
        let cell_points = scatterv(comm, root_rank);
        let cell_types = scatterv(comm, root_rank);
        let cell_degrees = scatterv(comm, root_rank);
        let cell_owners = scatterv(comm, root_rank);

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
    ) -> ParallelGridImpl<'a, C, Self::Grid>
    where
        Self::Grid: Sync,
    {
        let rank = comm.rank() as usize;
        // First we want to reorder everything so that owned data comes first in the arrays.

        // Reorder the vertices.
        let (vertex_indices, vertex_owners) = {
            let mut new_indices = Vec::<usize>::with_capacity(vertex_indices.len());
            let mut new_owners = Vec::<usize>::with_capacity(vertex_owners.len());

            for (&v, &o) in izip!(&vertex_indices, &vertex_owners).filter(|&(_, &o)| o == rank) {
                new_indices.push(v);
                new_owners.push(o);
            }

            for (&v, &o) in izip!(&vertex_indices, &vertex_owners).filter(|&(_, &o)| o != rank) {
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
            for (o, ct, cd) in izip!(cell_owners.iter(), cell_types.iter(), cell_degrees.iter()) {
                let npts = self.npts(*ct, *cd);
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
                    .filter(|&(_, _, _, &o)| o == rank)
            {
                new_indices.push(*cell_index);
                new_types.push(*cell_type);
                new_degrees.push(*cell_degree);
                new_owners.push(*cell_owner);
            }

            for (cell_index, cell_type, cell_degree, cell_owner) in
                izip!(&cell_indices, &cell_types, &cell_degrees, &cell_owners,)
                    .filter(|&(_, _, _, &o)| o != rank)
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
                .filter(|&(_, &o)| o == rank)
                .map(|(p, _)| p)
            {
                new_points.push(p);
            }

            for &p in izip!(&cell_points, &cell_point_owners)
                .filter(|&(_, &o)| o != rank)
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
            synchronize_entities(comm, 1, &vertex_indices, &vertex_owners);
        let (cell_global_indices, cell_owners) =
            synchronize_entities(comm, 1, &cell_indices, &cell_owners);

        // For later we need a map from vertex global indices to vertex owners.
        let mut vertex_global_indices_to_owners = HashMap::<usize, usize>::new();
        for (global_index, owner) in izip!(&vertex_global_indices, &vertex_owners) {
            let owner_rank = match *owner {
                Ownership::Owned => rank,
                Ownership::Ghost(p, _) => p,
                _ => panic!("Unsupported ownership: {owner:?}"),
            };
            vertex_global_indices_to_owners.insert(*global_index, owner_rank);
        }

        let point_type = serial_grid.entity_types(0)[0];
        let mut owners = HashMap::new();
        owners.insert(point_type, vertex_owners);
        let mut global_indices = HashMap::new();
        global_indices.insert(point_type, vertex_global_indices.clone());
        for dim in 1..self.tdim() {
            for t in serial_grid.entity_types(dim) {
                // We now iterate through the entities of varying dimensions and assign global indices and ownership.
                // First we need to get out the vertices assigned with each entity and figure out the ownership of the entity.
                // Ownership is determined by the minimum process that owns one of the vertices of the entity.

                // We get out the chunk length by probing the first entity. This only works because the grid has a single
                // element type.
                let chunk_length = serial_grid
                    .entity(*t, 0)
                    .unwrap()
                    .topology()
                    .sub_entity_iter(point_type)
                    .count();

                // We have the chunk length. So iterate through to get all the vertices.
                let number_of_entities = serial_grid.entity_iter(*t).count();
                let mut entities = Vec::with_capacity(number_of_entities * chunk_length);
                let mut entity_ranks = Vec::with_capacity(number_of_entities * chunk_length);
                for entity in serial_grid.entity_iter(*t) {
                    // We iterate through the vertices of the entity and get the global indices of the vertices.
                    // This works because the topology returns positions into the vertex array.
                    // Hence, we can use those indices to map to the global indices.
                    // We also sort everything since later we use the sorted global indices as keys in a hash map
                    // within the `synchronize_entities` function.
                    let vertices = entity
                        .topology()
                        .sub_entity_iter(point_type)
                        .map(|v| vertex_global_indices[v])
                        .sorted()
                        .collect_vec();
                    let owner = *vertices
                        .iter()
                        .map(|v| vertex_global_indices_to_owners.get(v).unwrap())
                        .min()
                        .unwrap();
                    entities.extend(vertices);
                    entity_ranks.push(owner);
                }

                // We now have the entities and their owners. We can synchronize the entities to get the global
                // indices and ownership information.

                let (entity_global_indices, entity_owners) =
                    synchronize_entities(comm, chunk_length, &entities, &entity_ranks);

                owners.insert(*t, entity_owners);
                global_indices.insert(*t, entity_global_indices);
            }
        }
        for (cell_type, owner, index) in izip!(cell_types, cell_owners, cell_global_indices) {
            owners.entry(cell_type).or_insert(vec![]).push(owner);
            global_indices
                .entry(cell_type)
                .or_insert(vec![])
                .push(index);
        }

        ParallelGridImpl::new(comm, serial_grid, owners, global_indices)
    }

    /// Partition the cells
    fn partition_cells(&self, nprocesses: usize, partitioner: GraphPartitioner) -> Vec<usize> {
        // Create an initial partitioning that roughly assigns the same
        // number of cells to each process.
        let mut partition = vec![];
        for i in 0..nprocesses {
            while partition.len() < self.cell_count() * (i + 1) / nprocesses {
                partition.push(i);
            }
        }

        match partitioner {
            GraphPartitioner::None => partition,
            GraphPartitioner::Manual(p) => p,
            #[cfg(feature = "coupe")]
            GraphPartitioner::Coupe => {
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
            #[cfg(feature = "scotch")]
            GraphPartitioner::Scotch => {
                let mut vertex_to_cells = HashMap::<usize, HashSet<i32>>::new();
                for i in 0..self.cell_count() {
                    for v in self.cell_vertices(i) {
                        vertex_to_cells.entry(*v).or_default().insert(i as i32);
                    }
                }

                let mut endpoints = vec![0];
                let mut connectivity = vec![];
                for i in 0..self.cell_count() {
                    let end = endpoints[endpoints.len() - 1] as usize;
                    for v in self.cell_vertices(i) {
                        for c in vertex_to_cells.get(v).unwrap() {
                            if !connectivity[end..].contains(c) {
                                connectivity.push(*c)
                            }
                        }
                    }
                    endpoints.push(connectivity.len() as i32);
                }

                let mut strategy = Strategy::new(); // use the default strategy.
                let arch = Architecture::complete(nprocesses as i32);

                let data = graph::Data::new(
                    0,
                    &endpoints[..endpoints.len() - 1],
                    &endpoints[1..],
                    &[],
                    &[],
                    &connectivity,
                    &[],
                );
                let mut graph = Graph::build(&data).expect("Could not build graph");
                let (vertnbr, _edgenbr) = graph.size();
                let mut parttab = vec![0; vertnbr as usize];
                let _ = graph
                    .mapping(&arch, &mut parttab)
                    .compute(&mut strategy)
                    .expect("");
                parttab.iter().map(|i| *i as usize).collect::<Vec<_>>()
            }
        }
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
                // Each set itself has a random ordering. This makes runs of the method
                // not reproducible. The easiest fix is to sort the set before adding it to `tmp`.
                let mut to_sort = Vec::<usize>::with_capacity(s.len());
                to_sort.extend(s.iter());
                to_sort.sort();
                tmp.extend(to_sort);
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

// Receive the scattered data from `root`.
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
// - Chunk length: An entity can have more than one index. This is the number of indices per entity.
// - indices: The set of all indices.
// - owners: The owning rank of each index.
// The output is a tuple consisting of the associated global indices and the ownership information.
fn synchronize_entities(
    comm: &impl Communicator,
    chunk_length: usize,
    entities: &[usize],
    owners: &[usize],
) -> (Vec<usize>, Vec<Ownership>) {
    let rank = comm.rank() as usize;
    let size = comm.size() as usize;

    // First count the local indices.
    let number_of_local_entities = owners.iter().filter(|&&o| o == rank).count();
    let mut current_global_index: usize = 0;

    // We now do an exclusive scan to get the right start point for the number of global indices.

    comm.exclusive_scan_into(
        &number_of_local_entities,
        &mut current_global_index,
        SystemOperation::sum(),
    );

    // Each process now has the start offset for its global indices. We now iterate through
    // the entities. We assign a global index if the index is owned by the current process.
    // If not we save the entity, its owning process and the local index into the `entities` array.
    let mut global_indices = vec![0; entities.len()];

    let mut ghosts = Vec::<Vec<usize>>::default();
    let mut local_indices_of_ghosts = Vec::<usize>::default();
    let mut owners_of_ghosts = Vec::<usize>::default();

    for (pos, (chunk, owner)) in izip!(&entities.iter().chunks(chunk_length), owners).enumerate() {
        let c = chunk.copied().collect_vec();
        if *owner == rank {
            global_indices[pos] = current_global_index;
            current_global_index += 1;
        } else {
            ghosts.push(c);
            local_indices_of_ghosts.push(pos);
            owners_of_ghosts.push(*owner);
        }
    }

    // We want to resort the ghost by process. This makes the communication a lot easier.
    let sort_indices = (0..ghosts.len())
        .sorted_by_key(|&i| owners_of_ghosts[i])
        .collect_vec();
    let ghosts = sort_indices
        .iter()
        .map(|&i| ghosts[i].clone())
        .collect::<Vec<_>>();
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

    let mut counts = vec![0; size];
    for o in owners_of_ghosts.iter() {
        counts[*o] += chunk_length;
    }

    // Now send the counts around via an all-to-all communication.

    let (recv_counts, recv_data) = all_to_all_varcount(
        comm,
        &counts,
        &ghosts.iter().flatten().copied().collect_vec(),
    );

    // Turn `recv_data` back into a chunked vector
    let recv_data = recv_data
        .iter()
        .chunks(chunk_length)
        .into_iter()
        .map(|c| c.copied().collect_vec())
        .collect_vec();

    // We now have the ghosts on the owning processes. We iterate through and assign the corresponding local indices
    // and send those back to the original processes.
    // The difficulty is to assign the local index to each received ghost index. We need a map from indices to local
    // indices. That is best done with a HashMap.
    let send_back_local_indices = {
        let mut map = HashMap::<Vec<usize>, usize>::new();
        for (pos, chunk) in entities.iter().chunks(chunk_length).into_iter().enumerate() {
            let c = chunk.copied().collect_vec();
            map.insert(c, pos);
        }

        recv_data.iter().map(|i| *map.get(i).unwrap()).collect_vec()
    };

    // Now we can actually send the indices back. We can almost use recv_counts as the new send_counts.
    // The issue is that the recv_counts are multiplied by the chunk length. We have to reverse this.

    let recv_counts = recv_counts.iter().map(|i| *i / chunk_length).collect_vec();

    let (_, remote_local_indices_of_ghosts) =
        all_to_all_varcount(comm, &recv_counts, &send_back_local_indices);

    // We now have to do exactly the same thing with the global indices so that the original process knows the
    // global indices of its ghosts.

    let send_back_global_indices = {
        let mut map = HashMap::<Vec<usize>, usize>::new();
        for (pos, chunk) in entities.iter().chunks(chunk_length).into_iter().enumerate() {
            map.insert(chunk.copied().collect_vec(), global_indices[pos]);
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

    let mut ownership = vec![Ownership::Owned; entities.len()];

    for (pos, ghost_local_index) in local_indices_of_ghosts.iter().enumerate() {
        ownership[*ghost_local_index] =
            Ownership::Ghost(owners_of_ghosts[pos], remote_local_indices_of_ghosts[pos]);
    }

    (global_indices, ownership)
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
    let mut recv_counts = vec![0; comm.size() as usize];
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
