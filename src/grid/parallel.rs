//! A parallel implementation of a grid
use crate::traits::{Grid, Topology, Builder, ParallelBuilder as ParallelBuilderTrait, ParallelGrid as ParallelGridTrait};
use crate::types::{CellLocalIndexPair, Ownership, Array2D};
use mpi::{
    request::{LocalScope, WaitGuard},
    topology::{Communicator, Process},
    traits::{Buffer, Destination, Equivalence, Source},
};
use ndelement::reference_cell;
use rlst::{rlst_dynamic_array2, Array, BaseArray, RandomAccessMut, VectorContainer};
use std::collections::HashMap;
use std::marker::PhantomData;

/// Grid local to a process
pub struct LocalGrid<G: Grid> {
    serial_grid: G,
    ownership: Vec<Vec<Ownership>>,
}

pub struct ParallelTopology<'a, T: Topology> {
    topology: &'a T,
}


impl<G: Grid> Grid for LocalGrid<G> {
    type T = G::T;
    type Entity<'a> = G::Entity<'a> where Self: 'a;
    type GeometryMap<'a> = G::GeometryMap<'a> where Self: 'a;
    type EntityDescriptor = G::EntityDescriptor;
    type EntityIter<'a> = G::EntityIter<'a> where Self: 'a;

    fn geometry_dim(&self) -> usize {
        self.serial_grid.geometry_dim()
    }
    fn topology_dim(&self) -> usize {
        self.serial_grid.topology_dim()
    }

    fn entity(&self, dim: usize, local_index: usize) -> Option<Self::Entity<'_>> {
        self.serial_grid.entity(dim, local_index)
    }

    fn entity_types(&self, dim: usize) -> &[Self::EntityDescriptor] {
        self.serial_grid.entity_types(dim)
    }

    fn entity_count(&self, entity_type: Self::EntityDescriptor) -> usize {
        self.serial_grid.entity_count(entity_type)
    }

    fn entity_iter(&self, dim: usize) -> Self::EntityIter<'_> {
        self.serial_grid.entity_iter(dim)
    }

    fn entity_from_id(&self, dim: usize, id: usize) -> Option<Self::Entity<'_>> {
        self.serial_grid.entity_from_id(dim, id)
    }

    fn geometry_map(&self, entity_type: Self::EntityDescriptor, points: &[Self::T]) -> Self::GeometryMap<'_> {
        self.serial_grid.geometry_map(entity_type, points)
    }
}

/// Parallel grid
pub struct ParallelGrid<'comm, C: Communicator, G: Grid> {
    pub(crate) comm: &'comm C,
    pub(crate) local_grid: LocalGrid<G>,
}

impl<'comm, C: Communicator, G: Grid> ParallelGrid<'comm, C, G> {
    /// Create new parallel grid
    pub fn new(
        comm: &'comm C,
        serial_grid: G,
        vertex_owners: HashMap<usize, usize>,
        edge_owners: HashMap<usize, usize>,
        cell_owners: HashMap<usize, usize>,
    ) -> Self {
        let rank = comm.rank() as usize;
        let size = comm.size() as usize;

        // Create cell ownership
        let mut cell_ownership = HashMap::new();
        let mut cells_to_query = vec![vec![]; size];

        for (id, owner) in cell_owners.iter() {
            if *owner != rank {
                cells_to_query[*owner].push(*id);
            }
        }
        mpi::request::scope(|scope| {
            for (p, cq) in cells_to_query.iter().enumerate() {
                if p != rank {
                    let _ =
                        WaitGuard::from(comm.process_at_rank(p as i32).immediate_send(scope, cq));
                }
            }
        });
        for p in 0..size {
            if p != rank {
                let (cells_queried, _status) =
                    comm.process_at_rank(p as i32).receive_vec::<usize>();
                let send_back =
                    cells_queried
                        .iter()
                        .map(|id| {
                            serial_grid.topology().face_index_to_flat_index(
                                Topology::cell_id_to_index(serial_grid.topology(), *id),
                            )
                        })
                        .collect::<Vec<_>>();
                mpi::request::scope(|scope| {
                    let _ = WaitGuard::from(
                        comm.process_at_rank(p as i32)
                            .immediate_send(scope, &send_back),
                    );
                });
            }
        }
        let mut cell_info = vec![vec![]; size];
        for (p, ci) in cell_info.iter_mut().enumerate() {
            if p != rank {
                (*ci, _) = comm.process_at_rank(p as i32).receive_vec::<usize>();
            }
        }

        let mut indices = vec![0; size];
        for (id, owner) in cell_owners.iter() {
            cell_ownership.insert(
                Topology::cell_id_to_index(serial_grid.topology(), *id),
                if *owner == rank {
                    Ownership::Owned
                } else {
                    indices[*owner] += 1;
                    Ownership::Ghost(*owner, cell_info[*owner][indices[*owner] - 1])
                },
            );
        }

        // Create vertex ownership
        let mut vertex_ownership = vec![Ownership::Owned; vertex_owners.len()];
        let mut vertices_to_query = vec![vec![]; size];

        for (id, owner) in vertex_owners.iter() {
            if *owner != rank {
                vertices_to_query[*owner].push(*id);
            }
        }
        mpi::request::scope(|scope| {
            for (p, vq) in vertices_to_query.iter().enumerate() {
                if p != rank {
                    let _ =
                        WaitGuard::from(comm.process_at_rank(p as i32).immediate_send(scope, vq));
                }
            }
        });
        for p in 0..size {
            if p != rank {
                let (vertices_queried, _status) =
                    comm.process_at_rank(p as i32).receive_vec::<usize>();
                let send_back = vertices_queried
                    .iter()
                    .map(|id| Topology::vertex_id_to_index(serial_grid.topology(), *id))
                    .collect::<Vec<_>>();
                mpi::request::scope(|scope| {
                    let _ = WaitGuard::from(
                        comm.process_at_rank(p as i32)
                            .immediate_send(scope, &send_back),
                    );
                });
            }
        }
        let mut vertex_info = vec![vec![]; size];
        for (p, vi) in vertex_info.iter_mut().enumerate() {
            if p != rank {
                (*vi, _) = comm.process_at_rank(p as i32).receive_vec::<usize>();
            }
        }

        let mut indices = vec![0; size];
        for (id, owner) in vertex_owners.iter() {
            vertex_ownership[Topology::vertex_id_to_index(serial_grid.topology(), *id)] =
                if *owner == rank {
                    Ownership::Owned
                } else {
                    indices[*owner] += 1;
                    Ownership::Ghost(*owner, vertex_info[*owner][indices[*owner] - 1])
                };
        }

        // Create edge ownership
        let mut edge_ownership = vec![Ownership::Owned; edge_owners.len()];
        let mut edges_to_query = vec![vec![]; size];

        for (id, owner) in edge_owners.iter() {
            if *owner != rank {
                edges_to_query[*owner].push(*id);
            }
        }
        mpi::request::scope(|scope| {
            for (p, eq) in edges_to_query.iter().enumerate() {
                if p != rank {
                    let _ =
                        WaitGuard::from(comm.process_at_rank(p as i32).immediate_send(scope, eq));
                }
            }
        });
        for p in 0..size {
            if p != rank {
                let (edges_queried, _status) =
                    comm.process_at_rank(p as i32).receive_vec::<usize>();
                let send_back = edges_queried
                    .iter()
                    .map(|id| Topology::edge_id_to_index(serial_grid.topology(), *id))
                    .collect::<Vec<_>>();
                mpi::request::scope(|scope| {
                    let _ = WaitGuard::from(
                        comm.process_at_rank(p as i32)
                            .immediate_send(scope, &send_back),
                    );
                });
            }
        }
        let mut edge_info = vec![vec![]; size];
        for (p, ei) in edge_info.iter_mut().enumerate() {
            if p != rank {
                (*ei, _) = comm.process_at_rank(p as i32).receive_vec::<usize>();
            }
        }

        let mut indices = vec![0; size];
        for (id, owner) in edge_owners.iter() {
            edge_ownership[Topology::edge_id_to_index(serial_grid.topology(), *id)] =
                if *owner == rank {
                    Ownership::Owned
                } else {
                    indices[*owner] += 1;
                    Ownership::Ghost(*owner, edge_info[*owner][indices[*owner] - 1])
                };
        }

        let local_grid = LocalGrid {
            serial_grid,
            vertex_ownership,
            edge_ownership,
            cell_ownership,
        };
        Self { comm, local_grid }
    }
}

impl<'comm, C: Communicator, G: Grid> Grid for ParallelGrid<'comm, C, G> {
    type T = G::T;
    type Entity<'a> = G::Entity<'a> where Self: 'a;
    type GeometryMap<'a> = G::GeometryMap<'a> where Self: 'a;
    type EntityDescriptor = G::EntityDescriptor;
    type EntityIter<'a> = G::EntityIter<'a> where Self: 'a;

    fn geometry_dim(&self) -> usize {
        self.local_grid.geometry_dim()
    }
    fn topology_dim(&self) -> usize {
        self.local_grid.topology_dim()
    }

    fn entity(&self, dim: usize, local_index: usize) -> Option<Self::Entity<'_>> {
        self.local_grid.entity(dim, local_index)
    }

    fn entity_types(&self, dim: usize) -> &[Self::EntityDescriptor] {
        self.local_grid.entity_types(dim)
    }

    fn entity_count(&self, entity_type: Self::EntityDescriptor) -> usize {
        self.local_grid.entity_count(entity_type)
    }

    fn entity_iter(&self, dim: usize) -> Self::EntityIter<'_> {
        self.local_grid.entity_iter(dim)
    }

    fn entity_from_id(&self, dim: usize, id: usize) -> Option<Self::Entity<'_>> {
        self.local_grid.entity_from_id(dim, id)
    }

    fn geometry_map(&self, entity_type: Self::EntityDescriptor, points: &[Self::T]) -> Self::GeometryMap<'_> {
        self.local_grid.geometry_map(entity_type, points)
    }
}
/*
struct ParallelBuilder<G: Grid> {
    grid: PhantomData<G>
}

impl<G: Grid> ParallelBuilderTrait for ParallelBuilder<G>
where
    Vec<G::T>: Buffer,
    G::T: Equivalence,
{
    type ParallelGrid<'a, C: Communicator + 'a> = ParallelGrid<'a, C, G>;

    fn create_parallel_grid<'a, C: Communicator>(
        self,
        comm: &'a C,
        cell_owners: &HashMap<usize, usize>,
    ) -> Self::ParallelGrid<'a, C> {
        let rank = comm.rank() as usize;
        let size = comm.size() as usize;

        let npts = self.point_indices_to_ids().len();
        let ncells = self.cell_indices_to_ids().len();

        let gdim = 3; // TODO

        // data used in computation
        let mut vertex_owners = vec![(-1, 0); npts];
        let mut vertex_counts = vec![0; size];
        let mut cell_indices_per_proc = vec![vec![]; size];
        let mut vertex_indices_per_proc = vec![vec![]; size];
        let mut point_indices_per_proc = vec![vec![]; size];
        let mut edge_owners = HashMap::new();
        let mut edge_ids = HashMap::new();
        let mut edge_counts = vec![0; size];
        let mut edges_included_per_proc = vec![vec![]; size];
        let mut edge_id = 0;

        // data to send to other processes
        let mut points_per_proc = vec![vec![]; size];
        let mut point_ids_per_proc = vec![vec![]; size];
        let mut cells_per_proc = vec![vec![]; size];
        let mut cell_owners_per_proc = vec![vec![]; size];
        let mut cell_ids_per_proc = vec![vec![]; size];
        let mut vertex_owners_per_proc = vec![vec![]; size];
        let mut edges_per_proc = vec![vec![]; size];
        let mut edge_owners_per_proc = vec![vec![]; size];
        let mut edge_ids_per_proc = vec![vec![]; size];
        let mut extra_cell_info_per_proc = vec![self.new_extra_cell_info(); size];

        for (index, id) in self.cell_indices_to_ids().iter().enumerate() {
            let owner = cell_owners[id];
            for pt in self.cell_points(index) {
                if !point_indices_per_proc[owner].contains(pt) {
                    point_indices_per_proc[owner].push(*pt);
                    for i in 0..gdim {
                        points_per_proc[owner].push(self.points()[pt * gdim + i])
                    }
                    point_ids_per_proc[owner].push(self.point_indices_to_ids()[*pt]);
                }
            }
            for v in self.cell_vertices(index) {
                if vertex_owners[*v].0 == -1 {
                    vertex_owners[*v] = (owner as i32, vertex_counts[owner]);
                }
                if !vertex_indices_per_proc[owner].contains(v) {
                    vertex_indices_per_proc[owner].push(*v);
                    vertex_owners_per_proc[owner].push(vertex_owners[*v].0 as usize);
                    vertex_counts[owner] += 1;
                }
            }
        }

        for (index, id) in self.cell_indices_to_ids().iter().enumerate() {
            let ref_conn = &reference_cell::connectivity(self.cell_type(index))[1];
            let owner = cell_owners[id];
            for e in ref_conn {
                let cell = self.cell_vertices(index);
                let mut v0 = cell[e[0][0]];
                let mut v1 = cell[e[0][1]];
                if v0 > v1 {
                    std::mem::swap(&mut v0, &mut v1);
                }
                if edge_owners.get_mut(&(v0, v1)).is_none() {
                    edge_owners.insert((v0, v1), (owner, edge_counts[owner]));
                    edge_ids.insert((v0, v1), edge_id);
                    edge_id += 1;
                    edges_included_per_proc[owner].push((v0, v1));
                    edges_per_proc[owner].push(v0);
                    edges_per_proc[owner].push(v1);
                    edge_owners_per_proc[owner].push(edge_owners[&(v0, v1)].0);
                    edge_ids_per_proc[owner].push(edge_ids[&(v0, v1)]);
                    edge_counts[owner] += 1;
                }
            }
        }

        for index in 0..ncells {
            for p in 0..size {
                for v in self.cell_points(index) {
                    if vertex_indices_per_proc[p].contains(v) {
                        cell_indices_per_proc[p].push(index);
                        break;
                    }
                }
            }
        }

        for p in 0..size {
            for index in &cell_indices_per_proc[p] {
                let id = self.cell_indices_to_ids()[*index];
                for pt in self.cell_points(*index) {
                    if !point_indices_per_proc[p].contains(pt) {
                        point_indices_per_proc[p].push(*pt);
                        for i in 0..gdim {
                            points_per_proc[p].push(self.points()[pt * gdim + i]);
                        }
                        point_ids_per_proc[p].push(self.point_indices_to_ids()[*pt]);
                    }
                }
                for v in self.cell_points(*index) {
                    if !vertex_indices_per_proc[p].contains(v) {
                        vertex_indices_per_proc[p].push(*v);
                        vertex_owners_per_proc[p].push(vertex_owners[*v].0 as usize);
                    }
                    cells_per_proc[p].push(
                        vertex_indices_per_proc[p]
                            .iter()
                            .position(|&r| r == *v)
                            .unwrap(),
                    );
                }
                let ref_conn = &reference_cell::connectivity(self.cell_type(*index))[1];

                for e in ref_conn {
                    let cell = self.cell_vertices(*index);
                    let mut v0 = cell[e[0][0]];
                    let mut v1 = cell[e[0][1]];
                    if v0 > v1 {
                        std::mem::swap(&mut v0, &mut v1);
                    }
                    if !edges_included_per_proc[p].contains(&(v0, v1)) {
                        edges_included_per_proc[p].push((v0, v1));
                        edges_per_proc[p].push(v0);
                        edges_per_proc[p].push(v1);
                        edge_owners_per_proc[p].push(edge_owners[&(v0, v1)].0);
                        edge_ids_per_proc[p].push(edge_ids[&(v0, v1)]);
                    }
                }

                cell_ids_per_proc[p].push(id);
                cell_owners_per_proc[p].push(cell_owners[&id]);
                self.push_extra_cell_info(&mut extra_cell_info_per_proc[p], id);
            }
        }

        mpi::request::scope(|scope| {
            for p in 1..size {
                let process = comm.process_at_rank(p as i32);
                let _ = WaitGuard::from(process.immediate_send(scope, &points_per_proc[p]));
                let _ = WaitGuard::from(process.immediate_send(scope, &point_ids_per_proc[p]));
                let _ = WaitGuard::from(process.immediate_send(scope, &cells_per_proc[p]));
                let _ = WaitGuard::from(process.immediate_send(scope, &cell_owners_per_proc[p]));
                let _ = WaitGuard::from(process.immediate_send(scope, &cell_ids_per_proc[p]));
                let _ = WaitGuard::from(process.immediate_send(scope, &vertex_owners_per_proc[p]));
                let _ = WaitGuard::from(process.immediate_send(scope, &edges_per_proc[p]));
                let _ = WaitGuard::from(process.immediate_send(scope, &edge_owners_per_proc[p]));
                let _ = WaitGuard::from(process.immediate_send(scope, &edge_ids_per_proc[p]));
                self.send_extra_cell_info(scope, &process, &extra_cell_info_per_proc[p]);
            }
        });

        self.create_internal(
            comm,
            &points_per_proc[rank],
            &point_ids_per_proc[rank],
            &cells_per_proc[rank],
            &cell_owners_per_proc[rank],
            &cell_ids_per_proc[rank],
            &vertex_owners_per_proc[rank],
            &point_ids_per_proc[rank],
            &edges_per_proc[rank],
            &edge_owners_per_proc[rank],
            &edge_ids_per_proc[rank],
            &extra_cell_info_per_proc[rank],
        )
    }

    fn receive_parallel_grid<C: Communicator>(
        self,
        comm: &C,
        root_rank: usize,
    ) -> ParallelGrid<'_, C, G> {
        let root_process = comm.process_at_rank(root_rank as i32);

        let (points, _status) =
            root_process.receive_vec::<Self::G::T>();
        let (point_ids, _status) = root_process.receive_vec::<usize>();
        let (cells, _status) = root_process.receive_vec::<usize>();
        let (cell_owners, _status) = root_process.receive_vec::<usize>();
        let (cell_ids, _status) = root_process.receive_vec::<usize>();
        let (vertex_owners, _status) = root_process.receive_vec::<usize>();
        let (edges, _status) = root_process.receive_vec::<usize>();
        let (edge_owners, _status) = root_process.receive_vec::<usize>();
        let (edge_ids, _status) = root_process.receive_vec::<usize>();
        let mut extra_cell_info = self.new_extra_cell_info();
        self.receive_extra_cell_info(&root_process, &mut extra_cell_info);

        self.create_internal(
            comm,
            &points,
            &point_ids,
            &cells,
            &cell_owners,
            &cell_ids,
            &vertex_owners,
            &point_ids,
            &edges,
            &edge_owners,
            &edge_ids,
            &extra_cell_info,
        )
    }
}
*/
