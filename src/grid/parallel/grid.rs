//! Parallel grid
use crate::{
    traits::{Entity, Grid, ParallelGrid as ParallelGridTrait},
    types::Ownership,
};
use mpi::traits::Communicator;

/// Parallel grid entity
pub struct ParallelGridEntity<E: Entity> {
    serial_entity: E,
    ownership: Ownership,
    global_index: usize,
}

impl<E: Entity> ParallelGridEntity<E> {
    /// Create new
    pub fn new(serial_entity: E, ownership: &[Ownership], global_indices: &[usize]) -> Self {
        let index = serial_entity.local_index();
        Self {
            serial_entity,
            ownership: ownership[index],
            global_index: global_indices[index],
        }
    }
}
impl<E: Entity> Entity for ParallelGridEntity<E> {
    type T = E::T;
    type EntityDescriptor = E::EntityDescriptor;
    type Topology<'a> = E::Topology<'a> where Self: 'a;
    type Geometry<'a> = E::Geometry<'a> where Self: 'a;
    fn entity_type(&self) -> E::EntityDescriptor {
        self.serial_entity.entity_type()
    }
    fn local_index(&self) -> usize {
        self.serial_entity.local_index()
    }
    fn global_index(&self) -> usize {
        self.global_index
    }
    fn geometry(&self) -> Self::Geometry<'_> {
        self.serial_entity.geometry()
    }
    fn topology(&self) -> Self::Topology<'_> {
        self.serial_entity.topology()
    }
    fn ownership(&self) -> Ownership {
        self.ownership
    }
}

/// Parallel grid entity iterator
pub struct ParallelGridEntityIter<'a, E: Entity, EntityIter: Iterator<Item = E>> {
    iter: EntityIter,
    ownership: &'a [Ownership],
    global_indices: &'a [usize],
}

impl<'a, E: Entity, EntityIter: Iterator<Item = E>> ParallelGridEntityIter<'a, E, EntityIter> {
    /// Create new
    pub fn new(iter: EntityIter, ownership: &'a [Ownership], global_indices: &'a [usize]) -> Self {
        Self {
            iter,
            ownership,
            global_indices,
        }
    }
}
impl<'a, E: Entity, EntityIter: Iterator<Item = E>> Iterator
    for ParallelGridEntityIter<'a, E, EntityIter>
{
    type Item = ParallelGridEntity<E>;

    fn next(&mut self) -> Option<ParallelGridEntity<E>> {
        let entity = self.iter.next();
        entity.map(|e| ParallelGridEntity::new(e, self.ownership, self.global_indices))
    }
}

/// Local grid on a process
pub struct LocalGrid<G: Grid> {
    serial_grid: G,
    ownership: Vec<Vec<Ownership>>,
    global_indices: Vec<Vec<usize>>,
}

impl<G: Grid> LocalGrid<G> {
    /// Create new
    pub fn new(
        serial_grid: G,
        ownership: Vec<Vec<Ownership>>,
        global_indices: Vec<Vec<usize>>,
    ) -> Self {
        Self {
            serial_grid,
            ownership,
            global_indices,
        }
    }
}
impl<G: Grid> Grid for LocalGrid<G> {
    type T = G::T;
    type Entity<'a> = ParallelGridEntity<G::Entity<'a>> where Self: 'a;
    type GeometryMap<'a> = G::GeometryMap<'a> where Self: 'a;
    type EntityDescriptor = G::EntityDescriptor;
    type EntityIter<'a> = ParallelGridEntityIter<'a, G::Entity<'a>, G::EntityIter<'a>>
    where
        Self: 'a;

    fn geometry_dim(&self) -> usize {
        self.serial_grid.geometry_dim()
    }
    fn topology_dim(&self) -> usize {
        self.serial_grid.topology_dim()
    }

    fn entity(&self, dim: usize, serial_index: usize) -> Option<Self::Entity<'_>> {
        self.serial_grid
            .entity(dim, serial_index)
            .map(|e| ParallelGridEntity::new(e, &self.ownership[dim], &self.global_indices[dim]))
    }

    fn entity_types(&self, dim: usize) -> &[Self::EntityDescriptor] {
        self.serial_grid.entity_types(dim)
    }

    fn entity_count(&self, entity_type: Self::EntityDescriptor) -> usize {
        self.serial_grid.entity_count(entity_type)
    }

    fn entity_iter(&self, dim: usize) -> Self::EntityIter<'_> {
        ParallelGridEntityIter::new(
            self.serial_grid.entity_iter(dim),
            &self.ownership[dim],
            &self.global_indices[dim],
        )
    }

    fn entity_from_id(&self, dim: usize, id: usize) -> Option<Self::Entity<'_>> {
        self.serial_grid
            .entity_from_id(dim, id)
            .map(|e| ParallelGridEntity::new(e, &self.ownership[dim], &self.global_indices[dim]))
    }

    fn geometry_map(
        &self,
        entity_type: Self::EntityDescriptor,
        points: &[Self::T],
    ) -> Self::GeometryMap<'_> {
        self.serial_grid.geometry_map(entity_type, points)
    }
}

/// Parallel grid
pub struct ParallelGrid<'a, C: Communicator, G: Grid> {
    comm: &'a C,
    local_grid: LocalGrid<G>,
}

impl<'a, C: Communicator, G: Grid> ParallelGrid<'a, C, G> {
    /// Create new
    pub fn new(
        comm: &'a C,
        serial_grid: G,
        ownership: Vec<Vec<Ownership>>,
        global_indices: Vec<Vec<usize>>,
    ) -> Self {
        Self {
            comm,
            local_grid: LocalGrid::new(serial_grid, ownership, global_indices),
        }
    }
}

impl<'g, C: Communicator, G: Grid> ParallelGridTrait<C> for ParallelGrid<'g, C, G> {
    type LocalGrid<'a> = LocalGrid<G> where Self: 'a;
    fn comm(&self) -> &C {
        self.comm
    }
    fn local_grid(&self) -> &LocalGrid<G> {
        &self.local_grid
    }
}
impl<'g, C: Communicator, G: Grid> Grid for ParallelGrid<'g, C, G> {
    type T = G::T;
    type Entity<'a> = <LocalGrid<G> as Grid>::Entity<'a> where Self: 'a;
    type GeometryMap<'a> = <LocalGrid<G> as Grid>::GeometryMap<'a> where Self: 'a;
    type EntityDescriptor = <LocalGrid<G> as Grid>::EntityDescriptor;
    type EntityIter<'a> = <LocalGrid<G> as Grid>::EntityIter<'a> where Self: 'a;

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

    fn geometry_map(
        &self,
        entity_type: Self::EntityDescriptor,
        points: &[Self::T],
    ) -> Self::GeometryMap<'_> {
        self.local_grid.geometry_map(entity_type, points)
    }
}
