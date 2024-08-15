//! Parallel grid
use crate::{
    traits::{Entity, Grid},
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
        entity.map(|e| ParallelGridEntity::new(e, &self.ownership, &self.global_indices))
    }
}

/// Parallel grid
pub struct ParallelGrid<'a, C: Communicator, G: Grid> {
    comm: &'a C,
    serial_grid: G,
    ownership: Vec<Vec<Ownership>>,
    global_indices: Vec<Vec<usize>>,
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
            serial_grid,
            ownership,
            global_indices,
        }
    }
}
impl<'g, C: Communicator, G: Grid> Grid for ParallelGrid<'g, C, G> {
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

    fn entity(&self, dim: usize, local_index: usize) -> Option<Self::Entity<'_>> {
        self.serial_grid
            .entity(dim, local_index)
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
