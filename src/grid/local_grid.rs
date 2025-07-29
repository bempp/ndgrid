//! Serial grids
#[cfg(feature = "serde")]
use crate::traits::ConvertToSerializable;
use crate::{
    traits::{Entity, Grid},
    types::Ownership,
};

mod single_element;
pub use single_element::{SingleElementGrid, SingleElementGridBuilder};

/// Grid entity
#[derive(Debug)]
pub struct GridEntity<E: Entity> {
    serial_entity: E,
    ownership: Ownership,
    global_index: usize,
}

impl<E: Entity> GridEntity<E> {
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
impl<E: Entity> Entity for GridEntity<E> {
    type T = E::T;
    type EntityDescriptor = E::EntityDescriptor;
    type Topology<'a>
        = E::Topology<'a>
    where
        Self: 'a;
    type Geometry<'a>
        = E::Geometry<'a>
    where
        Self: 'a;
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
    fn id(&self) -> Option<usize> {
        self.serial_entity.id()
    }
}

/// Grid entity iterator
#[derive(Debug)]
pub struct GridEntityIter<'a, E: Entity, EntityIter: Iterator<Item = E>> {
    iter: EntityIter,
    ownership: &'a [Ownership],
    global_indices: &'a [usize],
}

impl<'a, E: Entity, EntityIter: Iterator<Item = E>> GridEntityIter<'a, E, EntityIter> {
    /// Create new
    pub fn new(iter: EntityIter, ownership: &'a [Ownership], global_indices: &'a [usize]) -> Self {
        Self {
            iter,
            ownership,
            global_indices,
        }
    }
}
impl<E: Entity, EntityIter: Iterator<Item = E>> Iterator for GridEntityIter<'_, E, EntityIter> {
    type Item = GridEntity<E>;

    fn next(&mut self) -> Option<GridEntity<E>> {
        let entity = self.iter.next();
        entity.map(|e| GridEntity::new(e, self.ownership, self.global_indices))
    }
}

/// Local grid on a process
#[derive(Debug)]
pub struct LocalGrid<G: Grid + Sync> {
    serial_grid: G,
    ownership: Vec<Vec<Ownership>>,
    global_indices: Vec<Vec<usize>>,
}

#[cfg(feature = "serde")]
#[derive(serde::Serialize, Debug, serde::Deserialize)]
#[serde(bound = "for<'de2> S: serde::Deserialize<'de2>")]
/// A serde serializable grid
pub struct SerializableLocalGrid<S: serde::Serialize>
where
    for<'de2> S: serde::Deserialize<'de2>,
{
    serial_grid: S,
    ownership: Vec<Vec<Ownership>>,
    global_indices: Vec<Vec<usize>>,
}

#[cfg(feature = "serde")]
impl<S: serde::Serialize, G: Grid + Sync + ConvertToSerializable<SerializableType = S>>
    ConvertToSerializable for LocalGrid<G>
where
    for<'de2> S: serde::Deserialize<'de2>,
{
    type SerializableType = SerializableLocalGrid<S>;
    fn to_serializable(&self) -> SerializableLocalGrid<S> {
        SerializableLocalGrid {
            serial_grid: self.serial_grid.to_serializable(),
            ownership: self.ownership.clone(),
            global_indices: self.global_indices.clone(),
        }
    }
    fn from_serializable(s: SerializableLocalGrid<S>) -> Self {
        Self {
            serial_grid: G::from_serializable(s.serial_grid),
            ownership: s.ownership,
            global_indices: s.global_indices,
        }
    }
}

impl<G: Grid + Sync> LocalGrid<G> {
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
impl<G: Grid + Sync> Grid for LocalGrid<G> {
    type T = G::T;
    type Entity<'a>
        = GridEntity<G::Entity<'a>>
    where
        Self: 'a;
    type GeometryMap<'a>
        = G::GeometryMap<'a>
    where
        Self: 'a;
    type EntityDescriptor = G::EntityDescriptor;
    type EntityIter<'a>
        = GridEntityIter<'a, G::Entity<'a>, G::EntityIter<'a>>
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
            .map(|e| GridEntity::new(e, &self.ownership[dim], &self.global_indices[dim]))
    }

    fn entity_types(&self, dim: usize) -> &[Self::EntityDescriptor] {
        self.serial_grid.entity_types(dim)
    }

    fn entity_count(&self, entity_type: Self::EntityDescriptor) -> usize {
        self.serial_grid.entity_count(entity_type)
    }

    fn entity_iter(&self, dim: usize) -> Self::EntityIter<'_> {
        GridEntityIter::new(
            self.serial_grid.entity_iter(dim),
            &self.ownership[dim],
            &self.global_indices[dim],
        )
    }

    fn entity_from_id(&self, dim: usize, id: usize) -> Option<Self::Entity<'_>> {
        self.serial_grid
            .entity_from_id(dim, id)
            .map(|e| GridEntity::new(e, &self.ownership[dim], &self.global_indices[dim]))
    }

    fn geometry_map(
        &self,
        entity_type: Self::EntityDescriptor,
        points: &[Self::T],
    ) -> Self::GeometryMap<'_> {
        self.serial_grid.geometry_map(entity_type, points)
    }
}
