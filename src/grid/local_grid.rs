//! Serial grids
#[cfg(feature = "serde")]
use crate::traits::ConvertToSerializable;
use crate::{
    traits::{Entity, Grid},
    types::Ownership,
};
#[cfg(feature = "serde")]
use itertools::izip;
#[cfg(feature = "serde")]
use std::hash::Hash;
use std::{collections::HashMap, fmt::Debug};

mod single_element;
pub use single_element::{SingleElementGrid, SingleElementGridBuilder};

mod mixed;
pub use mixed::{MixedGrid, MixedGridBuilder};

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
    ownership: HashMap<G::EntityDescriptor, Vec<Ownership>>,
    global_indices: HashMap<G::EntityDescriptor, Vec<usize>>,
}

#[cfg(feature = "serde")]
#[derive(serde::Serialize, Debug, serde::Deserialize)]
#[serde(bound = "for<'de2> S: serde::Deserialize<'de2>")]
/// A serde serializable grid
pub struct SerializableLocalGrid<EntityDescriptor: serde::Serialize, S: serde::Serialize>
where
    for<'de2> S: serde::Deserialize<'de2>,
    for<'de2> EntityDescriptor: serde::Deserialize<'de2>,
{
    serial_grid: S,
    ownership_keys: Vec<EntityDescriptor>,
    ownership_values: Vec<Vec<Ownership>>,
    global_indices_keys: Vec<EntityDescriptor>,
    global_indices_values: Vec<Vec<usize>>,
}

#[cfg(feature = "serde")]
impl<
    S: serde::Serialize,
    EntityDescriptor: Debug + PartialEq + Eq + Clone + Copy + Hash + serde::Serialize,
    G: Grid<EntityDescriptor = EntityDescriptor> + Sync + ConvertToSerializable<SerializableType = S>,
> ConvertToSerializable for LocalGrid<G>
where
    for<'de2> S: serde::Deserialize<'de2>,
    for<'de2> EntityDescriptor: serde::Deserialize<'de2>,
{
    type SerializableType = SerializableLocalGrid<G::EntityDescriptor, S>;
    fn to_serializable(&self) -> SerializableLocalGrid<G::EntityDescriptor, S> {
        let mut ownership_keys = vec![];
        let mut ownership_values = vec![];
        for (i, j) in &self.ownership {
            ownership_keys.push(*i);
            ownership_values.push(j.clone());
        }
        let mut global_indices_keys = vec![];
        let mut global_indices_values = vec![];
        for (i, j) in &self.global_indices {
            global_indices_keys.push(*i);
            global_indices_values.push(j.clone());
        }
        SerializableLocalGrid {
            serial_grid: self.serial_grid.to_serializable(),
            ownership_keys,
            ownership_values,
            global_indices_keys,
            global_indices_values,
        }
    }
    fn from_serializable(s: SerializableLocalGrid<G::EntityDescriptor, S>) -> Self {
        Self {
            serial_grid: G::from_serializable(s.serial_grid),
            ownership: {
                let mut m = HashMap::new();
                for (i, j) in izip!(s.ownership_keys, s.ownership_values) {
                    m.insert(i, j);
                }
                m
            },
            global_indices: {
                let mut m = HashMap::new();
                for (i, j) in izip!(s.global_indices_keys, s.global_indices_values) {
                    m.insert(i, j);
                }
                m
            },
        }
    }
}

impl<G: Grid + Sync> LocalGrid<G> {
    /// Create new
    pub fn new(
        serial_grid: G,
        ownership: HashMap<G::EntityDescriptor, Vec<Ownership>>,
        global_indices: HashMap<G::EntityDescriptor, Vec<usize>>,
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

    fn entity(
        &self,
        entity_type: Self::EntityDescriptor,
        serial_index: usize,
    ) -> Option<Self::Entity<'_>> {
        self.serial_grid.entity(entity_type, serial_index).map(|e| {
            GridEntity::new(
                e,
                &self.ownership[&entity_type],
                &self.global_indices[&entity_type],
            )
        })
    }

    fn entity_types(&self, dim: usize) -> &[Self::EntityDescriptor] {
        self.serial_grid.entity_types(dim)
    }

    fn entity_count(&self, entity_type: Self::EntityDescriptor) -> usize {
        self.serial_grid.entity_count(entity_type)
    }

    fn entity_iter(&self, entity_type: Self::EntityDescriptor) -> Self::EntityIter<'_> {
        GridEntityIter::new(
            self.serial_grid.entity_iter(entity_type),
            &self.ownership[&entity_type],
            &self.global_indices[&entity_type],
        )
    }

    fn entity_from_id(
        &self,
        entity_type: Self::EntityDescriptor,
        id: usize,
    ) -> Option<Self::Entity<'_>> {
        self.serial_grid.entity_from_id(entity_type, id).map(|e| {
            GridEntity::new(
                e,
                &self.ownership[&entity_type],
                &self.global_indices[&entity_type],
            )
        })
    }

    fn geometry_map(
        &self,
        entity_type: Self::EntityDescriptor,
        geometry_degree: usize,
        points: &[Self::T],
    ) -> Self::GeometryMap<'_> {
        self.serial_grid
            .geometry_map(entity_type, geometry_degree, points)
    }
}
