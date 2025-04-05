use crate::core::{Entity, Model};

/// An instance pairs an entity with its representation in a specific model
pub struct Instance<E: Entity, M: Model> {
    entity: E,
    model: M,
}

impl<E: Entity, M: Model> Instance<E, M> {
    pub fn new(entity: E, model: M) -> Self {
        Self { entity, model }
    }

    pub fn entity(&self) -> &E {
        &self.entity
    }

    pub fn model(&self) -> &M {
        &self.model
    }
}

/// Operations connect instances
pub trait Operation<E: Entity, M1: Model, M2: Model> {
    type Error;
    
    fn connect(&self, source: &Instance<E, M1>) -> Result<Instance<E, M2>, Self::Error>;
} 