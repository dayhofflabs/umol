//! Plugin system for extending umol's functionality
//!
//! This module provides a flexible plugin architecture that allows extending umol with:
//! - New molecular properties
//! - New model types
//! - Model-to-model conversions
//! - File format handlers
//! - Ontology definitions
//!
//! # Plugin Architecture
//!
//! The plugin system is built around these core concepts:
//! - [`Plugin`]: The main trait that plugins must implement
//! - [`Registry`]: Central registry for all plugin components
//! - [`Capability`]: Namespaced identifiers for model capabilities
//! - Provider traits: Marker traits indicating what a plugin provides
//!   - [`ModelProvider`]: Provides model implementations
//!   - [`PropertyProvider`]: Provides property calculations
//!   - [`ConversionProvider`]: Provides model conversions
//!   - [`FormatProvider`]: Provides format handlers
//!   - [`OntologyProvider`]: Provides ontology features
//!

use crate::core::{
    error::{
        ModelError, PropertyError, ConversionError, PluginError, FormatError, Result
    },
    ConversionMetadata, Entity, Instance, Model, Capability,
};
use semver::Version;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex};

/// A feature that a plugin provides
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Feature {
    /// The name of the feature
    pub name: String,
}

impl Feature {
    /// Create a new feature
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
        }
    }
}

/// Base trait for property definitions provided by plugins
pub trait PropertyDefinition: Send + Sync {
    /// Get the name of this property
    fn name(&self) -> &str;

    /// Get a description of this property
    fn description(&self) -> &str;

    /// Get the units of this property, if applicable
    fn units(&self) -> Option<&str>;

    /// Get the capabilities required to compute this property
    fn required_capabilities(&self) -> HashSet<Capability>;
}

/// Trait for computing property values
pub trait PropertyCompute<M>: PropertyDefinition
where
    M: Model + 'static,
{
    /// Compute the property for a given instance
    fn compute(&self, instance: &Instance<M>) -> Result<f64>;
}

/// Base trait for model definitions provided by plugins
pub trait ModelDefinition: Send + Sync {
    /// Get the name of this model type
    fn name(&self) -> &str;

    /// Get a description of this model type
    fn description(&self) -> &str;

    /// Get the capabilities provided by this model type
    fn capabilities(&self) -> HashSet<Capability>;

    /// Create a new instance of this model type
    fn create(&self) -> Result<Box<dyn Model<Data = Entity>>>;
}

/// Base trait for conversion definitions provided by plugins
pub trait ConversionDefinition: Send + Sync {
    /// Get the source model type
    fn source_type(&self) -> &str;

    /// Get the target model type
    fn target_type(&self) -> &str;
}

/// Trait for converting between specific model types
pub trait ConversionCompute<M1, M2>: ConversionDefinition
where
    M1: Model + 'static,
    M2: Model + 'static,
{
    /// Convert between model types
    fn convert(&self, source: &Instance<M1>, params: &ConversionMetadata) -> Result<Instance<M2>>;
}

/// Trait for format handlers provided by plugins
pub trait FormatHandler: Send + Sync {
    /// Get the name of this format
    fn name(&self) -> &str;

    /// Get a description of this format
    fn description(&self) -> &str;

    /// Get the file extensions handled by this format
    fn extensions(&self) -> Vec<&str>;

    /// Read a model from a file
    fn read(&self, path: &str) -> Result<Box<dyn Model<Data = Entity>>>;

    /// Write a model to a file
    fn write(&self, model: &dyn Model<Data = Entity>, path: &str) -> Result<()>;
}

/// Base trait for ontology relation definitions provided by plugins
pub trait RelationDefinition: Send + Sync {
    /// Get the name of this relation
    fn name(&self) -> &str;

    /// Get a description of this relation
    fn description(&self) -> &str;

    /// Get the source entity type
    fn source_type(&self) -> &str;

    /// Get the target entity type
    fn target_type(&self) -> &str;
}

/// Trait for checking relations between specific entity types
pub trait RelationCompute: RelationDefinition {
    /// Check if this relation holds between two entities
    fn holds(&self, source: &Entity, target: &Entity) -> bool;
}

/// Marker trait for plugins that provide models
pub trait ModelProvider {}

/// Marker trait for plugins that provide properties
pub trait PropertyProvider {}

/// Marker trait for plugins that provide conversions
pub trait ConversionProvider {}

/// Marker trait for plugins that provide formats
pub trait FormatProvider {}

/// Marker trait for plugins that provide ontology features
pub trait OntologyProvider {}

/// Requirements for a plugin
#[derive(Debug, Default)]
pub struct PluginRequirements {
    /// Required plugins
    pub plugins: Vec<(String, Version)>,
    /// Required features
    pub features: Vec<Feature>,
}

/// A lazy-initialized component that is created on first use
pub struct LazyComponent<T: ?Sized> {
    /// The initialization function
    init: Box<dyn Fn() -> Result<Arc<T>> + Send + Sync>,
    /// The cached value
    value: Mutex<Option<Arc<T>>>,
}

impl<T: ?Sized> LazyComponent<T> {
    /// Create a new lazy component
    pub fn new(init: impl Fn() -> Result<Arc<T>> + Send + Sync + 'static) -> Self {
        Self {
            init: Box::new(init),
            value: Mutex::new(None),
        }
    }

    /// Get the component, initializing it if necessary
    pub fn get(&self) -> Result<Arc<T>> {
        let mut value = self.value.lock().unwrap();
        if let Some(v) = value.as_ref() {
            Ok(v.clone())
        } else {
            let v = (self.init)()?;
            *value = Some(v.clone());
            Ok(v)
        }
    }
}

/// A plugin that provides additional functionality to the framework
pub trait Plugin: Send + Sync + 'static {
    /// Get the name of this plugin
    fn name(&self) -> &str;

    /// Get the version of this plugin
    fn version(&self) -> &str;

    /// Get the features provided by this plugin
    fn features(&self) -> Vec<Feature>;

    /// Get the requirements of this plugin
    fn requires(&self) -> PluginRequirements {
        PluginRequirements::default()
    }

    /// Register this plugin's components
    fn register(&self, _registry: &mut Registry) -> Result<()> {
        todo!()
    }

    /// Check if this plugin provides model functionality
    fn is_model_provider(&self) -> bool {
        std::any::TypeId::of::<Self>() == std::any::TypeId::of::<dyn ModelProvider>()
    }

    /// Check if this plugin provides property functionality
    fn is_property_provider(&self) -> bool {
        std::any::TypeId::of::<Self>() == std::any::TypeId::of::<dyn PropertyProvider>()
    }

    /// Check if this plugin provides conversion functionality
    fn is_conversion_provider(&self) -> bool {
        std::any::TypeId::of::<Self>() == std::any::TypeId::of::<dyn ConversionProvider>()
    }

    /// Check if this plugin provides format functionality
    fn is_format_provider(&self) -> bool {
        std::any::TypeId::of::<Self>() == std::any::TypeId::of::<dyn FormatProvider>()
    }

    /// Check if this plugin provides ontology functionality
    fn is_ontology_provider(&self) -> bool {
        std::any::TypeId::of::<Self>() == std::any::TypeId::of::<dyn OntologyProvider>()
    }
}

/// Registry for managing plugins and their components
pub struct Registry {
    /// Registered plugins
    plugins: HashMap<String, (Version, Box<dyn Plugin>)>,
    /// Plugins that provide models
    model_providers: HashSet<String>,
    /// Plugins that provide properties
    property_providers: HashSet<String>,
    /// Plugins that provide conversions
    conversion_providers: HashSet<String>,
    /// Plugins that provide formats
    format_providers: HashSet<String>,
    /// Plugins that provide ontology features
    ontology_providers: HashSet<String>,
    /// Registered features
    features: HashSet<Feature>,
    /// Registered properties
    properties: HashMap<String, LazyComponent<dyn PropertyDefinition>>,
    /// Registered models
    models: HashMap<String, LazyComponent<dyn ModelDefinition>>,
    /// Registered conversions
    conversions: HashMap<(String, String), LazyComponent<dyn ConversionDefinition>>,
    /// Registered formats
    formats: HashMap<String, LazyComponent<dyn FormatHandler>>,
    /// Registered relations
    relations: HashMap<String, LazyComponent<dyn RelationDefinition>>,
}

impl Registry {
    /// Create a new registry
    pub fn new() -> Self {
        Self {
            plugins: HashMap::new(),
            model_providers: HashSet::new(),
            property_providers: HashSet::new(),
            conversion_providers: HashSet::new(),
            format_providers: HashSet::new(),
            ontology_providers: HashSet::new(),
            features: HashSet::new(),
            properties: HashMap::new(),
            models: HashMap::new(),
            conversions: HashMap::new(),
            formats: HashMap::new(),
            relations: HashMap::new(),
        }
    }

    /// Check if requirements are satisfied
    fn satisfies_requirements(&self, reqs: PluginRequirements) -> bool {
        // Check plugin requirements
        for (name, version) in reqs.plugins {
            if let Some((current, _)) = self.plugins.get(&name) {
                if current < &version {
                    return false;
                }
            } else {
                return false;
            }
        }

        // Check feature requirements
        for feature in reqs.features {
            if !self.features.contains(&feature) {
                return false;
            }
        }

        true
    }

    /// Register a plugin
    pub fn register_plugin(&mut self, plugin: Box<dyn Plugin>) -> Result<()> {
        let name = plugin.name().to_string();
        let version = plugin.version().parse::<Version>().unwrap();

        // Check requirements
        if !self.satisfies_requirements(plugin.requires()) {
            return Err(PluginError::DependencyNotFound(name).into());
        }

        // Track what the plugin provides
        if plugin.is_model_provider() {
            self.model_providers.insert(name.clone());
        }
        if plugin.is_property_provider() {
            self.property_providers.insert(name.clone());
        }
        if plugin.is_conversion_provider() {
            self.conversion_providers.insert(name.clone());
        }
        if plugin.is_format_provider() {
            self.format_providers.insert(name.clone());
        }
        if plugin.is_ontology_provider() {
            self.ontology_providers.insert(name.clone());
        }

        // Register the plugin's components
        plugin.register(self)?;

        // Store the plugin
        self.plugins.insert(name, (version, plugin));

        Ok(())
    }

    /// Register a relation
    pub fn register_relation(
        &mut self,
        name: String,
        initializer: impl Fn() -> Result<Arc<dyn RelationDefinition>> + Send + Sync + 'static,
    ) {
        self.relations.insert(name, LazyComponent::new(initializer));
    }

    /// Get a relation by name
    pub fn get_relation(&self, name: &str) -> Result<Arc<dyn RelationDefinition>> {
        self.relations
            .get(name)
            .ok_or_else(|| PluginError::ComponentInit(format!("Relation not found: {}", name)).into())
            .and_then(|component| component.get())
    }

    pub fn register_feature(&mut self, feature: Feature) {
        self.features.insert(feature);
    }

    pub fn register_property(
        &mut self,
        name: String,
        initializer: impl Fn() -> Result<Arc<dyn PropertyDefinition>> + Send + Sync + 'static,
    ) {
        self.properties.insert(name, LazyComponent::new(initializer));
    }

    pub fn register_model(
        &mut self,
        name: String,
        initializer: impl Fn() -> Result<Arc<dyn ModelDefinition>> + Send + Sync + 'static,
    ) {
        self.models.insert(name, LazyComponent::new(initializer));
    }

    pub fn register_conversion(
        &mut self,
        source: String,
        target: String,
        initializer: impl Fn() -> Result<Arc<dyn ConversionDefinition>> + Send + Sync + 'static,
    ) {
        self.conversions
            .insert((source, target), LazyComponent::new(initializer));
    }

    pub fn register_format(
        &mut self,
        name: String,
        initializer: impl Fn() -> Result<Arc<dyn FormatHandler>> + Send + Sync + 'static,
    ) {
        self.formats.insert(name, LazyComponent::new(initializer));
    }

    // Accessor methods
    pub fn get_property(&self, name: &str) -> Result<Arc<dyn PropertyDefinition>> {
        self.properties
            .get(name)
            .ok_or_else(|| PropertyError::NotFound(name.to_string()))?
            .get()
    }

    pub fn get_model(&self, name: &str) -> Result<Arc<dyn ModelDefinition>> {
        self.models
            .get(name)
            .ok_or_else(|| ModelError::NotFound(name.to_string()))?
            .get()
    }

    pub fn get_conversion(&self, source: &str, target: &str) -> Result<Arc<dyn ConversionDefinition>> {
        self.conversions
            .get(&(source.to_string(), target.to_string()))
            .ok_or_else(|| ConversionError::NotFound(source.to_string(), target.to_string()))?
            .get()
    }

    pub fn get_format(&self, name: &str) -> Result<Arc<dyn FormatHandler>> {
        self.formats
            .get(name)
            .ok_or_else(|| FormatError::NotFound(name.to_string()))?
            .get()
    }

    pub fn has_feature(&self, feature: &Feature) -> bool {
        self.features.contains(feature)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::error::{Error, PluginError, PropertyError};
    use std::sync::Arc;

    #[test]
    fn test_plugin_example() {
        struct MyPlugin;

        impl Plugin for MyPlugin {
            fn name(&self) -> &str { "my_plugin" }
            fn version(&self) -> &str { "1.0.0" }
            
            fn features(&self) -> Vec<Feature> {
                vec![Feature::new("has_feature")]
            }
            
            fn requires(&self) -> PluginRequirements {
                PluginRequirements {
                    plugins: vec![("core".to_string(), "1.0.0".parse().unwrap())],
                    features: vec![Feature::new("has_atoms")],
                }
            }
            
            fn register(&self, registry: &mut Registry) -> Result<()> {
                registry.register_feature(Feature::new("has_feature"));
                Ok(())
            }
        }

        impl ModelProvider for MyPlugin {}
        impl PropertyProvider for MyPlugin {}

        let plugin = MyPlugin;
        let mut registry = Registry::new();
        assert!(plugin.register(&mut registry).is_ok());
        assert!(registry.has_feature(&Feature::new("has_feature")));
    }

    #[test]
    fn test_plugin_version_compatibility() {
        struct TestPlugin;

        impl Plugin for TestPlugin {
            fn name(&self) -> &str { "test_plugin" }
            fn version(&self) -> &str { "2.0.0" }
            fn features(&self) -> Vec<Feature> { vec![] }
            fn requires(&self) -> PluginRequirements {
                PluginRequirements {
                    plugins: vec![("core".to_string(), "1.0.0".parse().unwrap())],
                    features: vec![],
                }
            }
            fn register(&self, _: &mut Registry) -> Result<()> { Ok(()) }
        }

        let mut registry = Registry::new();
        
        // Register core plugin first
        struct CorePlugin;
        impl Plugin for CorePlugin {
            fn name(&self) -> &str { "core" }
            fn version(&self) -> &str { "0.9.0" }  // Lower version than required
            fn features(&self) -> Vec<Feature> { vec![] }
            fn register(&self, _: &mut Registry) -> Result<()> { Ok(()) }
        }

        // Register core plugin
        registry.register_plugin(Box::new(CorePlugin)).unwrap();

        // Try to register test plugin - should fail due to version mismatch
        let result = registry.register_plugin(Box::new(TestPlugin));
        assert!(result.is_err());
        if let Err(Error::Plugin(PluginError::DependencyNotFound(_))) = result {
            // Expected error
        } else {
            panic!("Unexpected error: {:?}", result);
        }
    }

    #[test]
    fn test_plugin_feature_dependencies() {
        struct TestPlugin;

        impl Plugin for TestPlugin {
            fn name(&self) -> &str { "test_plugin" }
            fn version(&self) -> &str { "1.0.0" }
            fn features(&self) -> Vec<Feature> { vec![] }
            fn requires(&self) -> PluginRequirements {
                PluginRequirements {
                    plugins: vec![],
                    features: vec![Feature::new("required_feature")],
                }
            }
            fn register(&self, _: &mut Registry) -> Result<()> { Ok(()) }
        }

        let mut registry = Registry::new();
        
        // Try to register plugin without required feature
        let result = registry.register_plugin(Box::new(TestPlugin));
        assert!(result.is_err());
        if let Err(Error::Plugin(PluginError::DependencyNotFound(_))) = result {
            // Expected error
        } else {
            panic!("Unexpected error: {:?}", result);
        }

        // Add required feature
        registry.register_feature(Feature::new("required_feature"));

        // Now registration should succeed
        assert!(registry.register_plugin(Box::new(TestPlugin)).is_ok());
    }

    #[test]
    fn test_lazy_component_initialization() {
        let mut registry = Registry::new();
        
        // Register a property that fails to initialize
        registry.register_property(
            "failing_prop".to_string(),
            || -> Result<Arc<dyn PropertyDefinition>> {
                Err(Error::Property(PropertyError::NotFound("test".to_string())))
            }
        );

        // Try to get the property - should fail
        let result = registry.get_property("failing_prop");
        assert!(result.is_err());
        if let Err(Error::Property(PropertyError::NotFound(_))) = result {
            // Expected error
        } else {
            panic!("Expected PropertyError::NotFound, got a different error");
        }

        // Register a successful property
        struct TestProperty;
        impl PropertyDefinition for TestProperty {
            fn name(&self) -> &str { "test" }
            fn description(&self) -> &str { "test" }
            fn units(&self) -> Option<&str> { None }
            fn required_capabilities(&self) -> HashSet<Capability> { HashSet::new() }
        }

        registry.register_property(
            "success_prop".to_string(),
            || Ok(Arc::new(TestProperty) as Arc<dyn PropertyDefinition>)
        );

        // First access should initialize
        let prop1 = registry.get_property("success_prop");
        assert!(prop1.is_ok());

        // Second access should return cached value
        let prop2 = registry.get_property("success_prop");
        assert!(prop2.is_ok());
    }
}
