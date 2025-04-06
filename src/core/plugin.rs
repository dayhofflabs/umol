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
//! # Example
//!
//! ```rust
//! use umol::core::{Plugin, Registry, Capability, PluginRequirements, ModelProvider, PropertyProvider};
//!
//! struct MyPlugin;
//!
//! impl Plugin for MyPlugin {
//!     fn name(&self) -> &str { "my_plugin" }
//!     fn version(&self) -> Version { "1.0.0".parse().unwrap() }
//!     
//!     fn requires(&self) -> PluginRequirements {
//!         PluginRequirements {
//!             plugins: [("core", "1.0.0")].into_iter()
//!                 .map(|(k, v)| (k.to_string(), v.parse().unwrap()))
//!                 .collect(),
//!             capabilities: [
//!                 Capability::new("core", "has_atoms", 1)
//!             ].into_iter().collect(),
//!         }
//!     }
//!     
//!     fn register(&self, registry: &mut Registry) {
//!         // Register your components here
//!         registry.register_capability(
//!             Capability::new("my_plugin", "has_feature", 1)
//!         );
//!     }
//! }
//!
//! // Indicate that this plugin provides models and properties
//! impl ModelProvider for MyPlugin {}
//! impl PropertyProvider for MyPlugin {}
//! ```

use crate::core::{ConversionMetadata, Entity, Error, Instance, Model, Property, Result};
use semver::Version;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex};

/// A capability that a plugin provides or requires
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct Capability {
    /// The namespace of the capability (usually the plugin name)
    pub namespace: String,
    /// The name of the capability
    pub name: String,
}

impl Capability {
    /// Create a new capability
    pub fn new(namespace: impl Into<String>, name: impl Into<String>) -> Self {
        Self {
            namespace: namespace.into(),
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
#[derive(Default)]
pub struct PluginRequirements {
    /// Required plugins and their versions
    pub plugins: HashMap<String, Version>,
    /// Required capabilities
    pub capabilities: HashSet<Capability>,
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

/// Core trait that all plugins must implement.
///
/// A plugin can provide various components:
/// - Properties for calculating molecular properties
/// - Models for representing molecular systems
/// - Conversions between different model types
/// - File format handlers for IO operations
///
/// # Implementation Notes
///
/// 1. The `name` method should return a unique identifier for your plugin
/// 2. Version numbers should follow semantic versioning
/// 3. Use `requires` to specify dependencies on other plugins
/// 4. Register all components in the `register` method
///
/// # Example
///
/// ```rust
/// use umol::core::{Plugin, Registry, Version};
///
/// struct MyPlugin;
///
/// impl Plugin for MyPlugin {
///     fn name(&self) -> &str { "my_plugin" }
///     fn version(&self) -> Version { "1.0.0".parse().unwrap() }
///     fn register(&self, registry: &mut Registry) {
///         // Register your components here
///     }
/// }
/// ```
pub trait Plugin: Send + Sync {
    /// Unique name of the plugin
    fn name(&self) -> &str;

    /// Semantic version of the plugin
    fn version(&self) -> Version;

    /// Specify plugin requirements
    fn requires(&self) -> PluginRequirements {
        PluginRequirements::default()
    }

    /// Register this plugin's components with the registry
    fn register(&self, registry: &mut Registry);

    /// Check if this plugin provides models
    fn provides_models(&self) -> bool {
        std::any::TypeId::of::<Self>() == std::any::TypeId::of::<dyn ModelProvider>()
    }

    /// Check if this plugin provides properties
    fn provides_properties(&self) -> bool {
        std::any::TypeId::of::<Self>() == std::any::TypeId::of::<dyn PropertyProvider>()
    }

    /// Check if this plugin provides conversions
    fn provides_conversions(&self) -> bool {
        std::any::TypeId::of::<Self>() == std::any::TypeId::of::<dyn ConversionProvider>()
    }

    /// Check if this plugin provides formats
    fn provides_formats(&self) -> bool {
        std::any::TypeId::of::<Self>() == std::any::TypeId::of::<dyn FormatProvider>()
    }

    /// Check if this plugin provides ontology features
    fn provides_ontology(&self) -> bool {
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
    /// Registered capabilities
    capabilities: HashSet<Capability>,
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
            capabilities: HashSet::new(),
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

        // Check capability requirements
        for cap in reqs.capabilities {
            if !self.capabilities.contains(&cap) {
                return false;
            }
        }

        true
    }

    /// Register a plugin
    pub fn register_plugin(&mut self, plugin: Box<dyn Plugin>) -> Result<()> {
        let name = plugin.name().to_string();
        let version = plugin.version();

        // Check requirements
        if !self.satisfies_requirements(plugin.requires()) {
            return Err(Error::MissingDependency(name));
        }

        // Track what the plugin provides
        if plugin.provides_models() {
            self.model_providers.insert(name.clone());
        }
        if plugin.provides_properties() {
            self.property_providers.insert(name.clone());
        }
        if plugin.provides_conversions() {
            self.conversion_providers.insert(name.clone());
        }
        if plugin.provides_formats() {
            self.format_providers.insert(name.clone());
        }
        if plugin.provides_ontology() {
            self.ontology_providers.insert(name.clone());
        }

        // Register the plugin's components
        plugin.register(self);

        // Store the plugin
        self.plugins.insert(name, (version, plugin));

        Ok(())
    }

    pub fn register_capability(&mut self, capability: Capability) {
        self.capabilities.insert(capability);
    }

    pub fn register_property(
        &mut self,
        name: String,
        initializer: impl Fn() -> Result<Box<dyn PropertyDefinition>> + Send + Sync + 'static,
    ) {
        self.properties
            .insert(name, LazyComponent::new(initializer));
    }

    pub fn register_model(
        &mut self,
        name: String,
        initializer: impl Fn() -> Result<Box<dyn ModelDefinition>> + Send + Sync + 'static,
    ) {
        self.models.insert(name, LazyComponent::new(initializer));
    }

    pub fn register_conversion(
        &mut self,
        source: String,
        target: String,
        initializer: impl Fn() -> Result<Box<dyn ConversionDefinition>> + Send + Sync + 'static,
    ) {
        self.conversions
            .insert((source, target), LazyComponent::new(initializer));
    }

    pub fn register_format(
        &mut self,
        name: String,
        initializer: impl Fn() -> Result<Box<dyn FormatHandler>> + Send + Sync + 'static,
    ) {
        self.formats.insert(name, LazyComponent::new(initializer));
    }

    // Accessor methods
    pub fn get_property(&self, name: &str) -> Result<Arc<Box<dyn PropertyDefinition>>> {
        self.properties
            .get(name)
            .ok_or_else(|| Error::PropertyNotFound(name.to_string()))?
            .get()
    }

    pub fn get_model(&self, name: &str) -> Result<Arc<Box<dyn ModelDefinition>>> {
        self.models
            .get(name)
            .ok_or_else(|| Error::ModelNotFound(name.to_string()))?
            .get()
    }

    pub fn get_conversion(
        &self,
        source: &str,
        target: &str,
    ) -> Result<Arc<Box<dyn ConversionDefinition>>> {
        self.conversions
            .get(&(source.to_string(), target.to_string()))
            .ok_or_else(|| Error::ConversionNotFound(source.to_string(), target.to_string()))?
            .get()
    }

    pub fn get_format(&self, name: &str) -> Result<Arc<Box<dyn FormatHandler>>> {
        self.formats
            .get(name)
            .ok_or_else(|| Error::FormatNotFound(name.to_string()))?
            .get()
    }

    pub fn has_capability(&self, capability: &Capability) -> bool {
        self.capabilities.contains(capability)
    }
}
