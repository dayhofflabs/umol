//! Plugin system for extending umol's functionality
//! 
//! This module provides a flexible plugin architecture that allows extending umol with:
//! - New molecular properties
//! - New model types
//! - Model-to-model conversions
//! - File format handlers
//! 
//! # Plugin Architecture
//! 
//! The plugin system is built around these core concepts:
//! - [`Plugin`]: The main trait that plugins must implement
//! - [`Registry`]: Central registry for all plugin components
//! - [`Capability`]: Namespaced identifiers for model capabilities
//! 
//! # Example
//! 
//! ```rust
//! use umol::core::{Plugin, Registry, Capability, PluginRequirements};
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
//!                 Capability::new("core", "has_atoms")
//!             ].into_iter().collect(),
//!         }
//!     }
//!     
//!     fn register(&self, registry: &mut Registry) {
//!         // Register your components here
//!         registry.register_capability(
//!             Capability::new("my_plugin", "has_feature")
//!         );
//!     }
//! }
//! ```

use std::collections::{HashMap, HashSet};
use std::sync::{Arc, OnceLock};
use semver::Version;
use crate::core::{Error, Result};

/// A capability identifier with namespace support.
/// 
/// Capabilities are organized in namespaces to avoid conflicts between plugins.
/// The format is `namespace::name`, where:
/// - `namespace`: Typically the plugin name or module (e.g., "quantum", "core")
/// - `name`: The specific capability (e.g., "has_atoms", "has_wavefunction")
/// 
/// # Examples
/// 
/// ```rust
/// use umol::core::Capability;
/// 
/// // Create a capability directly
/// let cap = Capability::new("quantum", "has_wavefunction");
/// assert_eq!(cap.to_string(), "quantum::has_wavefunction");
/// 
/// // Parse from string
/// let cap = Capability::parse("core::has_atoms").unwrap();
/// assert_eq!(cap.namespace(), "core");
/// assert_eq!(cap.name(), "has_atoms");
/// ```
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct Capability {
    namespace: String,
    name: String,
}

impl Capability {
    pub fn new(namespace: impl Into<String>, name: impl Into<String>) -> Self {
        Self {
            namespace: namespace.into(),
            name: name.into(),
        }
    }

    pub fn parse(s: &str) -> Result<Self> {
        let parts: Vec<&str> = s.split("::").collect();
        match parts.as_slice() {
            [namespace, name] => Ok(Self::new(*namespace, *name)),
            _ => Err(Error::InvalidCapability(s.to_string()))
        }
    }

    pub fn namespace(&self) -> &str {
        &self.namespace
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn to_string(&self) -> String {
        format!("{}::{}", self.namespace, self.name)
    }
}

/// Requirements that must be satisfied before a plugin can be registered.
/// 
/// This includes:
/// - Other plugins that must be present with minimum versions
/// - Capabilities that must be available
/// 
/// # Example
/// 
/// ```rust
/// use umol::core::{PluginRequirements, Capability};
/// 
/// let reqs = PluginRequirements {
///     plugins: [("core", "1.0.0")].into_iter()
///         .map(|(k, v)| (k.to_string(), v.parse().unwrap()))
///         .collect(),
///     capabilities: [
///         Capability::new("core", "has_atoms"),
///         Capability::new("core", "has_coordinates")
///     ].into_iter().collect(),
/// };
/// ```
#[derive(Default)]
pub struct PluginRequirements {
    /// Required plugins and their minimum versions
    pub plugins: HashMap<String, Version>,
    /// Required capabilities
    pub capabilities: HashSet<Capability>,
}

/// A lazily initialized component that can be registered by a plugin.
/// 
/// This wrapper ensures that expensive component initialization only happens
/// when the component is first accessed, not when it's registered.
/// 
/// The component is wrapped in an Arc for thread-safe sharing.
pub struct LazyComponent<T> {
    initializer: Box<dyn Fn() -> Result<T> + Send + Sync>,
    value: OnceLock<Arc<T>>,
}

impl<T> LazyComponent<T> {
    pub fn new(initializer: impl Fn() -> Result<T> + Send + Sync + 'static) -> Self {
        Self {
            initializer: Box::new(initializer),
            value: OnceLock::new(),
        }
    }

    pub fn get(&self) -> Result<Arc<T>> {
        self.value
            .get_or_try_init(|| Ok(Arc::new((self.initializer)()?)))
            .cloned()
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
}

/// Central registry for all plugin-provided components.
/// 
/// The registry manages:
/// - Plugin registration and dependencies
/// - Capability tracking
/// - Lazy component initialization
/// - Component access and lookup
/// 
/// # Component Types
/// 
/// - Properties: Calculations that can be performed on models
/// - Models: Different representations of molecular systems
/// - Conversions: Transformations between model types
/// - Formats: File format handlers for IO operations
/// 
/// # Example
/// 
/// ```rust
/// use umol::core::{Registry, Plugin};
/// 
/// let mut registry = Registry::new();
/// 
/// // Register a plugin
/// let my_plugin = MyPlugin;
/// registry.register_plugin(Box::new(my_plugin))?;
/// 
/// // Access components
/// let property = registry.get_property("some_property")?;
/// let model = registry.get_model("some_model")?;
/// ```
pub struct Registry {
    plugins: HashMap<String, (Version, Box<dyn Plugin>)>,
    properties: HashMap<String, LazyComponent<Box<dyn PropertyDefinition>>>,
    models: HashMap<String, LazyComponent<Box<dyn ModelDefinition>>>,
    capabilities: HashSet<Capability>,
    conversions: HashMap<(String, String), LazyComponent<Box<dyn ConversionDefinition>>>,
    formats: HashMap<String, LazyComponent<Box<dyn FormatHandler>>>,
}

impl Registry {
    pub fn new() -> Self {
        Self {
            plugins: HashMap::new(),
            properties: HashMap::new(),
            models: HashMap::new(),
            capabilities: HashSet::new(),
            conversions: HashMap::new(),
            formats: HashMap::new(),
        }
    }

    pub fn register_plugin(&mut self, plugin: Box<dyn Plugin>) -> Result<()> {
        // Validate requirements
        let reqs = plugin.requires();
        
        // Check plugin dependencies
        for (req_name, req_version) in reqs.plugins {
            match self.plugins.get(&req_name) {
                Some((version, _)) if version >= &req_version => (),
                Some((version, _)) => return Err(Error::PluginVersionMismatch {
                    plugin: req_name,
                    required: req_version,
                    found: version.clone(),
                }),
                None => return Err(Error::MissingPlugin(req_name)),
            }
        }
        
        // Check capability dependencies
        for cap in reqs.capabilities {
            if !self.capabilities.contains(&cap) {
                return Err(Error::MissingCapability(cap));
            }
        }
        
        // Register the plugin
        let name = plugin.name().to_string();
        let version = plugin.version();
        
        // Allow plugin to register its components
        plugin.register(self);
        
        self.plugins.insert(name, (version, plugin));
        Ok(())
    }

    pub fn register_capability(&mut self, capability: Capability) {
        self.capabilities.insert(capability);
    }

    pub fn register_property(
        &mut self,
        name: String,
        initializer: impl Fn() -> Result<Box<dyn PropertyDefinition>> + Send + Sync + 'static
    ) {
        self.properties.insert(name, LazyComponent::new(initializer));
    }

    pub fn register_model(
        &mut self,
        name: String,
        initializer: impl Fn() -> Result<Box<dyn ModelDefinition>> + Send + Sync + 'static
    ) {
        self.models.insert(name, LazyComponent::new(initializer));
    }

    pub fn register_conversion(
        &mut self,
        source: String,
        target: String,
        initializer: impl Fn() -> Result<Box<dyn ConversionDefinition>> + Send + Sync + 'static
    ) {
        self.conversions.insert((source, target), LazyComponent::new(initializer));
    }

    pub fn register_format(
        &mut self,
        name: String,
        initializer: impl Fn() -> Result<Box<dyn FormatHandler>> + Send + Sync + 'static
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
        target: &str
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