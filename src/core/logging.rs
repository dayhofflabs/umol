//! Logging utilities for the umol project

use slog::{o, Drain, Logger};

/// Sets up the global logger with the specified log level
pub fn setup_logger(level: slog::Level) -> Logger {
    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::FullFormat::new(decorator)
        .build()
        .filter_level(level)
        .fuse();
    let drain = slog_async::Async::new(drain).build().fuse();
    
    Logger::root(drain, o!("version" => env!("CARGO_PKG_VERSION")))
}
/// Performs any necessary cleanup of the logging system
pub fn teardown_logger() {}

/// Creates a new logger with the specified component name
#[macro_export]
macro_rules! with_logger {
    ($log:expr, $component:expr) => {
        $log.new(o!("component" => $component))
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use slog::info;

    #[test]
    fn test_logger_macro() {
        let root = setup_logger(slog::Level::Info);
        let component_log = with_logger!(root, "test_component");
        
        // This should work and log with the component field
        info!(component_log, "Test message");
    }
}
