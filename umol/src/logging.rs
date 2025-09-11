//! Logging utilities for the umol project.

use std::env::{set_var, var_os};

use slog::{o, Drain, Level, Logger};
use slog_async::{Async, OverflowStrategy};
use slog_envlogger::EnvLogger;
use slog_term::{FullFormat, TermDecorator};

/// Sets up the global logger with the specified log level
pub fn setup_logger(level: Level) -> Logger {
    let decorator = TermDecorator::new().build();
    let fmt = FullFormat::new(decorator).build();
    if var_os("RUST_LOG").is_none() {
        let lv = match level {
            Level::Critical => "error",
            Level::Error => "error",
            Level::Warning => "warn",
            Level::Info => "info",
            Level::Debug => "debug",
            Level::Trace => "trace",
        };
        set_var("RUST_LOG", lv);
    }
    let filtered = EnvLogger::new(fmt).fuse();
    let drain = Async::new(filtered)
        .chan_size(1024)
        .overflow_strategy(OverflowStrategy::Block)
        .build()
        .fuse();
    Logger::root(drain, o!())
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
