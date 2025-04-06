//! Input/output operations and utilities.
//! 
//! This module provides basic I/O operations for molecular data:
//! - File format detection
//! - File system operations

use std::io::{self, Read, Seek};
use std::path::Path;
use crate::core::{Result, Error};

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Self {
        Error::Other(Box::new(err))
    }
}

/// A trait that combines Read and Seek capabilities
pub trait ReadSeek: Read + Seek {}
impl<T: Read + Seek> ReadSeek for T {}

/// Trait for detecting file formats
pub trait FormatDetector {
    /// Detect the format of a file by examining its contents
    /// 
    /// # Arguments
    /// 
    /// * `reader` - A reader that supports seeking (to allow peeking at the file contents)
    /// 
    /// # Returns
    /// 
    /// The detected format name as a string (e.g., "molden", "xyz")
    fn detect_format<R: Read + Seek>(&self, reader: &mut R) -> Result<String>;
}

/// Basic file system operations for molecular data
pub struct FileSystem;

impl FileSystem {
    /// Read a file and return a reader
    /// 
    /// # Arguments
    /// 
    /// * `path` - Path to the file to read
    /// 
    /// # Returns
    /// 
    /// A boxed reader that implements Read + Seek
    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<Box<dyn ReadSeek>> {
        let file = std::fs::File::open(path)?;
        Ok(Box::new(file))
    }
    
    /// Create a new file for writing
    /// 
    /// # Arguments
    /// 
    /// * `path` - Path where the file should be created
    /// 
    /// # Returns
    /// 
    /// A boxed writer that implements Write
    pub fn write_file<P: AsRef<Path>>(path: P) -> Result<Box<dyn std::io::Write>> {
        let file = std::fs::File::create(path)?;
        Ok(Box::new(file))
    }
}
