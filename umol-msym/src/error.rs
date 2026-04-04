use std::ffi::CStr;
use umol_msym_sys as ffi;

#[derive(Debug, thiserror::Error)]
#[error("{message}")]
pub struct Error {
    pub code: i32,
    pub message: String,
}

impl Error {
    pub(crate) fn from_code(code: ffi::msym_error_t) -> Self {
        let message = unsafe {
            let ptr = ffi::msymErrorString(code);
            if ptr.is_null() {
                "unknown error".to_string()
            } else {
                CStr::from_ptr(ptr).to_string_lossy().into_owned()
            }
        };
        let detail = unsafe {
            let ptr = ffi::msymGetErrorDetails();
            if ptr.is_null() {
                None
            } else {
                let s = CStr::from_ptr(ptr).to_string_lossy();
                if s.is_empty() { None } else { Some(s.into_owned()) }
            }
        };
        let message = match detail {
            Some(d) => format!("{message}: {d}"),
            None => message,
        };
        Self { code, message }
    }
}

pub(crate) fn check(code: ffi::msym_error_t) -> Result<(), Error> {
    if code == ffi::MSYM_SUCCESS {
        Ok(())
    } else {
        Err(Error::from_code(code))
    }
}
