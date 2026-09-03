//! Render one graph-IR molecule or reaction DSL file as SVG.

use std::error::Error;
use std::ffi::OsString;
use std::path::PathBuf;
use std::{env, fs, io, process};

use umol_graph_ir::dsl::{MoleculeDsl, ReactionDsl};
use umol_io::depict::Depict;

const USAGE: &str = "usage: depict_dsl <molecule|reaction> <input.dsl> <output.svg>";

fn main() {
    if let Err(error) = run() {
        eprintln!("depict_dsl: {error}");
        process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let mut args = env::args_os().skip(1);
    let Some(kind) = args.next() else {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, USAGE).into());
    };
    if kind == "-h" || kind == "--help" {
        println!("{USAGE}");
        return Ok(());
    }
    let input = required_path(args.next())?;
    let output = required_path(args.next())?;
    if args.next().is_some() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, USAGE).into());
    }

    let source = fs::read_to_string(&input).map_err(|error| {
        io::Error::new(
            error.kind(),
            format!("cannot read {}: {error}", input.display()),
        )
    })?;
    let kind = kind.into_string().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "input kind must be valid UTF-8",
        )
    })?;
    let depiction = match kind.as_str() {
        "molecule" => {
            let molecule = source.parse::<MoleculeDsl>()?;
            molecule.molecule().depict()?
        }
        "reaction" => {
            let reaction = source.parse::<ReactionDsl>()?;
            reaction.reaction().depict()?
        }
        _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, USAGE).into()),
    };

    fs::write(&output, depiction.render_svg()).map_err(|error| {
        io::Error::new(
            error.kind(),
            format!("cannot write {}: {error}", output.display()),
        )
    })?;
    Ok(())
}

fn required_path(value: Option<OsString>) -> Result<PathBuf, io::Error> {
    value
        .map(PathBuf::from)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, USAGE))
}
