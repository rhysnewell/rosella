use std::io::Write;

use anyhow::Result;
use clap::{CommandFactory, crate_version};

use crate::cli::Cli;

/// Read off the command line before clap parses, so a missing required argument cannot
/// stop the manual printing.
pub struct Request {
    pub subcommand: String,
    pub roff_only: bool,
}

pub fn requested() -> Option<Request> {
    let mut arguments = std::env::args().skip(1);
    let subcommand = arguments.next()?;
    let mut request = None;
    for argument in arguments {
        match argument.as_str() {
            "-H" | "--full-help" => {
                request.get_or_insert(false);
            }
            "--full-help-roff" => {
                request = Some(true);
            }
            _ => {}
        }
    }
    request.map(|roff_only| Request {
        subcommand,
        roff_only,
    })
}

pub fn render(subcommand: &str) -> Result<Vec<u8>> {
    let root = Cli::command();
    let command = root
        .get_subcommands()
        .find(|candidate| candidate.get_name() == subcommand)
        .ok_or_else(|| anyhow!("rosella has no {} subcommand", subcommand))?
        .clone()
        .name(String::leak(format!("rosella-{subcommand}")) as &str)
        .version(crate_version!())
        .author(crate::AUTHOR_AND_EMAIL);

    let mut rendered = Vec::new();
    clap_mangen::Man::new(command).render(&mut rendered)?;
    Ok(rendered)
}

pub fn print(request: &Request) -> Result<()> {
    let rendered = render(&request.subcommand)?;
    if request.roff_only {
        std::io::stdout().write_all(&rendered)?;
        return Ok(());
    }

    let mut page = tempfile::NamedTempFile::new()?;
    page.write_all(&rendered)?;
    page.flush()?;
    let shown = std::process::Command::new("man")
        .arg(page.path())
        .status()
        .map(|status| status.success())
        .unwrap_or(false);
    if !shown {
        std::io::stdout().write_all(&rendered)?;
    }
    Ok(())
}
