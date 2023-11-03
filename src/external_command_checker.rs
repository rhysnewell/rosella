use std;
use std::io::Read;

use anyhow::Result;

use log::{debug, error};

pub fn check_for_bwa() -> Result<()> {
    self::check_for_external_command_presence("BWA", "which bwa")
}

pub fn check_for_flight() -> Result<()>  {
    self::check_for_external_command_presence("flight", "flight bin -h")
}

pub fn check_for_minimap2() -> Result<()>  {
    self::check_for_external_command_presence("minimap2", "which minimap2")
}

pub fn check_for_ngmlr() -> Result<()>  {
    self::check_for_external_command_presence("ngmlr", "which ngmlr")
}


pub fn check_for_samtools() -> Result<()>  {
    self::check_for_external_command_presence("samtools", "which samtools")
}

fn check_for_external_command_presence(executable_name: &str, testing_cmd: &str) -> Result<()> {
    debug!("Checking for {} ..", executable_name);
    let mut cmd = std::process::Command::new("bash");
    cmd.arg("-c")
        .arg(testing_cmd)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    let mut process = cmd.spawn().expect("Unable to execute bash");
    let es = process.wait().expect(&format!(
        "Failed to glean exitstatus while checking for presence of {}",
        executable_name
    ));
    if !es.success() {
        error!(
            "Could not find an available {} executable.",
            executable_name
        );
        let mut err = String::new();
        process
            .stderr
            .expect("Failed to grab stderr from failed executable finding process")
            .read_to_string(&mut err)
            .expect("Failed to read stderr into string");
        error!("The STDERR was: {:?}", err);
        bail!(
            "Cannot continue without {}. Testing for presence with `{}` failed",
            executable_name, testing_cmd
        );
    }

    Ok(())
}