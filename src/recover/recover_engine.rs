use anyhow::Result;


pub fn run_recover(m: &clap::ArgMatches) -> Result<()> {
    let mut recover_engine = RecoverEngine::new(m)?;
    recover_engine.run()?;
    Ok(())
}

struct RecoverEngine {

}

impl RecoverEngine {
    pub fn new(m: &clap::ArgMatches) -> Result<Self> {
        // create output directory
        let output_directory = m.get_one::<String>("output-directory").unwrap();
        // create the output directory but do not fail if it already exists
        std::fs::create_dir_all(&output_directory)?;

        Ok(
            Self {
                
            }
        )
    }

    pub fn run(&mut self) -> Result<()> {
        unimplemented!();
    }
}