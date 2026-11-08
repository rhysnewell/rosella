use rosella::data_dir::{remember, remembered};

// One test per binary, since it points CONDA_PREFIX and ROSELLA_DATA_DIR at its own folders.
#[test]
fn a_remembered_directory_is_kept_in_the_active_environment_and_the_variable_wins() {
    let environment = tempfile::tempdir().unwrap();
    let data = tempfile::tempdir().unwrap();
    unsafe {
        std::env::remove_var("ROSELLA_DATA_DIR");
        std::env::set_var("CONDA_PREFIX", environment.path());
    }
    assert_eq!(remembered(), None);
    remember(data.path()).unwrap();
    assert!(environment.path().join("etc/rosella/data_dir").exists());
    assert_eq!(
        remembered(),
        Some(std::fs::canonicalize(data.path()).unwrap())
    );
    unsafe {
        std::env::set_var("ROSELLA_DATA_DIR", "/elsewhere");
    }
    assert_eq!(remembered(), Some("/elsewhere".into()));
}
