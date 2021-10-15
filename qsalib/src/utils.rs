use std::fs;
use std::path::PathBuf;

pub fn expand_dir(path: &str, extension: &str) -> Vec<PathBuf> {
    let mut ext_files: Vec<PathBuf> = vec![];

    for entry in fs::read_dir(path).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();

        if extension != "" {
            if let Some(ext) = path.extension() {
                if ext.to_str().unwrap().to_lowercase() != extension.to_string() {
                    continue
                }
            } else {
                continue
            }
        }

        ext_files.push(path);
    }

    ext_files
}