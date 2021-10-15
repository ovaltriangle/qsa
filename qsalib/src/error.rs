use std::fmt;

pub type Result<T> = std::result::Result<T, QSAError>;

#[derive(Debug)]
pub enum QSAError {
    BamChecksFailed,
    BAMNotFound,
    DirNotFound,
    CoverageHole,
}

impl fmt::Display for QSAError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            QSAError::BamChecksFailed =>
                write!(f, "Reference sequence length differs between BAM files"),
            QSAError::BAMNotFound =>
                write!(f, "One of the supplied BAM files were not found"),
            QSAError::DirNotFound =>
                write!(f, "One of the supplied directories were not found"),
            QSAError::CoverageHole =>
                write!(f, "One of the supplied BAM files has a coverage hole inside"),
        }
    }
}