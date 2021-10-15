// #[warn(missing_docs)]  // enable this when checking everything got documented

/// Matrices useful to the study of quasispecies viruses, like efficiency (normalized entropy).
pub mod matrices;

/// Entry point for quasispecies analysis starting with BAM files.
pub mod bam;

/// Functions which do not fall in a specific category and can be used wherever in the crate.
mod utils;

/// Custom error definitions for the `qsalib` crate.
pub mod error;

/// `qsalib` prelude, useful to explore the library without having to import everything manually.
pub mod prelude {
    pub use crate::bam::{BamDataBuilder, BamData, Bam};
    pub use crate::matrices::Matrices;
    pub use crate::error::{Result, QSAError};
}